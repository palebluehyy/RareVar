#!/usr/bin/env python

import os
import sys
import pysam
import multiprocessing as mp

def baseq_thres_to_char(thres,encoding):
	thres_char = ','
	if encoding == 'S' or encoding == 'L':
		thres_char = str(unichr(thres+33))
	if encoding == 'X':
		thres_char = str(unichr(thres+64+5))
	if encoding == 'I':
		thres_char = str(unichr(thres+64))
	if encoding == 'J':
		thres_char = str(unichr(thres+64-3))
	return thres_char
	
def strand_wdis(bam_fname, chrom, loc, alt,baseq_char_thres, max_depth=100000):
	l_dis = 0 # distance to leftmost aligned position
	r_dis = 0 # distance to rightmost aligned position
	fwd_wdis = 0.00
	rev_wdis = 0.00
	re = ""
	samfile = pysam.Samfile(bam_fname,"rb")
	for pileupcolumn in samfile.pileup(chrom,loc-1,loc,max_depth=max_depth):
		if pileupcolumn.pos == loc-1:
			fwd_dict ={}
			rev_dict ={}
			for pileupread in pileupcolumn.pileups:
				if pileupread.alignment.seq[pileupread.qpos] == alt:
					if pileupread.alignment.qual[pileupread.qpos] > baseq_char_thres:
						l_dis = loc - pileupread.alignment.pos
						r_dis = pileupread.alignment.aend - loc
						if l_dis < r_dis:
							dis_val = l_dis
						else:
							dis_val = r_dis
						if pileupread.alignment.is_reverse:
							if dis_val in rev_dict.keys():
								rev_dict[dis_val] = rev_dict[dis_val] + 1
							else:
								rev_dict[dis_val] = 1
						else:
							if dis_val in fwd_dict.keys():
								fwd_dict[dis_val] = fwd_dict[dis_val] + 1
							else:
								fwd_dict[dis_val] = 1
	fwd_sum = sum(fwd_dict.values()) * 1.00
	rev_sum = sum(rev_dict.values()) * 1.00
	#sys.stdout.write('fwd_sum: ' + str(fwd_sum) + ' rev_sum: ' + str(rev_sum) + '\n')
	for key in fwd_dict.keys():
		fwd_wdis = fwd_wdis + fwd_dict[key]/fwd_sum * key
	for key in rev_dict.keys():
		rev_wdis = rev_wdis + rev_dict[key]/rev_sum * key
	#sys.stdout.write('fwd_wdis: ' + str(fwd_wdis) + ' rev_wdis: ' + str(rev_wdis) + '\n')
	re=('\t'.join(["{0:.2f}".format(round(fwd_wdis,2)), "{0:.2f}".format(round(rev_wdis,2))]))
	#sys.stdout.write('re str: ' + re + '\n')
	return re
	
def worker(var_file_nm, out_file_nm, pid, nproc, bam_fname, baseq_char_thres):
	var_file = open(var_file_nm)
	out_file_nm = out_file_nm + '.' + str(pid)
	out_file = open(out_file_nm, "w")
	ln = 0
	line = var_file.next()
	for line in var_file:
		if ln % nproc == pid:
			line = line.rstrip("\n")
			cols = []
			cols = line.split("\t") # Input file should contain "chrom pos_srt pos_end".
			chrom = cols[0]
			loc = int(cols[1])
			alt = cols[3]
			re = strand_wdis(bam_fname, chrom, loc, alt,baseq_char_thres)
			out_file.write(str(ln)+'\t'+ line + '\t' + re +'\n')
		ln = ln + 1
	var_file.close()
	out_file.close()
									
def main():
	opt = sys.argv[1:] # input_fname[0] is the bam file, input_fname[1] is the variant file.

	if len(opt) < 7:
		print("""Usage: ./multiproc_strand_wdis_read.py -p # var.bam var.txt baseq_thres base_qual_encoding out.txt
-p is the number of processes used

var.bam file should be indexed by samtools.

var.txt should contain 4 columns: chrom pos(1-based) var_ref var_alt

baseq_thres should be an integer 

base_qual_encoding can take one of the following values: 
				S - Sanger        Phred+33,  raw reads typically (0, 40)
				X - Solexa        Solexa+64, raw reads typically (-5, 40)
				I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
				J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
				with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
				(Note: See discussion above).
				L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
				
out.txt: output file name

The output will attach 2 columns to the original region file: fwd_wdis rev_wdis""")
		sys.exit(1)

	nproc = int(opt[1])
	bam_fname = opt[2]
	var_file_nm = opt[3]
	var_file = open(opt[3])
	baseq_thres = int(opt[4])
	base_qual_encoding = opt[5]
	out_file_nm = opt[6]
	baseq_char_thres = baseq_thres_to_char(baseq_thres,base_qual_encoding)
	out_file_fh = open(out_file_nm, "w")
	out_file_fh.write(var_file.readline().rstrip("\n")+"\tfwd_wdis\trev_wdis\n")
	var_file.close()
	out_file_fh.close()
	#Setup a list of processes 
	processes = [mp.Process(target=worker, args=(var_file_nm, out_file_nm, x, nproc, bam_fname, baseq_char_thres)) for x in range(nproc)]
	# Run processes
	for p in processes:
		p.start()
		# Exit the completed processes
	for p in processes:
			p.join()
	os.system('cat ' + out_file_nm + '.* |sort -k1n |cut -f2- >> ' + out_file_nm )
	os.system('rm ' + out_file_nm + '.*')

if __name__ == "__main__":
	main()

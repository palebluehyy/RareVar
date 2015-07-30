#!/usr/bin/env python

import sys
import os
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
	
def MAPQ_altavg_diffratio(bam_fname, chrom, loc, ref, alt, baseq_char_thres, max_depth=100000):
	alt_avg = 0.00
	ref_avg = 0.00
	re = ""
	samfile = pysam.Samfile(bam_fname,"rb")
	for pileupcolumn in samfile.pileup(chrom,loc-1,loc,max_depth=max_depth):
		if pileupcolumn.pos == loc-1:
			alt_dict ={}
			ref_dict ={}
			for pileupread in pileupcolumn.pileups:
				if pileupread.alignment.seq[pileupread.qpos] == alt:
					if pileupread.alignment.qual[pileupread.qpos] > baseq_char_thres:
						if pileupread.alignment.mapq in alt_dict.keys():
							alt_dict[pileupread.alignment.mapq] = alt_dict[pileupread.alignment.mapq] + 1
						else:
							alt_dict[pileupread.alignment.mapq] = 1
				elif pileupread.alignment.seq[pileupread.qpos] == ref:
					if pileupread.alignment.qual[pileupread.qpos] > baseq_char_thres:
						if pileupread.alignment.mapq in ref_dict.keys():
							ref_dict[pileupread.alignment.mapq] = ref_dict[pileupread.alignment.mapq] + 1
						else:
							ref_dict[pileupread.alignment.mapq] = 1
	alt_sum = sum(alt_dict.values()) * 1.00
	ref_sum = sum(ref_dict.values()) * 1.00
	#sys.stdout.write('fwd_sum: ' + str(fwd_sum) + ' rev_sum: ' + str(rev_sum) + '\n')
	for key in alt_dict.keys():
		alt_avg = alt_avg + alt_dict[key]/alt_sum * key
	for key in ref_dict.keys():
		ref_avg = ref_avg + ref_dict[key]/ref_sum * key
	alt_refdiff = (alt_avg - ref_avg)
	#sys.stdout.write('alt_avg: ' + str(alt_avg) + ' alt_refdiff_ratio: ' + str(alt_refdiff_ratio) + '\n')
	re=('\t'.join(["{0:.2f}".format(round(alt_avg,2)),"{0:.2f}".format(round(ref_avg,2)), "{0:.2f}".format(round(alt_refdiff,2))]))
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
			ref = cols[2]
			alt = cols[3]
			re = MAPQ_altavg_diffratio(bam_fname, chrom, loc, ref, alt,baseq_char_thres)
			out_file.write(str(ln)+'\t'+ line + '\t' + re +'\n')
		ln = ln + 1
	var_file.close()
	out_file.close()					
									
def main():
	opt = sys.argv[1:] # input_fname[0] is the bam file, input_fname[1] is the variant file.

	if len(opt) < 7:
		print("""Usage: ./MAPQ_by_alt_avg_diffratio.py -p # var.bam var.txt baseq_thres base_qual_encoding out.txt

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

The output will attach 3 columns to the original region file: avg_alt alt_refdiff alt_refdiff_ratio
avg_alt: average alternative reads MAPQ from both strands
avg_ref: average reference reads MAPQ from both strands
alt_refdiff: alt_avg - ref_avg
""")
		sys.exit(1)
	nproc = int(opt[1])
	bam_fname = opt[2]
	var_file_nm = opt[3]
	var_file = open(var_file_nm)
	baseq_thres = int(opt[4])
	base_qual_encoding = opt[5]
	out_file_nm = opt[6]
	out_file_fh = open(out_file_nm, "w")
	baseq_char_thres = baseq_thres_to_char(baseq_thres,base_qual_encoding)
	out_file_fh.write(var_file.readline().rstrip("\n")+"\tavg_alt\tavg_ref\talt_refdiff\n")
	out_file_fh.close()
	var_file.close()
    
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

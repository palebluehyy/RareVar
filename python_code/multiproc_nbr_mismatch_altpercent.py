#!/usr/bin/env python

import sys
import os
import pysam
import re
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

def baseq_char_ASCII(letter,encoding):
	val = 0
	if encoding == 'S' or encoding == 'L':
		val = ord(letter) - 33
	if encoding == 'X' or encoding == 'I' or encoding == 'J':
		val = ord(letter) - 64
	return val	

def count_mismatch_inwindow(read, rpos, window_size):
	win_srt = rpos - window_size
	win_end = rpos + window_size
	# count number of insertions before win_srt and win_end
	nins_wsrt = count_nins(read.alignment.cigar, rpos)
	nins_wend = count_nins(read.alignment.cigar, rpos)
	for i in range(len(read.alignment.tags)):
		if read.alignment.tags[i][0] == 'MD':
			md_str = read.alignment.tags[i][1]
	md = parse_md(md_str)
	nsc = 0 ## number of softclipping bases
	if read.alignment.cigar[0][0] == 4: ## start with a softclipping
		nsc = read.alignment.cigar[0][1]
	srt = win_srt - nins_wsrt - nsc
	if srt < 0:
		srt = 0
	end = win_end - nins_wend - nsc
	nmm = 0 ## number of mismatches in the window
	nmm = nmismatch_perread(md, srt, end)
	return nmm

def count_nins(cigar, loc):
	tpos = 0
	i = 0
	nins = 0
	while i < len(cigar) and tpos + cigar[i][1] < loc:
		if cigar[i][0] == 1: # insertions
			nins = nins +1
		if cigar[i][0] != 2: # deletions
			tpos = tpos + cigar[i][1]
		i = i + 1
	return nins

def parse_md(md_str):
	splitter = re.compile("\\^[ACGT]+")
	nodel_str = splitter.split(md_str)
	md = [re.findall("([0-9]+)|([ATGC]+)", s) for s in nodel_str]
	md = sum(md, [])
	return md

def nmismatch_perread(md, srt, end):
	## md is a list of tuples, each tuple follows the pattern('num', 'str')), deletions are already deleted
	loc = srt
	tpos = 0
	i = 0
	nmm = 0
	gtsrt_tag = 0 
	while i < len(md) and tpos <= loc:
		if md[i][0] != '': # has numbers
			tpos = tpos + int(md[i][0])
		else:
			tpos = tpos + len(md[i][1])
			if gtsrt_tag == 1:
				nmm = nmm +1
		if tpos >= loc:
			loc = end
			if md[i][0] == '':
				nmm = nmm + 1
			gtsrt_tag = 1
		i = i + 1
	return nmm

def nbr_mismatch_altpercent(bam_fname, chrom, loc, ref, alt, baseq_char_thres,encoding, window_size, max_depth):
	nread_alt = 0
	nread_mismatch_alt = 0
	nmm = 0
	mm_percent = 0.00
	ans = ""
	samfile = pysam.Samfile(bam_fname,"rb")
	for pileupcolumn in samfile.pileup(chrom,loc-1,loc,max_depth=max_depth):
		if pileupcolumn.pos == loc-1:
			for pileupread in pileupcolumn.pileups:
				if pileupread.is_del != 1 and pileupread.alignment.seq[pileupread.qpos] == alt: ## note: is this place is a deletion, then the base is not correct
					if pileupread.alignment.qual[pileupread.qpos] >= baseq_char_thres:
						nread_alt = nread_alt + 1
						tnmm = count_mismatch_inwindow(pileupread, pileupread.qpos, window_size) - 1
						#print pileupread.alignment.qname + '\t number of mismatches' + str(tnmm)
						if tnmm >0:
							nread_mismatch_alt = nread_mismatch_alt + 1
						nmm = nmm + tnmm
	mm_percent = float(nread_mismatch_alt)/nread_alt
	if nread_mismatch_alt > 0:
		avg_mm_perread = float(nmm)/nread_mismatch_alt
	else:
		avg_mm_perread = 0
	ans=('\t'.join([str(nread_alt), str(nread_mismatch_alt), "{0:.2f}".format(round(mm_percent,2)),"{0:.2f}".format(round(avg_mm_perread,2))]))
	return ans

def worker(var_file_nm, out_file_nm, pid, nproc, bam_fname, baseq_char_thres,base_qual_encoding, window_size, max_depth):
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
			re = nbr_mismatch_altpercent(bam_fname, chrom, loc, ref, alt,baseq_char_thres, base_qual_encoding, window_size, max_depth)
			out_file.write(str(ln)+'\t'+ line + '\t' + re +'\n')
		ln = ln + 1
	var_file.close()
	out_file.close()

def main():
	opt = sys.argv[1:] # input_fname[0] is the bam file, input_fname[1] is the variant file.

	if len(opt) < 9:
		print("""Usage: ./nbr_mismatch_altpercent.py -p # var.bam var.txt baseq_thres base_qual_encoding window_size max_depth out.txt

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
window_size should be an integer, determines the upstream and downstream context to look at

max_depth is the maximum depth to look at

The output will attach 4 columns to the original region file: nread_alt nread_mismatch_alt mm_percent
nread_alt: number of alternative reads from both strands
nread_mismatch_alt: number of reads containing mismatches within the window_size from both strands
mm_percent: nread_mismatch_alt/ nread_alt
avg_mm_perread: sum(number of mismatches in each alternative read other than the alternative base mismatch) / nread_mismatch_alt
""")
		sys.exit(1)
	
	nproc = int(opt[1])
	bam_fname = opt[2]
	var_file_nm = opt[3]
	var_file = open(var_file_nm)
	baseq_thres = int(opt[4])
	base_qual_encoding = opt[5]
	window_size = int(opt[6])
	max_depth = int(opt[7])
	out_file_nm = opt[8]
	out_file_fh = open(out_file_nm, "w")
	baseq_char_thres = baseq_thres_to_char(baseq_thres,base_qual_encoding)
	out_file_fh.write(var_file.readline().rstrip("\n")+"\tnread_alt\tnread_mismatch_alt\tmm_percent\tavg_mm_perread\n")
	var_file.close()
	out_file_fh.close()
	
	processes = [mp.Process(target=worker, args=(var_file_nm, out_file_nm, x, nproc, bam_fname, baseq_char_thres, base_qual_encoding, window_size, max_depth)) for x in range(nproc)]
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


#!/usr/bin/env python

import argparse
import os
import pysam
import math
import multiprocessing as mp
import re

def GC(seq):
	"This function computes the GC percentage of a DNA sequence."
	nbases = seq.count('n') + seq.count('N')
	gc_content = float(seq.count('c')+seq.count('C')+seq.count('g')+seq.count('G'))*100.0/(len(seq)-nbases)
	return gc_content
	
def reverse_complement(seq): 
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}    
	bases = list(seq) 
	bases = reversed([complement.get(base,base) for base in bases])
	bases = ''.join(bases)
	return bases
	
def hrun_info(ref_fname, chrom, pos):
	length = 15
	min_hlen = 3
	GC_length = 50
	pos = int(pos)
	ref = pysam.Fastafile(ref_fname)
	seq_ref_srt = GC_length
	start = (pos - 1) - GC_length 
	if start < 0:
		start = 0
		seq_ref_srt = pos - 1 
	end = pos + GC_length
	GC_seq = ref.fetch(reference=chrom, start=start, end=end).upper()
	gc_content = GC(GC_seq)
	seq_ref_srt = length
	start = (pos - 1) - length 
	if start < 0:
		start = 0
		seq_ref_srt = pos - 1 
	end = pos + length
	seq = ref.fetch(reference=chrom, start=start, end=end).upper()
	tripat = seq[length-1:(length+2)]
	up_base = tripat[0]
	down_base = tripat[2]
	pre2bases = seq[(length-2) : length]
	rev_pre2bases = reverse_complement(seq[(length+1) : (length+3)])
	ref_allele = seq[seq_ref_srt:seq_ref_srt + 1]
	# info: list of list, or R data.frame structure, with 3 columnsreverse_complement(tripat), reverse_complement(up_base), 
	reverse_complement(down_base), rev_pre2bases,
	# length of homopolymer, type(A,C,G,T), distance to variant(0 ~ len)
	pattern='A{'+ re.escape(str(min_hlen)) +',}|C{'+ re.escape(str(min_hlen)) +',}|G{'+ re.escape(str(min_hlen)) +',}|T{'+ re.escape(str(min_hlen)) +',}'
	info = [(m.group(),m.start(0), m.end(0)) for m in re.finditer(pattern, seq,flags=re.IGNORECASE)]
	nhlen = 0
	hrun_op = 0
	lhlen = 0
	lhtype = 'N'
	up_chdist = length
	up_chtype = 'N'
	up_chlen = 0
	down_chdist = length
	down_chtype = 'N'
	down_chlen = 0
	dist = 0
	hrun_op_dist = length
	hrun_op_len = 0
	for i in range(len(info)):
		hlen = info[i][2] - info[i][1]
		nhlen = nhlen + hlen
		if info[i][2] <= length:
			dist = length - info[i][2] + 1 
			if dist < up_chdist:
				up_chdist = dist
				up_chtype = info[i][0][0].upper()
				up_chlen = hlen
		elif info[i][1] > length:
			dist = info[i][1] - length
			if dist < down_chdist:
				down_chdist = dist
				down_chtype = info[i][0][0].upper()
				down_chlen = hlen
		else:
			hrun_op = 1
			hrun_op_dist = min(info[i][2]-length-1, length - info[i][1])
			hrun_op_len = info[i][2] - info[i][1]			
		if i == 0:
			lhlen = hlen
			lhtype = info[i][0][0].upper()
		elif hlen > lhlen:
			lhlen = hlen
			lhtype = info[i][0][0]
	hden = float(nhlen) / (2*length +1)
	dist_list = sorted([(up_chdist,0),(down_chdist,1),(hrun_op_dist,2)])
	len_list = [up_chlen,down_chlen,hrun_op_len]
	hmer_dist = dist_list[0][0]
	hmer_len = len_list[dist_list[0][1]]
	if dist_list[0][0] == dist_list[1][0]:
		if len_list[dist_list[1][1]] > len_list[dist_list[0][1]]:
			hmer_len = len_list[dist_list[1][1]]
	anno = ['\t'.join([chrom, str(pos-1), str(pos), ref_allele, seq, tripat, up_base, down_base, pre2bases, str(lhlen), lhtype, 
							   str(up_chdist), up_chtype,str(up_chlen),str(down_chdist), down_chtype,str(down_chlen), "{:.4}".format(gc_content),
							   "{:.4f}".format(hden), str(hrun_op), str(hrun_op_dist), str(hrun_op_len),str(hmer_dist), str(hmer_len)]), 
			'\t'.join([chrom, str(pos-1), str(pos), reverse_complement(ref_allele), reverse_complement(seq), reverse_complement(tripat), reverse_complement(up_base), 
								reverse_complement(down_base), rev_pre2bases,str(lhlen), reverse_complement(lhtype),  str(down_chdist), reverse_complement(down_chtype),str(down_chlen),str(up_chdist), reverse_complement(up_chtype),
								str(up_chlen), "{:.4f}".format(gc_content), "{:.4f}".format(hden), str(hrun_op), str(hrun_op_dist), str(hrun_op_len),str(hmer_dist), str(hmer_len)])]
	ref.close()
	return anno

def hrun_info_worker(ref, bed, output_suffix):
	bed_file = open(bed)
	fwd_out_file = open(bed + output_suffix + 'fwd', 'w')
	rev_out_file = open(bed + output_suffix + 'rev', 'w')
	for line in bed_file:
		line = line.rstrip("\n")
		cols = []
		cols = line.split("\t") # Input file should contain "chrom pos_srt pos_end".
		chrom = cols[0]
		srt = int(cols[1])
		end = int(cols[2])
		for p in range(srt, end+1):
			anno = hrun_info(ref, chrom, p) ## 0-based or 1-based
			fwd_out_file.write(anno[0] +'\n')
			rev_out_file.write(anno[1] +'\n')
	bed_file.close()
	fwd_out_file.close()
	rev_out_file.close()	

def main():
	## Parse arguments
	parser = argparse.ArgumentParser(description='Add genome contexts on given bed formated regions.', epilog= """
	The output will attach 9 columns to the original region file: \n
	longest homopolymer length: lhlen,
	longest homopolymer type: lhtype,
	the 3 mer pattern: tripat
	up_base: immediate upper base
	down_base: immediate down_base
	closest upstream homopolymer distance to variant: up_chdist,
	closest upstream homopolymer type: up_chtype,
	closest upstream homopolymer length: up_chlen,
	closest downstream homopolymer distance to variant: down_chdist,
	closest downstream homopolymer type: down_chtype,
	closest downstream homopolymer length: down_chlen,
	homopolymer density: hden, 	
	a tag, 1 for snp in a homopolymer, 0 not: hrun_op,
	if in hmer, distance to nearest boundary: hrun_op_dist,
	if in hmer, length of the homopolymer: hrun_op_len,
	closest homopolymer distance: hmer_dist,
	closest homopolymer length: hmer_len
	""")
	parser.add_argument("--ref", type=str, help="the reference genome fasta file name", required=True)
	parser.add_argument("--bed", type=str, help="the bed file name for targeted regions", required=True)
	parser.add_argument("--np", type=int, help="number of processes", default=1)
	args = parser.parse_args()
	print args
	## check parameters
	if os.path.isfile(args.ref):
		if os.path.isfile(args.bed):
			if args.np > 1:
				num_lines = sum(1 for line in open(args.bed))
				nline_pernt = int(math.ceil(num_lines / float(args.np)))
				cmd = "split -d -l %d %s %s " % (nline_pernt, args.bed, args.bed) 
				os.system(cmd)
				print "Parameters ready. Generating genome context annotation file ..."
		else: 
			print "Bed file %s is missing." % args.bed
	else:
		print "Ref file %s is missing." % args.ref
							
	## add annotation
	out_suffix = '.hrun.anno'
	if args.np == 1:
		hrun_info_worker(args.ref, args.bed, out_suffix)
	else:
		bed_fnm = []
		out_suffix_vec = []
		for i in range(args.np):
			if i < 10:
				bed_fnm.append(args.bed + '0' + str(i))
				out_suffix_vec.append(out_suffix + '_0' + str(i))
			else:
				bed_fnm.append(args.bed + str(i))
				out_suffix_vec.append(out_suffix + "_" + str(i))
		#Setup a list of processes 
		processes = [mp.Process(target=hrun_info_worker, args=(args.ref, bed_fnm[f], out_suffix_vec[f])) for f in range(args.np)]
		# Run processes
		for p in processes:
			p.start()
			# Exit the completed processes
		for p in processes:
			p.join()
		bed_path = os.path.dirname(args.bed)
		cmd = "cat " + bed_path + "/*hrun.anno_*fwd > " + args.bed + ".hrun.anno.fwd"
		os.system(cmd)
		cmd = "rm " + bed_path + "/*bed*.hrun.anno_*fwd"
		#os.system(cmd)
		cmd = "cat " + bed_path + "/*hrun.anno_*rev > " + args.bed + ".hrun.anno.rev"
		os.system(cmd)
		cmd = "rm " + bed_path + "/*bed*.hrun.anno_*rev"
		#os.system(cmd)		

			
if __name__ == "__main__":
    main()

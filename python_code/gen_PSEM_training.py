#!/usr/bin/env python

import argparse
import os
import pysam
import math
import multiprocessing as mp
import re

def strand_bcf_to_strand_PSEM(bcf, anno):
	# first expand strand specific bcf to all 4 letters, fwd and rev in seperate files
	# strand specific bcf columns: chrom srt end depth fwd_depth fwd ACGT unknown rev_depth rev ACGT
	bcf_file = open(bcf)
	fwd_fnm = bcf + '.fwd'
	fwd_out_file = open(fwd_fnm, 'w')
	rev_fnm = bcf + '.rev'
	rev_out_file = open(rev_fnm, 'w')
	bases = ['A', 'C', 'G', 'T']
	rev_bases = ['T', 'G', 'C', 'A']
	for line in bcf_file:
		line = line.rstrip("\n")
		cols = []
		cols = line.split("\t") # Input file should contain "chrom pos_srt pos_end".
		chrom = cols[0]
		srt = cols[1]
		end = cols[2]
		fwd_depth = cols[4]
		for b in range(5,9):
			fwd_out_file.write("\t".join([chrom, srt, end, cols[b], fwd_depth, bases[b-5]]) + "\n")
		rev_depth = cols[10]
		for b in range(11,15):
			rev_out_file.write("\t".join([chrom, srt, end, cols[b], rev_depth, rev_bases[b-11]]) + "\n")		
	fwd_out_file.close()
	rev_out_file.close()		
	bcf_file.close()
	## overlap with genome context annotations
	# process forward strand file
	cmd = "bedtools intersect -wo -a " + fwd_fnm + " -b " + anno + ".fwd |cut -f 1,3-6,10,12-30 | awk 'BEGIN{OFS=\"\t\"}{if($5!=$6) print $0;}' > tmp " 
	os.system(cmd)
	cmd = "awk 'BEGIN{OFS=\"\t\"}{$6=$6\"_\"$5; if($5==$8 || $5==$9){$25=$25\"\t\"1} else $25=$25\"\t\"0} {print $0}' tmp | cut -f 1-4,6-26  > " + 	bcf + ".strand_PSEM.fwd"
	os.system(cmd)
	cmd = "rm " + tmp
	os.system(cmd)
	# process reverse strand file
	cmd = "bedtools intersect -wo -a " + rev_fnm + " -b " + anno + ".rev |cut -f 1,3-6,10,12-30 | awk 'BEGIN{OFS=\"\t\"}{if($5!=$6) print $0;}' > tmp " 
	os.system(cmd)
	cmd = "awk 'BEGIN{OFS=\"\t\"}{$6=$6\"_\"$5; if($5==$8 || $5==$9){$25=$25\"\t\"1} else $25=$25\"\t\"0} {print $0}' tmp | cut -f 1-4,6-26  > " + 	bcf + ".strand_PSEM.rev"
	os.system(cmd)
	cmd = "rm " + tmp
	os.system(cmd)

def main():
	## Parse arguments
	parser = argparse.ArgumentParser(description='Run strand specific pileup of the bam file on given bed formated regions.')
	parser.add_argument("--bcf", type=str, help="the strand specific basecall file name", required=True)
	parser.add_argument("--anno", type=str, help="the genome context annotation file name base, without fwd or rev suffix", required=True)
	args = parser.parse_args()
	## check parameters
	if os.path.isfile(args.bcf):
		if os.path.isfile(args.anno):
			print "Parameters ready. Generating pileup file ..."
		else: 
			print "anno file %s is missing." % args.bed
	else:
		print "bcf file %s is missing." % args.ref
	strand_bcf_to_strand_PSEM(bcf, anno)


if __name__ == "__main__":
	main()

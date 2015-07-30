#!/usr/bin/env python

import argparse
import os

def PSEM_filter(fnm, min_dep, max_dep, min_af, max_af, min_nalt):
	"Filter based on depth and allele frequency"
	cmd = "awk 'BEGIN{OFS=\"\t\"}{if($4>="+min_dep+" && $4<="+max_dep+" && $3>="+min_nalt+" && $3/($3+$4)>="+min_af+" && $3/($3+$4)<="+max_af+")}' " + fnm + " > " + fnm + ".filtered"
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
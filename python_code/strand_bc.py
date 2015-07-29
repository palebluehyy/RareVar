#!/usr/bin/env python

import argparse
import os
import math
import multiprocessing as mp
import re


def strand_bc(ref, bam, bed, out_suffix):
	bam_prefix = re.sub('.bam$', '', bam)
	cmd = "samtools mpileup -d 1000000 -s -l " + bed + " -f " + ref + " " + bam + " > " + bam_prefix + ".pileup" + out_suffix
	print "Samtools mpileup on bed region: %s ..." % bed
	os.system(cmd)
	print "Strand Specific Call on bed region: %s ..." % bed
	cmd = "java -cp /home/bin/htsvipr pileup2vipr " + bam_prefix + ".pileup" + out_suffix + " " + bam_prefix + ".pileup_parsed" + out_suffix
	os.system(cmd)
	cmd = "rm " + bam_prefix + ".pileup" + out_suffix

def main():
	## Parse arguments
	parser = argparse.ArgumentParser(description='Run strand specific pileup of the bam file on given bed formated regions.')
	parser.add_argument("--ref", type=str, help="the reference genome fasta file name", required=True)
	parser.add_argument("--bam", type=str, help="the bam file name", required=True)
	parser.add_argument("--bed", type=str, help="the bed file name for targeted regions", required=True)
	parser.add_argument("--np", type=int, help="number of processes (or chuncks to do pileup seperately)", default=1)
	args = parser.parse_args()
	## check parameters
	if os.path.isfile(args.ref):
		if os.path.isfile(args.bam):
			if os.path.isfile(args.bed):
				if args.np > 1:
					num_lines = sum(1 for line in open(args.bed))
					nline_pernt = int(math.ceil(num_lines / float(args.np)))
					cmd = "split -d -l %d %s %s " % (nline_pernt, args.bed, args.bed) 
					os.system(cmd)
				print "Parameters ready. Generating pileup file ..."
			else: 
				print "Bed file %s is missing." % args.bed
		else:
			print "Bam file %s is missing." % args.bam
	else:
		print "Ref file %s is missing." % args.ref
	## pileup
	if args.np == 1:
		strand_bc(args.ref, args.bam, args.bed)
	else:
		bed_fnm = []
		out_suffix_vec = []
		for i in range(args.np):
			if i < 10:
				bed_fnm.append(args.bed + '0' + str(i))
				out_suffix_vec.append('0' + str(i))
			else:
				bed_fnm.append(args.bed + str(i))
				out_suffix_vec.append(str(i))
	#Setup a list of processes 
	processes = [mp.Process(target=strand_bc, args=(args.ref, args.bam, bed_fnm[f], out_suffix_vec[f])) for f in range(args.np)]
	# Run processes
	for p in processes:
		p.start()
		# Exit the completed processes
	for p in processes:
		p.join()
			
	
	

if __name__ == "__main__":
	main()

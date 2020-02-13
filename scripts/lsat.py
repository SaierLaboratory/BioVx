#!/usr/bin/env python
#not to be confused with some expensive entrance exam for lawyers

from __future__ import print_function, division, unicode_literals

import argparse
import subprocess
import shlex
#import matplotlib
#matplotlib.use('Qt4Agg')
import io

import Bio.SeqIO

import shuffle_sequences_within_class as shuffleseq

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('--ssearch', default='ssearch36', help='Name/location of the SSEARCH binary (if not in $PATH as "ssearch36") (default: ssearch36)')
	parser.add_argument('--flags', default=[], nargs='+', help='Flags to pass on to SSEARCH. -E will be set to the size of the shuffle library for the sake of completeness')
	parser.add_argument('--prefix', default='lsatout', help='Prefix for output files (default:lsatout)')
	parser.add_argument('-c', '--count', type=int, default=1000, help='Number of shuffles to perform (default:1000)')
	parser.add_argument('--buffersize', default=1000, help='Number of sequences to place in buffer before writing to file')

	parser.add_argument('query', help='Query sequence. Use "asis:" to enter raw sequences')
	parser.add_argument('target', help='Target sequence. Use "asis:" to enter raw sequences  (will be shuffled)')
	#parser.add_argument('background', help='Library of sequences to compare the unshuffled alignments to')

	parser.add_argument('-s', action='store_true', help='Split by in/out')
	parser.add_argument('-t', action='store_true', help='Split by loop/tail') 
	args = parser.parse_args()

	queryfasta = ''
	if args.query.startswith('asis:'): queryfasta = args.query[5:]
	else: 
		with open(args.query) as f: queryfasta += f.read()
	if not queryfasta.startswith('>'): queryfasta = '>query\n{}'.format(queryfasta)
	queryseq = Bio.SeqIO.read(io.StringIO(queryfasta), 'fasta')

	targetfasta = ''
	if args.target.startswith('asis:'): targetfasta = args.target[5:]
	else:
		with open(args.target) as f: targetfasta += f.read()
	if not targetfasta.startswith('>'): targetfasta = '>target\n{}'.format(targetfasta)
	targetseq = Bio.SeqIO.read(io.StringIO(targetfasta), 'fasta')



	with open('{}_query.faa'.format(args.prefix), 'w') as f:
		f.write('>{}\n{}\n'.format(queryseq.name, queryseq.seq))

	f = open('{}_targetshuffled.faa'.format(args.prefix), 'w')

	f.write('>{}\n{}\n'.format(targetseq.name, targetseq.seq))

	mergedict = shuffleseq.get_mergedict(args.s, args.t)

	remaining = args.count
	for i in range(0, args.count, args.buffersize):
		if remaining < args.buffersize: count = remaining % args.buffersize
		else: count = args.buffersize

		sequences = shuffleseq.shuffle_seq(targetseq, count=count, prefix=args.prefix, mergedict=mergedict)
		f.write('\n'.join(sequences))
		remaining -= count

	f.close()
	f = open('{}_shuffled.ssearch'.format(args.prefix), 'w')
	cmd = [args.ssearch, '-E', str(args.count)]
	cmd.extend(args.flags)
	cmd.append('{}_query.faa'.format(args.prefix))
	cmd.append('{}_targetshuffled.faa'.format(args.prefix))
	out = subprocess.check_output(cmd)
	f.write(out)

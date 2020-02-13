#!/usr/bin/env python
from __future__ import print_function

import argparse
import random
import subprocess
import re
import Bio.SeqIO

def hmmtop(seq):
	
	topo = None
	indices = []

	p = subprocess.Popen(['hmmtop', '-if=--', '-is=pseudo', '-pi=spred', '-sf=FAS'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	out, err = p.communicate(input='>seq\n{}'.format(seq))

	relevant = re.findall('(?:IN|OUT)\s+(?:[0-9]+\s*)+$', out)[0]
	if relevant.startswith('IN'): topo = 1
	elif relevant.startswith('OUT'): topo = 0
	else: topo = None
	x = relevant.split()
	indices = [[int(x[i])-1, int(x[i+1])-1] for i in range(2, len(x), 2)]


	return topo, indices

def gen_states(topo, indices, seqlen):

	if topo == 0:
		start = 'O'
		inverse = 'I'
	elif topo == 1:
		start = 'I'
		inverse = 'O'
	else: raise ValueError

	states = [start for i in range(seqlen)]

	for i, tms in enumerate(indices):
		for j in range(tms[0], tms[1]+1): 
			states[j] = 'H'

	for i, tms in enumerate(indices):
		for j in range(tms[0]-1, 0, -1): 
			if states[j] == 'H': break
			elif (i % 2): #before an even (note 1 vs. 0-indexing) TMS: inverse
				states[j] = inverse
			elif not (i % 2): #before an odd (note 1 vs. 0-indexing) TMS:normal 
				states[j] = start
		for j in range(tms[1]+1, seqlen, 1):
			if states[j] == 'H': break
			elif (i % 2): #before an even (note 1 vs. 0-indexing) TMS: inverse
				states[j] = start 
			elif not (i % 2): #before an odd (note 1 vs. 0-indexing) TMS:normal 
				states[j] = inverse

	for i, tms in enumerate(indices):
		for j in range(tms[1]+1, min(seqlen, tms[1]+1+15), 1):
			if states[j] == 'H': break
			else: states[j] = states[j].lower()

		for j in range(tms[0]-1, max(0, tms[0]-1-15), -1): 
			if states[j] == 'H': break
			else: states[j] = states[j].lower()

	return states


def get_shuffle(template, mergedict=None):
	shuffleme = {}
	for k in template:
		if not mergedict:
			shuffleme[k] = template[k][:]
		else:
			if k in mergedict:
				shuffleme[k] = shuffleme[mergedict[k]]
			else:
				shuffleme[k] = template[k][:]
	return shuffleme


def get_pieces(seq, mergedict=None):
	topo, indices = hmmtop(seq.seq)
	states = gen_states(topo, indices, len(seq))
	pieces = {}
	for k in 'HIiOo': pieces[k] = ''
	for resn, state in zip(seq, states):
		pieces[state] += resn

	template = {}
	for k in 'HIiOo':
		if mergedict: 
			if k in mergedict:
				try: template[mergedict[k]].extend(list(pieces[k]))
				except KeyError: template[mergedict[k]] = list(pieces[k])
				template[k] = template[mergedict[k]]
			else: template[k] = list(pieces[k])
		else: template[k] = list(pieces[k])
		#template[k] = list(pieces[k])
	return template, states

def shuffle_seq(seq, count=1, prefix='shuffled', mergedict=None):
	sequences = []
	template, states = get_pieces(seq, mergedict=mergedict)
	for i in range(count):
		shuffleme = get_shuffle(template, mergedict=mergedict)
		for k in shuffleme:
			if k in mergedict: continue
			random.shuffle(shuffleme[k])
		newseq = ''

		#for k in shuffleme: print(k, len(template[k]), len(shuffleme[k]), states.count(k))

		for state in states:
			newseq += shuffleme[state].pop(0)
		sequences.append('>{}_{:016X}\n{}'.format(prefix, abs(hash(newseq)), newseq))
		del shuffleme
	return sequences

def get_mergedict(inout=False, looptail=False):
	#TODO: optimize away all these annoyingly hardcoded lines
	#split both by tail/loop and in/out
	if inout and looptail: mergedict = {}
	#split by in/out but not tail/loop
	elif inout: mergedict = {'i':'I', 'o':'O'}
	#split by tail/loop but not in/out
	elif looptail: mergedict = {'O':'I', 'o':'i'}
	#everything outside a TMS is a loop
	else: mergedict = {'i':'I', 'O':'I', 'o':'I'}
	return mergedict


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-e', action='store_true', help='Echo seed sequence in output as first sequence')
	parser.add_argument('-n', default=1, type=int, help='Number of shuffles (default:1)')
	parser.add_argument('-o', default='/dev/stdout', help='Where to write the shuffled sequences (default:stdout)')
	parser.add_argument('-p', default='shuffled', help='Sequence label prefix (default:shuffled)')
	parser.add_argument('-s', action='store_true', help='Split I/i from O/o when shuffling to keep loop topology')
	parser.add_argument('-t', action='store_true', help='Split I from i and O from o when shuffling to retain some information on distance from the membrane')
	parser.add_argument('infile', nargs='?', default='/dev/stdin', help='File to read in (default:stdin)')
	args = parser.parse_args()


	mergedict = get_mergedict(args.s, args.t)

	f = open(args.o, 'w')

	fasta = Bio.SeqIO.parse(args.infile, 'fasta')

	bufsize = 100

	for seq in fasta:
		remaining = args.n
		if args.e: f.write('>{}\n{}\n'.format(seq.name, seq.seq))

		for i in range(0, args.n, bufsize):
			if remaining < bufsize: count = remaining % bufsize
			else: count = bufsize

			sequences = shuffle_seq(seq, count=count, prefix=args.p, mergedict=mergedict)
			f.write('\n'.join(sequences))
			remaining -= count


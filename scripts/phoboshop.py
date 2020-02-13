#!/usr/bin/env python

from __future__ import print_function, division
import argparse

import os
if 'MPLBACKEND' not in os.environ: os.environ['MPLBACKEND'] = 'Qt4Agg'
import quod
import numpy as np

import sys
import subprocess
import re
import shlex


class Protein(object):
	fasta = ''
	entities = {}
	def __init__(self, fasta):
		if not fasta.startswith('>'): raise ValueError('fasta is not in FASTA format')
		elif '\n' not in fasta: raise ValueError('fasta is not in FASTA format')
		self.fasta = fasta
		self.header = self.get_header()
		self.sequence = self.get_sequence()

	def get_entropy(self, a=None, b=None):
		seq = self.get_sequence(a, b)

		#seq = seq[int(np.floor(a)):int(np.ceil(b))]
		total = 0
		counts = {}
		for resn in seq:
			try: counts[resn] += 1
			except KeyError: counts[resn] = 1
			total += 1
		entropy = 0
		for k in counts: entropy -= counts[k]/total*np.log2(counts[k]/total)
		return entropy
			

	def get_sequence(self, a=None, b=None):
		seq = re.sub('[^A-Za-z\-]', '', self.fasta[self.fasta.find('\n'):])

		a = 0 if a is None else a
		a = max(0, a-1)
		b = len(seq) if b is None else b
		seq = seq[int(np.floor(a)):int(np.floor(b))]
		return seq


	def get_header(self):
		return self.fasta[:self.fasta.find('\n')]

	def render(self, style=0):
		what = quod.What(self.get_sequence(), style=style)
		self.entities['what'] = what

		entlist = []
		for k in self.entities: entlist.append(self.entities[k])
		return entlist

class Phoboshop(object):

	proteins = []
	mode = 'normal'
	fig = None
	plot = None

	outfile = '/dev/stdout'

	pipeto = None

	testline = None
	selections = {}
	queue = {}

	def __init__(self):
		pass

	def add_protein(self, protein):
		self.proteins.append(protein)

	def run(self):
		print('Use [?] to get help on phoboshop shortcuts', file=sys.stderr)
		self.fig = quod.plt.figure()
		self.plot = quod.Plot(fig=self.fig)
		for i, p in enumerate(self.proteins): 
			for e in p.render(style=i): 
				self.plot.add(e)

		self.update_title()
		self.plot.ax.figure.canvas.mpl_connect('button_press_event', self.onmousedown)
		self.plot.ax.figure.canvas.mpl_connect('button_release_event', self.onmouseup)
		self.plot.ax.figure.canvas.mpl_connect('motion_notify_event', self.onmousemove)
		self.plot.ax.figure.canvas.mpl_connect('key_press_event', self.onkeydown)
		#self.testline = self.plot.ax.plot([], [])
		self.testline, = self.plot.ax.plot([], [])

		self.plot.render()
		quod.plt.show()

	def onmousedown(self, event):
		if not event.inaxes: return

		if self.mode.startswith('cut'):
			self.mode = 'cutdrag'
			if 'cut' in self.selections and self.selections['cut'] is not None:
				vspan = self.selections['cut']
			else: 
				vspan = self.plot.ax.axvspan(event.xdata, event.xdata, fc='red', alpha=0.4, hatch='/')
				self.selections['cut'] = vspan
			xy = vspan.get_xy()
			xy[:,0] = event.xdata
			xy[2:4,0] = event.xdata + 1
			vspan.set_xy(xy)
			vspan.figure.canvas.draw()

		elif self.mode == 'delete':
			if 'cut' in self.queue:
				popme = []
				for i, vspan in enumerate(self.queue['cut']):
					xy = vspan.get_xy()
					xmin = min(xy[0,0], xy[2,0])
					xmax = max(xy[0,0], xy[2,0])
					if xmin <= event.xdata <= xmax: 
						vspan.remove()
						popme.append(i)
						self.plot.ax.figure.canvas.draw()
				for i in popme[::-1]: self.queue['cut'].pop(i)
			if 'cut' in self.selections and self.selections['cut'] is not None:
				vspan = self.selections['cut']
				xy = vspan.get_xy()
				xmin = xy[0,0]
				xmax = xy[2,0]
				if xmin <= event.xdata <= xmax: 
					vspan.remove()
					del vspan
					self.plot.ax.figure.canvas.draw()

		elif self.mode == 'entropy':
			vspan = self.plot.ax.axvspan(event.xdata, event.xdata+1, fc='g', alpha=0.2, hatch='\\')
			self.mode = 'entropydrag'
			self.selections['entropy'] = vspan
			vspan.figure.canvas.draw()

		elif self.mode == 'move':
			if 'cut' in self.queue:
				for vspan in self.queue['cut']:
					xy = vspan.get_xy()
					
					xmin = min(xy[0,0], xy[2,0])
					xmax = max(xy[0,0], xy[2,0])
					#edges
					if xmax-2 <= event.xdata <= xmax+2:
						self.mode = 'moveresizedrag'
						self.selections['resize'] = (vspan, 'r')
						return
					elif xmin-2 <= event.xdata <= xmin+2:
						self.mode = 'moveresizedrag'
						self.selections['resize'] = (vspan, 'l')
						return
					elif xmin <= event.xdata <= xmax:
						self.mode = 'movedrag'
						self.selections['move'] = (vspan, event.xdata)
						return
			
	def onmouseup(self, event):
		#if self.mode.startswith('cut'): self.mode = 'cut'
		if self.mode == 'cutdrag':
			try: self.queue['cut']
			except KeyError: self.queue['cut'] = []

			vspan = self.selections['cut']
			xy = vspan.get_xy()
			xmin = min(xy[:,0])
			xmax = max(xy[:,0])
			xy[:,0] = xmin
			xy[2:4,0] = xmax

			self.queue['cut'].append(self.selections['cut'])
			self.selections['cut'] = None
			self.mode = 'cut'
		elif self.mode == 'entropydrag':
			vspan = self.selections['entropy']

			xy = vspan.get_xy()
			xmin = int(np.floor(min(xy[:,0])))
			xmax = int(np.ceil(max(xy[:,0])))
			for p in self.proteins:
				e = p.get_entropy(xmin, xmax)
				length = len(p.get_sequence(xmin, xmax))
				if length == 0:
					print('{}: {:0.3f} bits/aa({}-{}: {} aa, max entropy: {:0.3f} bits/aa)'.format(
						p.header,
						e, 
						xmin, xmax, 
						length,
						min(-np.log2(1/20), 0)
					))
				else:
					print('{}: {:0.3f} bits/aa({}-{}: {} aa, max entropy: {:0.3f} bits/aa)'.format(
						p.header,
						e, 
						xmin, xmax, 
						length,
						min(-np.log2(1/20), -np.log2(1/(length)))
					))
				#print(p.get_sequence(xmin, xmax))
				

			vspan.remove()
			del vspan
			self.selections['entropy'] = None
			self.mode = 'entropy'
			self.plot.ax.figure.canvas.draw()

		elif self.mode == 'movedrag':
			self.selections['move'] = None
			self.mode = 'move'
		elif self.mode == 'moveresizedrag':

			vspan, side = self.selections['resize']
			xy = vspan.get_xy()
			xmin = min(xy[:,0])
			xmax = max(xy[:,0])
			xy[:,0] = xmin
			xy[2:4,0] = xmax

			self.selections['resize'] = None
			self.mode = 'move'
		#if self.mode.startswith('cut'): self.mode = 'cut'


	def onmousemove(self, event):
		if self.mode == 'cutdrag':
			if event.xdata is None: return
			if 'cut' in self.selections:
				vspan = self.selections['cut']
				xy = vspan.get_xy()
				xy[2:4,0] = event.xdata
				vspan.set_xy(xy)
				vspan.figure.canvas.draw()

		elif self.mode == 'entropydrag':
			if event.xdata is None: return
			vspan = self.selections['entropy']
			xy = vspan.get_xy()
			xy[2:4,0] = event.xdata
			vspan.set_xy(xy)
			vspan.figure.canvas.draw()

		elif self.mode == 'movedrag':
			if event.xdata is None: return
			vspan, orig = self.selections['move']
			xy = vspan.get_xy()
			xy[:,0] += (event.xdata - orig)
			vspan.set_xy(xy)
			vspan.figure.canvas.draw()
			self.selections['move'] = (vspan, event.xdata)

		elif self.mode == 'moveresizedrag':
			if event.xdata is None: return
			vspan, side = self.selections['resize']
			xy = vspan.get_xy()
			if side == 'l':
				xy[:2,0] = event.xdata
				xy[4:,0] = event.xdata
			elif side == 'r':
				xy[2:4,0] = event.xdata
			vspan.set_xy(xy)
			vspan.figure.canvas.draw()
		else: return


	def _merge(self, spanlist):
		def _overlap(span1, span2):
			# |---|
			#  |-|
			if (span1[0] <= span2[0]) and (span2[1] <= span1[1]): return True
			#  |-|
			# |---|
			elif (span2[0] <= span1[0]) and (span1[1] <= span2[1]): return True
			#  |---|
			# |---|
			elif (span2[0] <= span1[0] <= span2[1]): return True
			elif (span1[0] <= span2[1] <= span1[1]): return True
			# |---|
			#  |---|
			elif (span1[0] <= span2[0] <= span1[1]): return True
			elif (span2[0] <= span1[1] <= span2[1]): return True
			else: return False
		def _check_overlaps(spanlist):
			for i, span1 in enumerate(spanlist):
				for span2 in spanlist[i+1:]:
					if _overlap(span1, span2): return True
			return False

		spanlist.sort()
		while _check_overlaps(spanlist):
			popme = []
			for i, span1 in enumerate(spanlist):
				for j, span2 in enumerate(spanlist):
					if j <= i: continue
					if _overlap(span1, span2): 
						span1[1] = max(span1[1], span2[1])
						popme.append(j)
			#for j in popme[::-1]: spanlist.pop(j)
			spanlist.pop(popme[-1])
		return spanlist
					
		

	def write_cuts(self):
		spans = []
		if 'cut' in self.queue:
			for vspan in self.queue['cut']:
				xy = vspan.get_xy()
				xmin = int(round(min(xy[:,0])))
				xmax = int(round(max(xy[:,0])))

				span = [xmin, xmax]

				spans.append(span)
		else: spans = [[-1,-1]]
		if not spans: spans = [[-1,-1]]
		
		spans = self._merge(spans)
		
		s = ''
		for p in self.proteins:
			s += p.get_header() + '\n'
			laststart = 0
			seq = p.get_sequence()
			for span in spans:
				s += seq[laststart:span[0]]
				laststart = span[1]
			if laststart != 0: s += seq[laststart:]
			s += '\n'
		s += '\n'
		if self.pipeto:
			p = subprocess.Popen(shlex.split(self.pipeto), stdin=subprocess.PIPE)
			out, err = p.communicate(input=s)
		with open(self.outfile, 'w') as f:
			f.write(s)
				

	def print_help(self):
		#p, = self.plot.ax.plot([0,10,20,30,40], [1,2,1,0,-1])
		#p.figure.canvas.draw()
		#p.remove()
		print('''
CONTROLS
========

c - Enter cutting mode. Click and drag to mark unwanted segments
d - Enter selection-deleting mode. As a double negative is a positive, clicking on segments marked for deletion marks them as wanted once again
e - Enter entropy mode. Click and drag to probe for local sequence complexity
m - Enter move/resize mode. Click and drag cores of selections to move them or edges to resize them
w - Write cut sequence(s) to disk/stdout and run the selected pipeto program if defined
? - Pull up this help page
ESC - Enter normal mode, which does absolutely nothing
	''', file=sys.stderr)

	def onkeydown(self, event):
		#prevent mode-switching while dragging stuff around
		if 'drag' in self.mode: return

		if event.key == 'c':
			self.mode = 'cut'
		elif event.key == 'd':
			self.mode = 'delete'
		elif event.key == 'e':
			self.mode = 'entropy'
		elif event.key == 'escape':
			self.mode = 'normal'

		elif event.key == 'm':
			self.mode = 'move'

		elif event.key == 'w':
			self.write_cuts()

		elif event.key == '?':
			self.print_help()

		#elif event.key == 'm':
		#	if event.inaxes != self.plot.ax: return
		#	print(event.xdata, event.ydata)
		#	xs = list(self.testline.get_xdata())
		#	xs.append(event.xdata)
		#	self.testline.set_xdata(xs)
		#	ys = list(self.testline.get_ydata())
		#	ys.append(event.ydata)
		#	self.testline.set_ydata(ys)
		#	self.testline.figure.canvas.draw()
		#	#self.plot.ax.

		self.update_title()

	def update_title(self):
		if self.mode == 'normal': self.fig.canvas.set_window_title('Normal mode')
		elif self.mode.startswith('cut'): self.fig.canvas.set_window_title('SELECT-FOR-DELETION mode')
		elif self.mode.startswith('move'): self.fig.canvas.set_window_title('MOVE-SELECTION mode')
		elif self.mode.startswith('entropy'): self.fig.canvas.set_window_title('ENTROPY mode')
		elif self.mode == 'delete': self.fig.canvas.set_window_title('DESELECT mode')

def split_fasta(f):
	firstline = True
	sequences = []
	for l in f:
		if firstline and not l.startswith('>'): sequences.append('>sequence\n')
		firstline = False
		if l.startswith('>'):
			sequences.append(l)
		else:
			sequences[-1] += l
	return sequences

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('infile', nargs='+', help='Sequence files to read')
	parser.add_argument('-o', metavar='OUTFILE', default='/dev/stdout', help="Where to save the cut sequence(s). (default: print to stdout)")
	parser.add_argument('-p', metavar='PROGNAME | "PROGNAME ARG1 ARG2..."', help="Where to pipe the output (in addition to saving the cut sequence(s). Specify `xclip' to send stuff to the clipboard on Linux and `pbcopy' to send stuff to the clipboard on Mac. If the target program requires multiple arguments, enclose PROGNAME in quotes")

	args = parser.parse_args()


	phoboshop = Phoboshop()
	phoboshop.outfile = args.o

	if args.p: phoboshop.pipeto = args.p

	for infn in args.infile:
		with open(infn) as f:
			for fasta in split_fasta(f):
				p = Protein(fasta=fasta)
				phoboshop.add_protein(p)
	phoboshop.run()

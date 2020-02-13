#!/usr/bin/env python
from __future__ import print_function, division, unicode_literals

import argparse
import subprocess
import tempfile
import re

import os

import numpy as np

from Bio import AlignIO, SeqIO

import matplotlib.transforms
import quod
if __name__ == '__main__':
	backend = os.environ.get('MPLBACKEND', 'Agg')
	quod.plt.switch_backend(backend)

def pickmax(a, b, c):
	maxfont = 0
	if a is None: pass
	else: maxfont = a
	if b is None: pass
	else: maxfont = max(maxfont, b)
	if c is None: pass
	else: maxfont = max(maxfont, c)

	if maxfont: return maxfont
	else: return None
	
def parse_hmmtop_for_centers(l, correction):
	indices = re.findall('(?: IN| OUT)((?:\s*(?:[0-9]+))+)', l.strip())[0].strip().split()
	indices = [int(i) for i in indices[1:]]
	centerlist = []
	for i in range(0, len(indices), 2):
		center = ((indices[i] - correction[i]) + (indices[i+1] - correction[i+1])) / 2
		centerlist.append(center)
	return centerlist

def get_tmcenters(msa):
	tf = tempfile.NamedTemporaryFile()
	SeqIO.write(msa, tf.name, 'fasta')
	tf.flush()

	corrlist = []
	for seq in msa:
		correction = np.zeros(len(seq))
		delta = 0
		for i, resn in enumerate(seq): 
			if resn not in 'ACDEFGHIKLMNPQRSTVWY': delta += 1
			correction[i] = delta
		corrlist.append(correction)

	out = subprocess.check_output(['hmmtop', '-if={}'.format(tf.name), '-pi=spred', '-is=pseudo', '-sf=FAS'])
	tmcenters = []
	for l, correction in zip(out.split('\n'), corrlist):
		if not l: continue

		tmcenters.extend(parse_hmmtop_for_centers(l, correction))
	
	return tmcenters

def get_similarities(msa, hydro_wt=1., gap_wt=1., index=quod.HYDROPATHY):
	gap_counts = np.zeros(msa.get_alignment_length())
	hydro_stdevs = np.zeros(msa.get_alignment_length())
	gap_costs = np.zeros(msa.get_alignment_length())

	for i in range(msa.get_alignment_length()):
		position = []
		gaps = 0
		for seq in msa:
			resn = seq[i]
			if resn != '-': position.append(index.get(resn, np.nan))
			else: gaps += 1
		gap_costs[i] = (len(msa) - gaps - 1) / (len(msa) - 1)

		score = -np.nanstd(position)
		hydro_stdevs[i] = score

	minstd, maxstd = min(hydro_stdevs), max(hydro_stdevs)

	norm_hydro_stdevs = (hydro_stdevs - minstd) / (maxstd - minstd)
	#print(norm_hydro_stdevs)

	scores = norm_hydro_stdevs**hydro_wt * gap_costs**gap_wt
	#print(scores)

	#exit()
	return scores

def get_average_hydropathies(msa, window=19, index=quod.HYDROPATHY, kernel='moving'):
	positions = np.zeros(msa.get_alignment_length())#-window+1)
	occupancies = np.zeros(msa.get_alignment_length())#-window+1)

	for seq in msa:
		for i, resn in enumerate(seq):
			if resn in index: 
				if index[resn]:
					positions[i] += index[resn]
					occupancies[i] += 1

	positions /= occupancies
	rawhydro, gaps = quod.collapse(positions)
	hydro = np.convolve(rawhydro, quod.get_kernel(kernel, window), 'valid')
	result = quod.uncollapse(list(hydro), gaps)
	result = np.hstack([[np.nan]*(window//2), result, [np.nan]*(window - (window//2))])
	return result

def get_average_amphipathicities(msa, window=19, index=quod.AMPHIPATHICITY, angle=100):
	positions = np.zeros(msa.get_alignment_length())
	occupancies = np.zeros(msa.get_alignment_length())

	cosines = np.zeros(msa.get_alignment_length())
	sines = np.zeros(msa.get_alignment_length())

	for seq in msa:
		for i in range(msa.get_alignment_length()):
			if seq[i] != '-':
				amphi = index.get(seq[i], 0)
				#amphipathicities.append(amphi)

				cosines[i] += amphi * np.cos(i * np.pi / 180 * angle)
				sines[i] += amphi * np.sin(i * np.pi / 180 * angle)
				occupancies[i] += 1

	amphipathicities = np.zeros(msa.get_alignment_length() - window)
	for i in range(msa.get_alignment_length() - window):
		amphipathicities[i] = (occupancies[i]/len(msa) * (np.mean(cosines[i:i+window])**2 + np.mean(sines[i:i+window])**2)**.5) / window
	
	amphipathicities *= 5/2 / len(msa) * window #somehow canceling out window produces a more faithful replication
	amphipathicities = np.hstack([[np.nan]*(window//2), amphipathicities])
	return amphipathicities

class TMcenter(quod.Vspans):
	def __init__(self, centerlist=None, ymin=-3, ymax=-2.5, linewidth=1.5, color='k', alpha=0.2):
		self.centerlist = [] if centerlist is None else centerlist
		self.spans = self.centerlist
		self.ymin = ymin
		self.ymax = ymax
		self.linewidth = linewidth
		self.color = color
		self.style = color
		self.alpha = alpha

	def get_bounding_box(self):
		if not self.centerlist: return None
		else: return np.array([min(self.centerlist), self.ymin, max(self.centerlist) - min(self.centerlist), self.ymax - self.ymin])

	def draw(self, plot):
		axvlines = []
		#for x in self.centerlist:
		#	axvlines.append(plot.ax.axvline(x=x, ymin=self.ymin, ymax=self.ymax-self.ymin, linewidth=self.linewidth, color=self.style, alpha=self.get_alpha()))
		#return axvlines
		return plot.ax.broken_barh([(x,0) for x in self.centerlist], (self.ymin, self.ymax-self.ymin), linewidth=self.linewidth, edgecolor=self.style, facecolor=self.style, alpha=self.get_alpha())

def plot_tmcenters(centerlist, ax, ymin, ymax, width=0.5):
	xranges = [[center-width/2, width] for center in centerlist]
	yrange = [ymin, ymax-ymin]
	#fc = (1.0, 0.5, 0.25, 0.25)
	fc = (0., 0., 0., 0.25)
	ec = (0., 0., 0., 0.)
	#return ax.broken_barh(xranges, yrange, facecolor=fc, edgecolor=ec)
	axvlines = []
	for x in centerlist:
		axvlines.append(ax.axvline(x=x, ymin=-2, ymax=0.1, linewidth=1.0, color='k', alpha=0.2))
	return axvlines

def plot_similarities(simlist, ax, ymin, ymax, window=10):
	scores = np.array([np.nanmean(simlist[i:i+window]) for i in range(0, len(simlist)-window)])
	X = np.arange(0, len(simlist)-window) + window/2
	color = 'gray'

	normscores = (scores - min(scores)) / (max(scores) - min(scores))
	scaledscores = normscores * (ymax - ymin)
	shiftedscores = scaledscores + ymin

	return ax.plot(X, shiftedscores, color, lw=1.0)

def plot_hydropathies(hydropathies, ax):
	return ax.plot(hydropathies, 'r', lw=1.0)

def plot_amphipathicities(amphipathicities, ax):
	return ax.plot(amphipathicities, 'g', lw=1.0)

class Avehas(quod.What):
	def __init__(self, f, **kwargs):

		self.window = kwargs.get('window', 19)
		self.kernel = kwargs.get('kernel', 'moving')
		self.simwindow = kwargs.get('simwindow', 10)

		self.grid = kwargs.get('grid', False)

		self.xlabel = kwargs.get('xlabel', 'Position')
		self.ylabel = kwargs.get('ylabel', 'Relative values')

		self.axisfont = kwargs.get('axisfont', None)

		self.xticks = kwargs.get('xticks', None)
		self.tickfont = kwargs.get('tickfont', None)

		self.ltitlefont = kwargs.get('ltitlefont', None)
		self.ltitle = kwargs.get('ltitle', None)
		self.ctitlefont = kwargs.get('ctitlefont', None)
		self.ctitle = kwargs.get('ctitle', None)
		self.rtitlefont = kwargs.get('rtitlefont', None)
		self.rtitle = kwargs.get('rtitle', None)

		self.amphi = kwargs.get('amphi', False)
		self.entdict = {'hydro': quod.Curve(style=0), 'tms': quod.Vspans(), 'amphi': quod.Curve(style=2), 'simil': quod.Curve(style='gray'), 'tmcenter': TMcenter([])}
		#self.entities = [quod.Curve(style=0), quod.Vspans(), quod.Curve(style=0), TMcenter([])]
		self.entities = [self.entdict['hydro'], self.entdict['tms'], self.entdict['amphi'], self.entdict['simil'], self.entdict['tmcenter']]
		self.linecolor = 'red'
		self.tmscolor = 'orange'

		#attribs subject to change!
		self.linecolorhydro = 'red'
		self.linecolorsimil = 'gray'
		self.linewidthsimil = 1.0
		self.tmscolorcenter = '#00000033'

		self.windowsimil = 10

		self.length = 0

		self.parse_alignment(f)

		self.entdict['hydro'].Y = self.hydropathies
		self.entdict['hydro'].X = np.arange(0, len(self.hydropathies))
		self.entdict['hydro'].style = self.linecolorhydro
		bounds = self.entdict['hydro'].get_bounding_box()
		bounds[3] = 6
		bounds[1] = -3
		self.entdict['hydro'].bounding_box = bounds

	#def get_bounding_box(self):
		#left = 0
		#right = len(self.hydropathies)
		#bottom, top = -3, 3
		#return np.array([left, bottom, right-left, top-bottom])

	def draw(self, plot):
		ax = plot.ax


		#plot_similarities(self.similarities, ax, ymin=-3, ymax=-1.5)
		#scores = np.array([np.nanmean(self.similarities[i:i+self.windowsimil]) for i in range(0, len(self.similarities)-self.windowsimil)])
		scores = quod.nanconvolve(self.similarities, quod.get_kernel(self.kernel, self.windowsimil), 'valid')
		X = np.arange(0, len(self.similarities)-self.windowsimil) + self.windowsimil/2
		color = 'gray'

		normscores = (scores - min(scores)) / (max(scores) - min(scores))
		scaledscores = normscores * (-1.5 - -3)
		shiftedscores = scaledscores + -3
		self.entdict['simil'].X = X
		self.entdict['simil'].Y = shiftedscores

		self.entdict['simil'].style = self.linecolorsimil
		self.entdict['simil'].linewidth = self.linewidthsimil

		#plot_tmcenters(self.tmcenters, ax, ymin=-3, ymax=-2.5)
		self.entdict['tmcenter'].centerlist = self.tmcenters
		self.entdict['tmcenter'].ymin = -3
		self.entdict['tmcenter'].ymax = -2.5
		self.entdict['tmcenter'].style = self.tmscolorcenter
		

		for e in self.entities: e.draw(plot)

		if self.grid: ax.grid(1)

		if self.xlabel: ax.set_xlabel(self.xlabel)
		if self.ylabel: ax.set_xlabel(self.ylabel)

		if self.axisfont:
			ax.set_xlabel(ax.get_xlabel(), fontsize=self.axisfont)
			ax.set_ylabel(ax.get_ylabel(), fontsize=self.axisfont)

		if self.xticks:
			xlim = ax.get_xlim()
			ticks = np.arange(xlim[0], xlim[1], self.xticks)
			ax.set_xticks(self.ticks)

		if self.tickfont: ax.tick_params(labelsize=self.tickfont)

		if self.ltitle is not None or self.ctitle is not None or self.rtitle is not None:
			titlepad = matplotlib.rcParams['axes.titlepad']
			transOffset = matplotlib.transforms.ScaledTranslation(0., titlepad/72., ax.figure.dpi_scale_trans)
			if self.ltitle is not None:
				t = ax.text(0.0, 1.0, self.ltitle,
					ha='left', va='baseline', fontsize=self.ltitlefont, transform=ax.transAxes)
				t.set_transform(ax.transAxes + transOffset)
			if self.ctitle is not None:
				t = ax.text(0.5, 1.0, self.ctitle,
					ha='center', va='baseline', fontsize=self.ctitlefont, transform=ax.transAxes)
				t.set_transform(ax.transAxes + transOffset)
			if self.rtitle is not None:
				t = ax.text(0.0, 1.0, self.rtitle,
					ha='left', va='baseline', fontsize=self.rtitlefont, transform=ax.transAxes)
				t.set_transform(ax.transAxes + transOffset)
			ax.set_title(' ', fontsize=pickmax(self.ltitlefont, self.ctitlefont, self.rtitlefont))

		#plot_hydropathies(self.hydropathies, ax)
		if self.amphi: plot_amphipathicities(self.amphipathicities, ax)

		#ax.grid(1)
		plot.grid = self.grid
		#ax.set_xlim(0, len(self.hydropathies))
		#plot.xlim = [0, len(self.hydropathies)]
		#ax.set_ylim(-3, 3)
		#ax.set_xlabel(self.xlabel)
		plot.axeslabels = [self.xlabel, self.ylabel]
		ax.axhline(0, c='k', lw=0.5)

	def __len__(self): return self.length
	
	def parse_alignment(self, f):
		alignments = AlignIO.parse(f, 'clustal')
		for msa in alignments:
			self.length = msa.get_alignment_length()
			self.hydropathies = get_average_hydropathies(msa, window=self.window, kernel=self.kernel)
			self.amphipathicities = get_average_amphipathicities(msa, window=self.window)
			self.similarities = get_similarities(msa)
			self.tmcenters = get_tmcenters(msa)

			#what's the worst that can happen???
			return msa

def avehas(f, ax, **kwargs):#, window=19, simwindow=10, grid=False):
	alignments = AlignIO.parse(f, 'clustal')

	grid = kwargs.get('grid', False)
	if grid: ax.grid(1)

	xlabel = kwargs.get('xlabel', 'Position')
	if xlabel: ax.set_xlabel(xlabel)

	ylabel = kwargs.get('ylabel', 'Relative values')
	if ylabel: ax.set_xlabel(xlabel)

	axisfont = kwargs.get('axisfont', None)
	if axisfont:
		ax.set_xlabel(ax.get_xlabel(), fontsize=axisfont)
		ax.set_ylabel(ax.get_ylabel(), fontsize=axisfont)

	xticks = kwargs.get('xticks', None)
	if xticks:
		xlim = ax.get_xlim()
		ticks = np.arange(xlim[0], xlim[1], xticks)
		ax.set_xticks(ticks)

	tickfont = kwargs.get('tickfont', None)
	if tickfont: ax.tick_params(labelsize=tickfont)

	ltitlefont = kwargs.get('ltitlefont', None)
	ltitle = kwargs.get('ltitle', None)
	ctitlefont = kwargs.get('ctitlefont', None)
	ctitle = kwargs.get('ctitle', None)
	rtitlefont = kwargs.get('rtitlefont', None)
	rtitle = kwargs.get('rtitle', None)
	if ltitle is not None or ctitle is not None or rtitle is not None:
		titlepad = matplotlib.rcParams['axes.titlepad']
		transOffset = matplotlib.transforms.ScaledTranslation(0., titlepad/72., fig.dpi_scale_trans)
		if ltitle is not None:
			t = ax.text(0.0, 1.0, ltitle,
				ha='left', va='baseline', fontsize=ltitlefont, transform=ax.transAxes)
			t.set_transform(ax.transAxes + transOffset)
		if ctitle is not None:
			t = ax.text(0.5, 1.0, ctitle,
				ha='center', va='baseline', fontsize=ctitlefont, transform=ax.transAxes)
			t.set_transform(ax.transAxes + transOffset)
		if rtitle is not None:
			t = ax.text(0.0, 1.0, rtitle,
				ha='left', va='baseline', fontsize=rtitlefont, transform=ax.transAxes)
			t.set_transform(ax.transAxes + transOffset)
		ax.set_title(' ', fontsize=pickmax(ltitlefont, ctitlefont, rtitlefont))

	for msa in alignments:
		hydropathies = get_average_hydropathies(msa, window=window)
		amphipathicities = get_average_amphipathicities(msa, window=window)
		similarities = get_similarities(msa)
		tmcenters = get_tmcenters(msa)

		plot_hydropathies(hydropathies, ax)
		plot_amphipathicities(amphipathicities, ax)

		plot_tmcenters(tmcenters, ax, ymin=-3, ymax=-2.5)
		plot_similarities(similarities, ax, ymin=-3, ymax=-1.5)

		#ax.grid(1)
		ax.set_xlim(0, len(hydropathies))
		ax.set_ylim(-3, 3)
		ax.set_xlabel('Position')
		ax.set_ylabel('Relative values')
		ax.axhline(0, c='k', lw=0.5)

		#what's the worst that can happen???
		return msa

def test_clustalio(f, window=19):
	alignments = AlignIO.parse(f, 'clustal')
	for msa in alignments:
		hydropathies = get_average_hydropathies(msa, window=window)
		amphipathicities = get_average_amphipathicities(msa, window=window)
		similarities = get_similarities(msa)
		tmcenters = get_tmcenters(msa)

		fig = quod.plt.figure()
		fig.canvas.set_window_title('AveHAS 3')
		ax = fig.add_subplot(111)

		plot_hydropathies(hydropathies, ax)
		plot_amphipathicities(amphipathicities, ax)

		plot_tmcenters(tmcenters, ax, ymin=-3, ymax=-2.5)
		plot_similarities(similarities, ax, ymin=-3, ymax=-1.5)

		#ax.grid(1)
		ax.set_xlim(0, len(hydropathies))
		ax.set_ylim(-3, 3)
		ax.set_xlabel('Position')
		ax.set_ylabel('Relative values')
		ax.axhline(0, c='k', lw=0.5)

		#what's the worst that can happen???
		return fig, ax, msa

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	#importants
	#hydor window
	#angle

	#cosmetics
	parser.add_argument('--tight', action='store_true', help='Use tight layout (may break for plots with large titles')
	parser.add_argument('--grid', action='store_true', help='Draw a grid')
	parser.add_argument('--dpi', default=300, type=int, help='DPI in case plot gets saved')
	parser.add_argument('--width', default=None, type=float, help='Plot width in inches (or system units?) (default:dynamic)')
	parser.add_argument('--height', default=5.5, type=float, help='Plot height in inches (or system units?) (default:5.5)')
	parser.add_argument('--xlabel', help='What to label the x axis')
	parser.add_argument('--ylabel', help='What to label the y axis')
	parser.add_argument('--axisfont', help='Font size for axis labels')

	parser.add_argument('--xticks', type=int, help='X tick spacing')
	parser.add_argument('--tickfont', type=float, help='X tick label size')

	parser.add_argument('--ltitle', help='Left title')
	parser.add_argument('--ltitlefont', help='Font size for left title')
	parser.add_argument('--ctitle', help='Centered title')
	parser.add_argument('--ctitlefont', help='Font size for center title')
	parser.add_argument('--rtitle', help='Right title')
	parser.add_argument('--rtitlefont', help='Font size for right title')

	parser.add_argument('-q', action='store_true', help='Run in noninteractive mode')

	parser.add_argument('infile')
	parser.add_argument('-o', help='Where to save the plot (if anywhere) (default: nowhere)')
	args = parser.parse_args()

	interactive = True
	if args.q: interactive = False

	if interactive: 
		backend = os.environ.get('MPLBACKEND', 'Qt4Agg')
		quod.plt.switch_backend(backend)

	with open(args.infile) as f: 
		fig, ax, msa = test_clustalio(f)

	#CFG stuff
	if args.grid: ax.grid(1)

	if args.tight: fig.set_tight_layout(1)

	if args.width is None: fig.set_figwidth(msa.get_alignment_length()/50)
	else: fig.set_figwidth(args.width)
	if args.height is None: fig.set_figheight(5.5)
	else: fig.set_figheight(args.height)

	if args.xlabel: ax.set_xlabel(args.xlabel)
	if args.ylabel: ax.set_ylabel(args.ylabel)

	if args.axisfont:
		ax.set_xlabel(ax.get_xlabel(), fontsize=args.axisfont)
		ax.set_ylabel(ax.get_ylabel(), fontsize=args.axisfont)

	if args.xticks:
		ax.set_xticks(np.arange(ax.get_xlim()[0], ax.get_xlim()[1], args.xticks))

	if args.tickfont:
		ax.tick_params(labelsize=args.tickfont)

	titlepad = matplotlib.rcParams['axes.titlepad']
	transOffset = matplotlib.transforms.ScaledTranslation(0., titlepad/72., fig.dpi_scale_trans)
	if args.ltitle:
		t = ax.text(0.0, 1.0, args.ltitle,
			ha='left', va='baseline', fontsize=args.ltitlefont, transform=ax.transAxes)
		t.set_transform(ax.transAxes + transOffset)
	if args.ctitle:
		t = ax.text(0.5, 1.0, args.ctitle,
			ha='center', va='baseline', fontsize=args.ctitlefont, transform=ax.transAxes)
		t.set_transform(ax.transAxes + transOffset)
	if args.rtitle:
		t = ax.text(1.0, 1.0, args.rtitle,
			ha='right', va='baseline', fontsize=args.rtitlefont, transform=ax.transAxes)
		t.set_transform(ax.transAxes + transOffset)
	if args.ltitle or args.ctitle or args.rtitle:
		ax.set_title(' ', fontsize=pickmax(args.ltitlefont, args.ctitlefont, args.rtitlefont))

	#Save the plot (if requested)
	if args.o: fig.savefig(args.o, dpi=args.dpi)

	#fig.canvas.draw()
	quod.plt.show()

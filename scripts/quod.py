#!/usr/bin/env python
from __future__ import division, print_function, unicode_literals

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO, Seq, AlignIO
import subprocess
import argparse
import sys
import libquod

def main(*args, **kwargs):
	'''
	positional arguments: Biopython Seq objects
	kwargs:
	'''

	###

	seqlist = list(args)
	angle = kwargs.get('angle', 100)
	kernels = kwargs.get('kernels', ['flat'])
	modes = kwargs.get('modes', ['hydro'])
	window = kwargs.get('window', None)
	windows = kwargs.get('windows', [])

	###

	notms = kwargs.get('notms', [])
	loadtms = kwargs.get('loadtms', {})

	#####

	for mode in modes:
		if not libquod.validate_mode(mode): libquod.error('Mode "{}" not implemented'.format(mode))
		else: pass

	###

	skiptopo = set([int(x) for x in notms])

	loadtopo = {}
	if len(loadtms) == 1: loadtopo[0] = loadtms[0]
	elif len(loadtms):
		try: 
			int(loadtms[0])
			for i in range(0, len(loadtms), 2):
				loadtopo[loadtms[i]] = loadtms[i+1]
		except ValueError:
			for i, fn in enumerate(loadtms):
				loadtopo[i] = fn

	###

	multi = kwargs.get('multi', 'stack')

	if kwargs.get('multi') in (None, 'stack'): pass

	elif kwargs.get('multi').startswith('frag'):
		pass
		#probably do something like flagging the relevant sequences

	###

	if len(seqlist) < len(modes): 
		#multiple modes per sequence
		#in this case, repeat sequences until lengths match
		for i in range(0, len(modes) - len(seqlist)): 
			seqlist.append(seqlist[i % len(seqlist)])

			#suppress redundant HMMTOPs:
			skiptopo.add(len(seqlist)-1)

	elif len(seqlist) > len(modes):
		#multiple sequences per mode
		#in this case, repeat modes until lengths match
		for i in range(0, len(seqlist) - len(modes)): modes.append(modes[i % len(seqlist)])
		#TODO: ...unless, of course, --multi is msa or frag

	if len(seqlist) < len(kernels):
		#more kernels than modes :(
		#why, user, why?
		#just dump the extras now that redundant sequences have already been introduced
		kernels = kernels[:len(seqlist)]

	elif len(seqlist) > len(kernels):
		#more modes than kernels
		#recycle them as before
		for i in range(0, len(seqlist) - len(kernels)): kernels.append(kernels[i % len(seqlist)])

	###

	if window is None:
		for mode in modes:
			if mode.startswith('entropy'): windows.append(21)
			elif mode.startswith('psipred'): windows.append(1)
			elif mode.startswith('conform'): windows.append(41)
			elif mode.startswith('charge'): windows.append(11)
			else: windows.append(19)
	else:
		if len(window) >= len(modes):
			for mode, win in zip(modes, window): windows.append(win)
		else:
			for i, mode in enumerate(modes):
				windows.append(window[i % len(window)])

	#####

	fig = plt.figure()
	ax = fig.gca()
	ax.axhline(0, color='black', lw=1)
	xlim, ylim = None, None
	skip = 0
	for i, (seq, mode, window, kernel) in enumerate(zip(seqlist, modes, windows, kernels)):
		if skip > 0: 
			skip -= 1
			continue

		if mode.startswith('ident'): 
			skip += 1
			curve = libquod.entities.Identity.compute(seq, seqlist[i+1], kernel=kernel, window=window)

		elif multi.startswith('frag'):
			skip += 1
			curve = libquod.entities.Parser.parse_mode(mode).compute(seqlist[i+1], kernel=kernel, window=window, fragment=seq)

		elif isinstance(seq, AlignIO.MultipleSeqAlignment):
			curve = libquod.entities.Parser.parse_mode(mode).compute_msa(seq, kernel=kernel, window=window)

			#expose
			occupancy = libquod.entities.Occupancy.compute_msa(seq, kernel=kernel)
			occupancy.plot(ax=ax)

			similarity = libquod.entities.AveHASimilarity.compute_msa(seq, kernel=kernel, window=11, ymin=-3, ymax=-2)
			similarity.edgecolor = 'gray'
			similarity.plot(ax=ax)


		else:
			curve = libquod.entities.Parser.parse_mode(mode).compute(seq, kernel=kernel, window=window)

		curve.edgecolor = libquod.entities.get_darkcolor(i)
		curve.plot(ax=ax)
		xlim, ylim = libquod.update_lims(xlim, ylim, curve.get_bounding_box())

		if i not in skiptopo:
			if isinstance(seq, AlignIO.MultipleSeqAlignment):
				hmmtop = libquod.entities.HMMTOP.compute_msa(seq)
				hmmtop.mode = 'density'
				hmmtop.plot(ax=ax)

				hmmtop.facecolor = 'k'
				hmmtop.mode = 'centers'
			elif multi.startswith('frag'):
				hmmtop = libquod.entities.HMMTOP.compute(seqlist[i+1], fragment=seq)
				hmmtop.facecolor = libquod.entities.get_lightcolor(i)
			else:
				hmmtop = libquod.entities.HMMTOP.compute(seq)
				hmmtop.facecolor = libquod.entities.get_lightcolor(i)
		else: hmmtop = libquod.entities.HMMTOP([], fc=libquod.entities.get_lightcolor(i))

		if i in loadtopo: 
			with open(loadtopo[i]) as fh:
				hmmtop.spans.extend(libquod.entities.HMMTOP.parse(fh))

		hmmtop.plot(ax=ax)

	kwargs['xlim'] = xlim
	kwargs['ylim'] = ylim
	xlim, ylim = libquod.draw_other(ax, **kwargs)

	#ax.set_xlim(xlim) #Matplotlib is good enough at xlims for simple plots
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)

	if kwargs.get('grid'): ax.grid(True)

	if kwargs.get('legend'): ax.legend(True)

	if 'title' in kwargs: ax.set_title(kwargs['title'])

	if 'xlabel' in kwargs: ax.set_xlabel(kwargs['xlabel'])

	if 'ylabel' in kwargs: ax.set_ylabel(kwargs['ylabel'])

	if ylim is not None:
		ax.set_yticks(np.arange(ylim[0], ylim[1]+1, 1))

	ax.set_title(kwargs.get('title', None))

	if not kwargs.get('quiet'): plt.show()
	plt.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	####################
	ioargs = parser.add_argument_group('IO arguments')
	####################

	ioargs.add_argument('infile', nargs='*', default=['/dev/stdin'], help='Sequence file(s) to read in (default: stdin)')

	if sys.platform == 'darwin':
		ioargs.add_argument('-a', '--viewer', 
			default='open',
			help='Viewer to use for opening plots (default:open)')
	else: #TODO: support BSDs and such
		ioargs.add_argument('-a', '--viewer', 
			default='xdg-open',
			help='Viewer to use for opening plots (default:xdg-open)')

	ioargs.add_argument('-o', '--outfile', 
		help='Where to write the resulting plot')

	ioargs.add_argument('-q', '--quiet', action='store_true', 
		help='Run QUOD in quiet mode (disables interactive plot)')

	ioargs.add_argument('-s', '--sequences', action='store_true', 
		help='Interpret all positional arguments as raw sequences')

	ioargs.add_argument('-v', action='store_true',
		help='Enable verbose output, unhide warnings')

	ioargs.add_argument('--multi', default='stack',
		help='''How to handle multiple sequences. Accepts the following:

			`\033[1mstack\033[0m': print over each other,

			`frag': extract subsets from stacked pairwise alignments,

			`msa': interpret input as a multiple sequence alignment (default: stack)''')

	ioargs.add_argument('--height', 
		help='Plot height (default:5.5)')

	ioargs.add_argument('--width', default=None, 
		help='Plot width (default:dynamic)')

	####################
	coreargs = parser.add_argument_group('Core plot arguments')
	####################

	coreargs.add_argument('--angle', type=float, default=100, 
		help='Angle for amphipathicity calculations')

	coreargs.add_argument('--kernel', default='flat',
		help='Smoothing kernel to use. If using multiple, separate with "+" signs, e.g. --kernel flat+hann (\033[1mflat\033[0m, gaussian, cosine, triangular, quadratic, quartic)')

	#coreargs.add_argument('--mode', nargs='+', default=['hydro'],
	coreargs.add_argument('--mode', default='hydro',
		help='What kinds of plots to produce. Accepts multiple arguments. If plotting things on the same axes, separate with "+", e.g. --mode hydro+amphi (\033[1mhydropathy\033[0m, compositional, conformational, solventaccess, psipred, hmmtopemission)')

	coreargs.add_argument('--window', default=None,
		help='Window width for smoothing curves. Defaults to 19 for hydropathy/amphipathicity/hmmtopemission, 40 for compositional/conformational, and 1 for psipred (default: 19)')

	####################
	tmsargs = parser.add_argument_group('TMS arguments')
	####################

	tmsargs.add_argument('--add-tms', metavar='start-end[:color]', nargs='+', default=[],
		help='Add TMSs and color them. Color specifications can take the form of named colors known to Matplotlib or #-delimited hex triplets/quartets')

	#-dt

	#-et

	tmsargs.add_argument('--load-tms', metavar='(+seqid) filename', nargs='+', default=[],
		help='Load TMSs from a file formatted like HMMTOP output. Quod looks grabs as many line-final integers (and intervening whitespace) as possible, so lines filled entirely with pairs of integers are also accepted. Does not overwrite existing TMSs; use -nt to suppress TMSs for relevant sequences.')

	tmsargs.add_argument('--no-tms', metavar='(+seqid)', nargs='+', default=[],
		help='Suppress HMMTOP predictions for specific sequences')

	#-rt

	####################
	regionargs = parser.add_argument_group('Region/domain arguments')
	####################

	regionargs.add_argument('--add-region', metavar='start-end[:color][:"label"]', nargs='+',
		help='Draw regions')

	regionargs.add_argument('--region-font', type=float, default=8, 
		help='Font size for region markers in points (default: 8)')

	regionargs.add_argument('--mark', metavar='(+id):(resn1,resn2,resn3|REGEXP(:color))',
		help='Add circular markers on specified curves. Accepts comma-separated lists of residues or regular expressions')

	####################
	miscargs = parser.add_argument_group('Miscellaneous plot features')
	####################

	miscargs.add_argument('-b', '--bars', nargs='+', 
		help='Draw vertical bars at these positions')

	miscargs.add_argument('-w', metavar='start-end[:y[:+-[:text]]]', nargs='+',
		help='Surround intervals in bars with wedges pointing in.')

	miscargs.add_argument('-W', metavar='x[:scale[:y]]', nargs='+',
		help='Draw single bars with wedges')

	miscargs.add_argument('--grid', action='store_true',
		help='Draw the major grid')

	miscargs.add_argument('--legend', action='store_true',
		help='Draw the legend')

	miscargs.add_argument('--xlim', type=float, nargs=2, default=None,
		help='X limits')

	miscargs.add_argument('--xticks', type=float, default=None,
		help='X tick spacing')

	miscargs.add_argument('--ylim', type=float, nargs=2, default=None,
		help='Y limits')

	miscargs.add_argument('--yticks', type=float, default=1,
		help='Y tick spacing')

	####################
	textargs = parser.add_argument_group('Text features')
	####################

	parser.add_argument('-l', '--title',
		help='Give the figure a title')

	parser.add_argument('--xlabel',
		help='Set x axis label')

	parser.add_argument('--ylabel',
		help='Set y axis label')

	parser.add_argument('--axis-font', type=float, 
		help='Axis label size (pt)')

	parser.add_argument('--tick-font', type=float,
		help='Tick label size (pt)')
	

	args = parser.parse_args()

	seqlist = []
	for fn in args.infile:
		if args.sequences: seqlist.append(str2seq(fn))
		elif libquod.onlycaps(fn): seqlist.append(str2seq(fn))
		else: seqlist.extend(list(SeqIO.parse(fn, 'fasta'))) #TODO: Autodetect non-FASTA sequences

	kwargs = {}

	#IO args
	kwargs['outfile'] = args.outfile
	kwargs['quiet'] = args.quiet
	kwargs['multi'] = args.multi
	kwargs['height'] = args.height
	kwargs['width'] = args.width
	kwargs['multi'] = args.multi

	#core args
	kwargs['angle'] = args.angle
	kwargs['kernels'] = args.kernel.split('+')
	kwargs['modes'] = args.mode.split('+')

	#tms args
	kwargs['addtms'] = args.add_tms
	kwargs['loadtms'] = args.load_tms
	kwargs['notms'] = args.no_tms

	#region args
	kwargs['addregion'] = args.add_region
	kwargs['regionfont'] = args.region_font
	kwargs['mark'] = args.mark

	#wall args
	kwargs['walls'] = args.w

	#display args
	kwargs['grid'] = args.grid
	kwargs['legend'] = args.legend
	kwargs['title'] = args.title
	kwargs['xlabel'] = args.xlabel
	kwargs['ylabel'] = args.ylabel

	if isinstance(args.window, str): kwargs['window'] = [int(x) for x in args.window.split('+')]
	else: kwargs['window'] = args.window

	main(*seqlist, **kwargs)

	if args.outfile and not args.quiet: subprocess.call([args.viewer, args.outfile])

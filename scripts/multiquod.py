#!/usr/bin/env python
from __future__ import print_function, division, unicode_literals
import os
#if 'MPLBACKEND' not in os.environ: os.environ['MPLBACKEND'] = 'Agg'
#import matplotlib
#matplotlib.use('Qt4Agg')
#from matplotlib.figure import Figure
from collections import defaultdict
import shlex
import argparse
import re
import sys
import io
import numpy as np

import quod as quod
import matplotlib.gridspec as gridspec

import matplotlib.transforms

import warnings
try: import avehas3
except ImportError:
	avehas3 = None
	warnings.warn('Could not find AveHAS3. AveHAS subplots will be disabled.')

def error(*things):
	print(*things, file=sys.stderr)
	exit(1)

def blank_config():
	#TODO: encapsulate stuff to plot-specific classes
	cfg = {}
	cfg['n_rows'] = 0
	cfg['rowscheme'] = []
	cfg['subplots'] = set()
	cfg['sequences'] = defaultdict(lambda: '', {})
	cfg['outfile'] = None

	cfg['frag'] = {}
	cfg['width'] = 8.5
	cfg['height'] = 11
	cfg['hres'] = 100
	cfg['hgap'] = 0.5
	cfg['vgap'] = -1.5
	#cfg['outdir'] = 'magnify_plots'
	cfg['dpi'] = 300
	cfg['baselength'] = defaultdict(lambda: None, {})
	cfg['weight'] = defaultdict(lambda: 1.0, {})
	cfg['length'] = {}
	cfg['color'] = {} #defaultdict(lambda: {'line':'red', 'tms':'orange'}))
	cfg['domains'] = {}
	cfg['label'] = {}
	cfg['linestyle'] = defaultdict(lambda: '-', {})
	cfg['xticks'] = {}
	cfg['title'] = {}
	cfg['ltitle'] = {}
	cfg['ctitle'] = {}
	cfg['rtitle'] = {}
	cfg['mode'] = {}
	cfg['xlim'] = {}
	cfg['ylim'] = {}

	cfg['xlabel'] = {}
	cfg['ylabel'] = {}

	cfg['amphi'] = {}

	cfg['text'] = {}

	cfg['walls'] = {}

	cfg['rawcfg'] = ''

	return cfg

def adj_font(args, cfg):
	name = args[0]
	if 'font' not in cfg: cfg['font'] = {}
	if name not in cfg['font']: cfg['font'][name] = {}
	try: cfg['font'][name][args[1]] = float(args[2])
	except ValueError: cfg['font'][name][args[1]] = args[2]

def set_prop(args, cfg):

	if len(args) < 1: raise TypeError('Not enough arguments')
	elif len(args) == 1:
		if False: pass
		else: raise TypeError('Unrecognized unary property {}'.format(args[0]))
	#figure-wide properties
	elif len(args) == 2:
		floatprops = (
			'width',
			'height',
			'hgap',
			'margin',
			'vgap',
		)
		intprops = (
			'hres',
			'dpi',
		)
		strprops = (
		#	'mode'
		)
		if args[0] in strprops: 
			cfg[args[0]] = args[1]
		elif args[0] in floatprops: 
			try: cfg[args[0]] = float(args[1])
			except ValueError as e: raise ValueError(e)
		elif args[0] in intprops: 
			try: cfg[args[0]] = int(args[1])
			except ValueError as e: raise ValueError(e)
		else: raise TypeError('Unrecognized property {}'.format(args[0]))
	#subplot-wide properties
	elif len(args) == 3:
		floatprops = (
			'xticks',
			'weight',
			'baselength',
		)
		strprops = (
			'title',
			'ltitle',
			'ctitle',
			'rtitle',
			'xlabel',
			'ylabel',
		)

		if args[0] in strprops:
			
			if args[0] not in cfg: cfg[args[0]] = {}
			try: cfg[args[0]][args[1]] = args[2]
			except ValueError as e: raise ValueError(e)
		elif args[0] in floatprops:
			try: cfg[args[0]][args[1]] = float(args[2])
			except ValueError as e: raise ValueError(e)
	#entity-specific properties
	#elif len(args) == 4:
	elif 'linestyle' in args[0] or 'linewidth' in args[0]:
		#unstable!
		strprops = (
			'linestyle',
			'linestyle.simil',
			'linestyle.amphi',
		)
		floatprops = (
			'linewidth',
			'linewidth.simil',
			'linewidth.amphi',
		)
		if args[0] in strprops:
			try: cfg[args[0]][args[1]][int(args[2])] = args[3]
			except KeyError: cfg[args[0]][args[1]] = {int(args[2]):args[3]}
		elif args[0] in floatprops:
			if args[0] not in cfg: cfg[args[0]] = {}
			if args[1] not in cfg[args[0]]: cfg[args[0]][args[1]] = {}

			try: cfg[args[0]][args[1]][int(args[2])] = float(args[3])
			except KeyError: cfg[args[0]][args[1]] = {int(args[2]):float(args[3])}
		elif args[0] == 'font': adj_font(args[1:], cfg)

	elif len(args) >= 4:
		if args[0] == 'color':
			color(args[1:], cfg)
		elif args[0].startswith('color'):
			subentity = args[0].split('.')[1]
			color(args[1:], cfg, subentity=subentity)
		elif args[0].endswith('lim'):
			name = args[1]
			lim = [float(args[2]), float(args[3])]
			if name in cfg[args[0]]: cfg[args[0]][name] = lim
			else: cfg[args[0]] = {name:lim}
			
		#	if args[1] not in cfg[args[0]]: cfg[args[0]][args[1]] = defaultdict(lambda: {'line':'red', 'tms':'orange'}, {})
		#	for colspec in args[3:]:
		#		if colspec.startswith('line:'): 
		#			cfg[args[0]][args[1]][int(args[2])]['line'] = colspec[5:]
		#		elif colspec.startswith('tms:'): 
		#			cfg[args[0]][args[1]][int(args[2])]['tms'] = colspec[4:]
		#		else: raise ValueError('Invalid color specification {}'.format(colspec))
		elif args[0].endswith('label'): cfg[args[0]][args[1]] = args[2]
		elif args[0].endswith('font'): adj_font(args[1:], cfg)
		else: raise TypeError('Too many arguments')
	else: raise TypeError('Too many arguments')

	#elif len(args) < 2: raise TypeError('Unrecognized binary property {}'.format(args[0]))

	return cfg

def def_prop(args, cfg, frag=False, mode='quod', amphi=False):
	if len(args) < 1: raise TypeError('Not enough arguments')
	else:
		name = args[0]
		cfg['subplots'].add(name)
		cfg['sequences'][name] = []
		cfg['label'][name] = name
		cfg['linestyle'][name] = defaultdict(lambda: '-', {})
		cfg['xticks']['name'] = None
		cfg['title'][name] = ''
		cfg['frag'][name] = frag
		cfg['mode'][name] = mode
		cfg['amphi'][name] = amphi

		if mode in ('empty',):
			cfg['baselength'][name] = int(args[1])
		else:
			for i, arg in enumerate(args[1:]):
				if arg.startswith('seq:'):
					cfg['sequences'][name].append('>seq_{}\n{}'.format(i+1, arg[4:]))
				elif arg.startswith('asis:'):
					cfg['sequences'][name].append('>seq_{}\n{}'.format(i+1, arg[4:]))
				elif arg.startswith('file:'):
					with open(arg[5:]) as f:
						cfg['sequences'][name].append(f.read())
				else:
					with open(arg) as f:
						cfg['sequences'][name].append(f.read())
	return cfg

def color(cmd, cfg, subentity=None):
	parser = argparse.ArgumentParser(prog='color')
	parser.add_argument('subplot', help='Subplot name')
	#parser.add_argument('seqid', nargs='?', default=0, type=int, help='Sequence index')
	parser.add_argument('colors', metavar='[SEQID] [line:COLOR1] [tms:COLOR2]', nargs='+', help='Color specification (line:red tms:orange)')

	args = parser.parse_args(cmd)

	try: 
		seqid = int(args.colors[0])
		colors = args.colors[1:]
	except ValueError: 
		seqid = 0
		colors = args.colors

	if args.subplot not in cfg['color']: cfg['color'][args.subplot] = {}
	if seqid not in cfg['color'][args.subplot]: cfg['color'][args.subplot][seqid] = {}

	#while len(cfg['color'][args.subplot]) <= seqid: 
	#	cfg['color'][args.subplot].append({'line':'red', 'tms':'orange'})

	for colspec in colors:
		if colspec.startswith('line:'): cfg['color'][args.subplot][seqid]['line'] = colspec[5:]
		if colspec.startswith('tms:'): cfg['color'][args.subplot][seqid]['tms'] = colspec[4:]

def add_domain(args, cfg): 
	#, spans=[], yspan=[], label='', style='orange', alpha=None, pos='above', size=8, center=True)
	name = args[0]

	if len(args) < 6: raise IndexError('Incomplete domain specification (subplot, x1, x2, y1, y2, color[, label, valign, halign, size]')
	lims = [float(x) for x in args[1:5]]
	xlim = [lims[:2]]
	ylim = lims[2:4]
	ylim[1] -= ylim[0]
	color = args[5]

	#b 30 200 -2.8 -2.5 green "A domain" 8 above center
	#0  1   2    3    4     5         6  7     8      9
	if len(args) >= 7: label = args[6]
	else: label = 'unnamed_domain'

	if len(args) >= 8: font = float(args[7])
	else: font = 8.

	if len(args) >= 9: 
		valign = args[8]
		try: valign = float(valign)
		except ValueError: pass
	else: valign = 'above'

	if len(args) >= 10: halign = args[9]
	else: halign = 'center'
	#FIXME: add right-align support to quod
	if halign == 'center': halign = True
	else: halign = False

	kwargs = {'xlim':xlim, 'ylim':ylim, 'color':color, 'label':label, 'fontsize':font, 'pos':valign, 'halign':halign}

	try: cfg['domains'][name].append(kwargs)
	except KeyError: cfg['domains'][name] = [kwargs]
	return cfg
	
def add_walls(args, cfg):
	name = args[0]
	start = int(args[1])
	end = int(args[2])
	if len(args) >= 4: y = float(args[3])
	else: y = 2.

	if len(args) >= 5:
		if args[4] == '+': lim = 1
		elif args[4] == '-': lim = -1
		else: lim = None
	else: lim = None

	if len(args) >= 6: thickness = float(args[5])
	else: thickness = 1.0

	try: cfg['walls'][name].append({'start':start, 'end':end, 'y':y, 'lim':lim, 'thickness':thickness})
	except KeyError: cfg['walls'][name] = [{'start':start, 'end':end, 'y':y, 'lim':lim, 'thickness':thickness}]
	return cfg

def add_text(args, cfg):
	if args[0] != 'xy': raise NotImplementedError('Error parsing "text {}": Coordinate system {} not implemented'.format(' '.join(args), args[0]))

	i = 1

	name = args[1]
	i += 1
	cfg['text'][name] = {}

	cfg['text'][name]['halign'] = 'l'
	if args[i].startswith('l'): 
		cfg['text'][name]['halign'] = 'l'
		i += 1
	elif args[i].startswith('c'): 
		cfg['text'][name]['halign'] = 'c'
		i += 1
	elif args[i].startswith('r'): 
		cfg['text'][name]['halign'] = 'r'
		i += 1
	else:
		try: float(args[i])
		except ValueError: raise ValueError('Invalid alignment parameter: {}'.format(args[i]))

	x = float(args[i])
	i += 1
	y = float(args[i])
	i += 1
	cfg['text'][name]['pos'] = [x, y]

	try: 
		fontsize = float(args[i])
		i += 1
	except ValueError: 
		fontsize = None

	cfg['text'][name]['font'] = fontsize
	cfg['text'][name]['text'] = args[i]
	#print(cfg['text'][name])
	

def parse_config(f, cfg=None):
	cfg = blank_config() if cfg is None else cfg

	for l in f: cfg['rawcfg'] += l.decode('utf-8')
	cfg['rawcfg'] += '\n'

	f.seek(0)

	for linenum, bl in enumerate(f):
		l = bl.decode('utf-8')
		if not l.strip(): continue
		elif l.strip().startswith('#'): continue
		cmd = shlex.split(bl)
		if cmd[0] == 'addrow':
			cfg['rowscheme'].append([])
			cfg['n_rows'] += 1
			for name in cmd[1:]:
				cfg['rowscheme'][-1].append(name)
				cfg['subplots'].add(name)

		elif cmd[0] == 'set':
			try: set_prop(cmd[1:], cfg)
			#except TypeError as e: raise TypeError('l. {}: set: {}'.format(linenum+1, e))
			#except ValueError as e: raise ValueError('l. {}: set: {}'.format(linenum+1, e))
			except ZeroDivisionError: pass

		elif cmd[0] == 'def':
			parser = argparse.ArgumentParser(prog='def')
			parser.add_argument('--frag', action='store_true', help='Treat these sequences as FRAG FULL pairs')
			parser.add_argument('--mode', default='quod', help='Plotting mode (avehas, \033[1mquod\033[0m, empty)')
			parser.add_argument('--amphi', action='store_true', help='Plot amphipathicity in addition to hydropathy where applicable (currently only implemented for mode:avehas')
			parser.add_argument('label', help='Label for subplot')
			parser.add_argument('infiles', nargs='+', help='Sequences to load in')
			args = parser.parse_args(cmd[1:])

			try: def_prop([args.label] + args.infiles, cfg, frag=args.frag, mode=args.mode, amphi=args.amphi)
			#except TypeError as e: raise TypeError('l. {}: def: {}'.format(linenum+1, e))
			except ValueError as e: raise ValueError('l. {}: def: {}'.format(linenum+1, e))

		elif cmd[0] == 'addrow': cfg['rowscheme'].append(cmd[1:])

		#elif cmd[0] == 'color': 
		#	try: color(cmd[1:], cfg)
		#	except TypeError as e: raise TypeError('l. {}: color: {}'.format(linenum+1, e))
		#	except ValueError as e: raise ValueError('l. {}: color: {}'.format(linenum+1, e))

		elif cmd[0] == 'add':
			if cmd[1] == 'domain': add_domain(cmd[2:], cfg)
			elif cmd[1] == 'walls': add_walls(cmd[2:], cfg)
			else: raise TypeError('Unrecognized addition: {}'.format(cmd[1]))

		elif cmd[0] == 'tms':
			if 'tms' not in cfg: cfg['tms'] = []
			cfg['tms'].append(cmd[1:])

		elif cmd[0] == 'save':
			if len(cmd) < 2: raise TypeError('l. {}: save needs a filename'.format(linenum+1))
			cfg['outfile'] = cmd[1]

		elif cmd[0] == 'text':
			if len(cmd) < 6: raise TypeError('l. {}: text: usage: text xy SUBPLOTNAME [HALIGN] X Y [FONTSIZE] "TEXT"')
			add_text(cmd[1:], cfg)

		else: raise TypeError('l. {}: Unrecognized directive: {}'.format(linenum+1, cmd[0]))

	#check that everything in rowscheme has been properly defined
	for row in cfg['rowscheme']:
		for name in row:
			if name not in cfg['subplots']:
				raise ValueError('Subplot {} not defined'.format(name))
	return cfg

def get_titlefont(cfg, name, subname='title', default=None):
	if 'font' in cfg:
		if name in cfg['font']:
			titlefont = cfg['font'][name].get(subname, default)
		else: titlefont = None
	else: titlefont = None
	return titlefont

def plot(cfg, name, fig, ax, entities, xlabel=' ', ylabel=None, yticks=None):
	quodplot = quod.Plot(fig=fig, ax=ax)
	quodplot.width = cfg['width']
	quodplot.height = cfg['height']

	if name in cfg['xlim']: quodplot.xlim = cfg['xlim'][name]
	if name in cfg['ylim']: quodplot.ylim = cfg['ylim'][name]

	for i, e in enumerate(entities):
		if quod.is_what(e):
			e.set_style(i)
			if (name in cfg['color']) and (i in cfg['color'][name]):
				e.set_tms_color(cfg['color'][name][i]['tms'])
				e.set_curve_color(cfg['color'][name][i]['line'])
			if 'linestyle' in cfg and name in cfg['linestyle']:
				if i in cfg['linestyle'][name]: e.set_linestyle(cfg['linestyle'][name][i])
			if 'linewidth' in cfg and name in cfg['linewidth']:
				if i in cfg['linewidth'][name]: e.set_linewidth(cfg['linewidth'][name][i])
		else:
			if avehas3 is not None:
				if type(e) is avehas3.Avehas:
					e.set_style(i)
					if (name in cfg['color']) and (i in cfg['color'][name]):
						e.entities[0].style = cfg['color'][name][i]['line']
						e.entities[1].style = cfg['color'][name][i]['tms']
					#if 'linestyle' in cfg and name in cfg['linestyle']:
					#	if i in cfg['linestyle'][name]: e.entities[0].linestyle = cfg['linestyle'][name][i]
					#if 'linewidth' in cfg and name in cfg['linewidth']:
					#	if i in cfg['linewidth'][name]: e.entities[0].linewidth = cfg['linewidth'][name][i]
					#	if i in cfg['linewidth'][name]: e.entities[1].linewidth = cfg['linewidth'][name][i]
					for prop in cfg:
						#UNSTABLE!!!
						if prop.startswith('linewidth') and name in cfg[prop]:
							if '.' not in prop: 
								if i in cfg['linewidth'][name]: e.entdict['hydro'].linewidth = cfg['linewidth'][name][i]
								if i in cfg['linewidth'][name]: e.entdict['amphi'].linewidth = cfg['linewidth'][name][i]
								if i in cfg['linewidth'][name]: e.entdict['simil'].linewidth = cfg['linewidth'][name][i]
							else:
								sp = prop.split('.')
								rootprop = sp[0]
								if sp[1] == 'simil': e.linewidthsimil = cfg[prop][name][i]
								elif sp[1] == 'amphi': e.linewidthamphi = cfg[prop][name][i]
								
						if prop.startswith('linestyle') and name in cfg[prop]:
							if '.' not in prop:
								if i in cfg['linestyle'][name]: e.entities[0].linestyle = cfg['linestyle'][name][i]
							else:
								sp = prop.split('.')
								rootprop = sp[0]
								if sp[1] == 'simil': e.linestylesimil = cfg[prop][name][i]
								elif sp[1] == 'amphi': e.linestyleamphi = cfg[prop][name][i]
		quodplot.add(e)
	quodplot.render()

	ax.set_xlabel(xlabel)
	#ax.set_title(cfg['title'].get(name, ''), fontsize=get_titlefont(cfg, name))
	ax.set_title('')
	if ylabel is not None: ax.set_ylabel(ylabel)
	if yticks is not None: ax.set_yticks(yticks)

	if name in cfg['xticks']:
		xlim = ax.get_xlim()
		ax.set_xticks(np.arange(xlim[0], xlim[1], cfg['xticks'][name]))

def do_tms_stuff(entities, cfg):

	if 'tms' not in cfg: return

	def list2line(l):
		s = ''
		for x in l: s += x + ' '
		return s.strip()
	for cmd in cfg['tms']:
		if cmd[0] == 'add':
			if len(cmd) < 5: raise TypeError('tms: Not enough arguments: {}'.format(list2line(cmd)))

			parser = argparse.ArgumentParser(prog='tms add')
			parser.add_argument('name', help='Which subplot to perform on')
			parser.add_argument('--color', help='Color of TMS')
			parser.add_argument('indices', nargs='+', type=int, help='Indices to draw TMSs for')
			args = parser.parse_args(cmd[1:])

			if len(args.indices) % 2: 
				seqi = args.indices[0]
				spans = args.indices[1:]
			else:
				seqi = None
				spans = args.indices


			spans = [spans[i:i+2] for i in range(0, len(spans), 2)]

			if seqi is None:
				if args.color is None: color = 'orange'
				else: color = args.color
				entities[args.name].append(quod.HMMTOP('', nohmmtop=True, style=color))
				entities[args.name][-1].spans = spans
			else:
				si = 0
				found = False
				for ei, e in enumerate(entities[args.name]):
					try: e.entities
					except AttributeError: continue
					try: e.entities[1]
					except IndexError: continue

					if si == seqi:
						e.entities[1].spans += spans
						#if args.color is not None: e.entities[1].style = args.color
						if args.color is not None: 
							if args.name not in cfg['color']: cfg['color'][args.name] = []
							while len(cfg['color'][args.name]) <= si: 
								cfg['color'][args.name].append({'line':'red', 'tms':'orange'})
							cfg['color'][args.name][si]['tms'] = args.color
						#	cfg['color'][args.name][ei]['tms'] = args.color
						found = True
						break
					else: si += 1
				if not found:
					raise ValueError('tms: Could not find subplot {} sequence {}: {}'.format(args.name, seqi, ' '.join(cmd)))

		elif cmd[0] in ('load', 'append'):
			if len(cmd) < 4: raise TypeError('tms: Not enough arguments: {}'.format(list2line(cmd)))
			name = cmd[1]
			seqi, color = None, None
			try: seqi = int(cmd[2])
			except ValueError: color = cmd[2]
			spans = []
			for rawfn in cmd[3:]:
				if rawfn.startswith('file:'): fn = rawfn[5:]
				else: fn = rawfn
				with open(fn) as f:
					for l in f:
						raw = re.findall('(?:[0-9]+\s*)+$', l.strip())
						if not raw: continue
						else:
							for s in raw: 
								indices = [int(x) for x in s.split()]
								for i in range(len(indices)%2, len(indices), 2):
									spans.append(indices[i:i+2])
						#print(re.findall('(?:[0-9]+\s*)+$', l.strip()))
			if seqi is None:
				entities[name].append(HMMTOP('', style=color, nohmmtop=True))
				entities[name][-1].spans = spans
			else:
				si = 0
				found = False
				for ei, e in enumerate(entities[name]):
					if si == seqi:
						if cmd[0] == 'load': e.entities[1].spans = spans
						elif cmd[0] == 'append': e.entities[1].spans += spans
						else: raise Exception('Impossible exception')
						found = True
						break
					if quod.is_what(e): si += 1
				if not found: 
					raise ValueError('tms: Could not find subplot {} sequence {}: {}'.format(name, seqi, list2line(cmd)))

		elif cmd[0] == 'erase':
			if len(cmd) < 3: raise TypeError('tms: Not enough arguments: {}'.format(list2line(cmd)))
			name = cmd[1]
			spans = []

			try: 
				for i in range(3, len(cmd), 2): spans.append([int(x) for x in cmd[i:i+2]])
			except ValueError:
				if cmd[-1] == 'all': spans.append([-9e9, 9e9])
				else: raise ValueError('tms: Malformed indices: {}'.format(list2line(cmd)))

			si = 0
			seqi = int(cmd[2])
				
			found = False

			#TODO: automerge eraser targets

			for ei, e in enumerate(entities[name]):
				if si == seqi:

					popme = []
					appendme = []
					try: assert e.iswhat
					except AttributeError: continue #not a WHAT-like object
					except AssertionError: continue #not sufficiently what-like to qualify

					for i, oldspan in enumerate(e.entities[1].spans):
						for eraseme in spans:
							#case 1: oldspan is a subset of eraseme
							if (eraseme[0] <= oldspan[0] and oldspan[1] <= eraseme[1]):
								popme.append(i)
							#case 2: eraseme is a subset of oldspan
							elif (eraseme[0] >= oldspan[0] and oldspan[1] >= eraseme[1]):
								popme.append(i)
								appendme.append([oldspan[0], eraseme[0]])
								appendme.append([oldspan[1], eraseme[1]])
							#case 3: eraseme overlaps oldspan
							elif (eraseme[0] <= oldspan[0] <= eraseme[1]): 
								e.entities[1].spans[i][0] = max(oldspan[0], eraseme[1])
							#case 4: oldspan overlaps eraseme 
							if (eraseme[0] <= oldspan[1] <= eraseme[1]): 
								e.entities[1].spans[i][0] = min(oldspan[1], eraseme[0])
					popme = sorted(set(popme))
					for i in popme[::-1]: e.entities[1].spans.pop(i)

					e.entities[1].spans += appendme
					found = True
				if quod.is_what(e): si += 1
			if not found: 
				raise ValueError('tms: Could not find subplot {} sequence {}: {}'.format(name, seqi, list2line(cmd)))

		elif cmd[0] == 'delete':
			if len(cmd) < 4: raise TypeError('tms: Not enough arguments: {}'.format(list2line(cmd)))
			name = cmd[1]
			seqi = int(cmd[2])
			si = 0
			found = False

			deleteme = sorted([int(x) for x in cmd[3:]])[::-1]
			for ei, e in enumerate(entities[name]):
				if si == seqi:
					for i in deleteme: e.entities[1].spans.pop(i)
					found = True
					break
				if quod.is_what(e): si += 1
			if not found: 
				raise ValueError('tms: Could not find subplot {} sequence {}: {}'.format(name, seqi, list2line(cmd)))


	pass

def run_quod(cfg):

	entities = {}
	maxlens = {}
	for name in cfg['subplots']:
		entities[name] = []
		maxlens[name] = 0
		mode = cfg['mode'][name]
		if mode in ('hvordan', 'quod'):
			if cfg['frag'][name]:
				for i in range(0, len(cfg['sequences'][name]), 2):
					entities[name].append(quod.FragmentWhat(cfg['sequences'][name][i], cfg['sequences'][name][i+1]))
					maxlens[name] = max(maxlens[name], len(entities[name][-1]))
			else:
				for seq in cfg['sequences'][name]:
					entities[name].append(quod.What(seq))
					maxlens[name] = max(maxlens[name], len(entities[name][-1]))
		elif mode in ('avehas', ):
			if avehas3 is None: raise ImportError('Could not find avehas3.py')
			for seq in cfg['sequences'][name]:
				entities[name].append(avehas3.Avehas(io.BytesIO(seq.encode('utf-8')), amphi=cfg['amphi'][name]))
				#print(name, seq, maxlens[name])
				#for msa in avehas3.AlignIO.parse(
				#print(len(entities[name][-1]))
				maxlens[name] = max(maxlens[name], len(entities[name][-1]))
		elif mode in ('empty', ):
			spacer = Spacer(length=cfg['baselength'][name])
			maxlens[name] = cfg['baselength'][name]
			entities[name].append(spacer)

		for seq in cfg['sequences'][name]:
			#walls
			if name in cfg['walls']:
				for walldef in cfg['walls'][name]:
					entities[name].append(quod.Wall(spans=[[walldef['start'], walldef['end']]], y=walldef['y'], ylim=walldef['lim'], thickness=walldef['thickness'], wedge=walldef['thickness']))

			#domains
			if name in cfg['domains']:
				for domaindef in cfg['domains'][name]:
					entities[name].append(quod.Region(
						spans=domaindef['xlim'], 
						yspan=domaindef['ylim'],
						label=domaindef['label'],
						style=domaindef['color'],
						pos=domaindef['pos'],
						size=domaindef['fontsize'],
						center=domaindef['halign']
					))
		if name in cfg['text']:
			label = cfg['text'][name]
			entities[name].append(quod.Text(label['pos'], label['text'], label.get('halign', 'l'), label.get('font', None)))

		if cfg['baselength'][name] is None: cfg['baselength'][name] = maxlens[name]
		#cfg['length'][name] = maxlens[name]
		#cfg['weight'][name] = maxlens[name]
		cfg['length'][name] = cfg['baselength'][name] * cfg['weight'][name]

	if 'tms' in cfg: do_tms_stuff(entities, cfg)

	fig = quod.Figure()
	fig.set_tight_layout(True)

	gs = gridspec.GridSpec(len(cfg['rowscheme']), cfg['hres'])
	hgap = cfg['hgap']/cfg['width']/2
	margin = 0.03 if 'margin' not in cfg else cfg['margin']
	gs.update(
		left=margin + hgap, 
		right=cfg['width'] - margin, 
		top=cfg['height'] - margin-hgap, 
		bottom=margin, wspace=0)

	def get_limit(wt1, wtlist, offset=0): 
		x = int((offset + wt1/sum(wtlist)) * cfg['hres'])
		return max(0, min(x, cfg['hres']))

	if True:

		lims = []
		cfg['lims'] = {}
		for namelist in cfg['rowscheme']:
			row = [cfg['length'][name] for name in namelist]
			last = -hgap
			if len(namelist) == 1: 
				cfg['lims'][namelist[0]] = [0, cfg['hres']]
			else:
				bounds = [0]

				s = 0
				rawweight = [cfg['length'][name] for name in namelist]
				weight = [w/sum(rawweight) for w in rawweight]

				for i in range(len(namelist) - 1):
					bounds.append((weight[i] - hgap) * cfg['hres'])
					bounds.append((weight[i] + hgap) * cfg['hres'])

				#for i in namelist


				bounds.append(cfg['hres'])
				bounds = [int(x) for x in bounds]
				
				for i in range(0, len(namelist)):
					name = namelist[i]
					cfg['lims'][name] = [bounds[2*i], bounds[2*i + 1]]



		axdict = {}
		for r, row in enumerate(cfg['rowscheme']):
			for name in row:
				#axdict[name] = fig.add_subplot(gs[r, cfg['lims'][name][0]:cfg['lims'][name][1]])
				axdict[name] = fig.add_subplot(gs[r,cfg['lims'][name][0]:cfg['lims'][name][1]])
			

		n = 0
		name2row = {}
		stretchlabel = {}
		firstcol = set()

		for r, row in enumerate(cfg['rowscheme']):
			for c, name in enumerate(row):
				if not c: firstcol.add(name)
				name2row[name] = r
				if 'ylabel' in cfg and name in cfg['ylabel']: ylabel = cfg['ylabel'][name]
				else: ylabel = '' if c else None
				yticks = None#[] if c else None
				plot(cfg, name, fig, axdict[name], entities[name], ylabel=ylabel, yticks=yticks)

				if c: axdict[name].set_yticklabels([] * len(axdict[name].get_yticklabels()))

				try: title = cfg['title'][name]
				except KeyError: title = chr(65 + n)

				titlefont = get_titlefont(cfg, name)
				#TODO: implement rowtitle separately
				#axdict[name].set_title(cfg['title'].get(name, ''), loc='left', fontsize=titlefont)

				titlepad = quod.matplotlib.rcParams['axes.titlepad']
				transOffset = matplotlib.transforms.ScaledTranslation(0., titlepad/72., fig.dpi_scale_trans)

				if name in cfg['ctitle']:
					subtitlefont = get_titlefont(cfg, name, subname='ctitle', default=titlefont)
					t = axdict[name].text(0.5, 1.0, cfg['ctitle'].get(name, ''), 
						ha='center', va='baseline', fontsize=subtitlefont, transform=axdict[name].transAxes)
					t.set_transform(axdict[name].transAxes + transOffset)
				if name in cfg['ltitle']:
					subtitlefont = get_titlefont(cfg, name, subname='ltitle', default=titlefont)
					t = axdict[name].text(0.0, 1.0, cfg['ltitle'].get(name, ''), 
						ha='left', va='baseline', fontsize=subtitlefont, transform=axdict[name].transAxes)
					t.set_transform(axdict[name].transAxes + transOffset)
				if name in cfg['rtitle']:
					subtitlefont = get_titlefont(cfg, name, subname='rtitle', default=titlefont)
					t = axdict[name].text(1.0, 1.0, cfg['rtitle'].get(name, ''), 
						ha='right', va='baseline', fontsize=subtitlefont, transform=axdict[name].transAxes)
					t.set_transform(axdict[name].transAxes + transOffset)

				if name in cfg['title']:
					axdict[name].set_title(cfg['title'].get(name, ''), loc='left', fontsize=titlefont)
					
				n += 1

			if len(row) > 1:
				ax = fig.add_subplot(gs[r,:])
				ax.set_xticks([0])
				ax.axes.get_yaxis().set_visible(False)
				#ax.axes.get_xaxis().set_visible(False)
				ax.set_frame_on(False)
				#TODO: expose customizing xlabels per-plot/per-row
				ax.set_xlabel('Position')
				ax.tick_params(labelcolor=(1,1,1,0), color=(1,1,1,0))
				for side in ['left', 'right', 'bottom', 'top']:
					ax.spines[side].set_color((1,1,1,0))
					ax.spines[side].set_linewidth(0)

				stretchlabel[r] = ax
			elif len(row) == 1:
				#ax = fig.add_subplot(gs[r,:])
				ax = axdict[row[0]]
				ax.set_xlabel('Position')
				stretchlabel[r] = ax

		if 'xlabel' in cfg:
			for name in cfg['xlabel']:
				if name.startswith('row:'): 
					stretchlabel[name2row[name[4:]]].set_xlabel(cfg['xlabel'][name])
				else:
					axdict[name].set_xlabel(cfg['xlabel'][name])
					stretchlabel[name2row[name]].set_xlabel(' ')
		if 'ylabel' in cfg:
			for name in cfg['ylabel']:
				stretchlabel[name2row[name]].set_ylabel(cfg['ylabel'][name])

		if 'font' in cfg:
			for name in cfg['font']:
				for target in cfg['font'][name]:
					if target.endswith('ticks'):
						axdict[name].tick_params(labelsize=cfg['font'][name][target])
						if name in firstcol:
							stretchlabel[name2row[name]].tick_params(labelsize=cfg['font'][name][target])
					elif target.endswith('xaxis'):
						axdict[name].set_xlabel(axdict[name].get_xlabel(), 
							fontsize=cfg['font'][name][target])
						stretchlabel[name2row[name]].set_xlabel(stretchlabel[name2row[name]].get_xlabel(), 
							fontsize=cfg['font'][name][target])
					elif target.endswith('yaxis'):
						axdict[name].set_ylabel(axdict[name].get_ylabel(), 
							fontsize=cfg['font'][name][target])

		#gs.tight_layout(fig, pad=cfg['margin'], w_pad=cfg['hgap'], h_pad=0.0)
		#gs.tight_layout(fig, pad=0, w_pad=cfg['hgap'], h_pad=cfg['vgap'])
		gs.tight_layout(fig, pad=cfg['margin'], w_pad=cfg['hgap'], h_pad=cfg['vgap'])
		fig.savefig(cfg['outfile'], dpi=cfg['dpi'])

	else: raise NotImplementedError('Unimplemented mode: {}'.format(cfg['mode']))

class Spacer(quod.Entity):
	def __init__(self, length=200, weight=1.0): 
		self.length = length
		self.weight = weight

	def __len__(self): return self.length

	def get_bounding_box(self):
		return np.array([0, 0, self.length * self.weight, 0])

	def draw(self, plot):
		plot.ax.set_frame_on(False)
		plot.ax.set_xticks([])
		plot.ax.set_yticks([])
		plot.ax.axes.get_yaxis().set_visible(False)
		plot.ax.axes.get_xaxis().set_visible(False)
		plot.ax.tick_params(labelcolor=(1,1,1,0), color=(1,1,1,0), width=0)
		for side in ['left', 'right', 'bottom', 'top']:
			plot.ax.spines[side].set_color((1,1,1,0))
			plot.ax.spines[side].set_linewidth(0)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('infile', nargs='*', default=['/dev/stdin'], help='configuration file to load')
	args = parser.parse_args()

	cfg = blank_config()
	for fn in args.infile:
		with open(fn) as f:
			parse_config(f, cfg=cfg)
	if not cfg['outfile']: raise ValueError('No outfile defined ("save FILENAME")')
	run_quod(cfg)



from . import entities
from . import indexing

import re
import subprocess
import numpy as np
import matplotlib.cm
import sys

def info(*things):
	print('[INFO]', *things, file=sys.stderr)
def warn(*things):
	print('[WARNING]', *things, file=sys.stderr)
def error(*things):
	print('[ERROR]', *things, file=sys.stderr)
	exit(1)

def draw_other(ax, **kwargs):
	xlim = kwargs.get('xlim')
	ylim = kwargs.get('ylim')
	if 'addtms' in kwargs and kwargs['addtms']:
		for batch in kwargs['addtms']:
			args = batch.split(':')
			if len(args) == 0: raise ValueError('Invalid TMS specification: "{}"'.format(batch))

			if len(args) >= 1: spans = entities.Parser.parse_ranges(args[0])
			else: spans = []

			if len(args) >= 2: color = args[1]
			else: color = 'orange'

			hmmtop = entities.HMMTOP(spans=spans, fc=color)
			hmmtop.plot(ax)
			xlim, ylim = update_lims(xlim, ylim, hmmtop.get_bounding_box())

	if 'addregion' in kwargs and kwargs['addregion']:
		for batch in kwargs['addregion']:
			args = batch.split(':')
			if len(args) == 0: raise ValueError('Invalid region specification: "{}"'.format(batch))

			if len(args) >= 1: spans = entities.Parser.parse_ranges(args[0])
			else: spans = []

			if len(args) >= 2: color = args[1]
			else: raise ValueError('Missing color: "{}"'.format(batch))

			if len(args) >= 4:
				yspan = entities.Parser.parse_ranges(args[2])[0]
				text = args[3]
			elif len(args) >= 3: 
				yspan = [-2.75, -2.5]
				text = args[2]
			else: 
				yspan = [-2.75, -2.5]
				text = ''

			region = entities.Region(spans=spans, yspan=yspan, fc=color, text=text, ha='c', va='t', fontsize=kwargs.get('fontsize'))
			region.plot(ax)
			xlim, ylim = update_lims(xlim, ylim,region.get_bounding_box())

	if 'walls' in kwargs and kwargs['walls']:
		for batch in kwargs['walls']:
			args = batch.split(':')
			spans = entities.Parser.parse_ranges(args[0])

			#if len(args) >= 2: scale = float(args[1])
			#else: scale = 1.0

			if len(args) >= 2: y = float(args[1])
			else: y = 1.732

			if len(args) >= 3: ypos = args[2]
			else: ypos = '+-'

			if len(args) >= 4: text = args[3]
			else: text = ''

			walls = entities.Wall(spans=spans, y=y, ypos=ypos, text=text)
			xlim, ylim = update_lims(xlim, ylim, walls.get_bounding_box())
			walls.plot(ax)
		
	return xlim, ylim

def onlycaps(text):
	if re.match('^[-A-Z\s]+$', text):
		return True

def update_lims(oldxlim, oldylim, newlims=None):
	if oldxlim is None and oldylim is None:
		if newlims is None: 
			xlim = None
			ylim = None
		else:
			xlim = newlims[0]
			ylim = newlims[1]
	else:
		if newlims[0] is not None:
			xlim = [min(oldxlim[0], newlims[0][0]), max(oldxlim[1], newlims[0][1])]
		else: xlim = oldxlim

		if newlims[1] is not None:
			ylim = [min(oldylim[0], newlims[1][0]), max(oldylim[1], newlims[1][1])]
		else: ylim = oldylim

	return xlim, ylim

def validate_mode(mode):
	if mode.startswith('hydro'): return True
	elif mode.startswith('amphi'): return True
	elif mode.startswith('charge'): return True
	elif mode.startswith('entropy'): return True
	elif mode.startswith('psipred'): return True
	elif mode.startswith('conform'): return True
	elif mode.startswith('ident'): return True

	else: error('Mode "{}" not implemented'.format(mode))

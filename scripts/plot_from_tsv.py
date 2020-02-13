#!/usr/bin/env python

from __future__ import print_function, division

import matplotlib
matplotlib.use('Qt4Agg')

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import argparse
import numpy
import os
import re

def numeric(x):
	try: return int(x)
	except ValueError: pass

	try: return float(x)
	except ValueError: pass

	if '%' in x:
		try: return float(x[:-1])/100
		except ValueError: pass

	return None

class AdHocFigure(object):
	def __init__(self, fig=None, force_dims=None, x=None, y=None, z=None, xlabel=None, ylabel=None, zlabel=None, xlim=None, ylim=None, zlim=None, title=None):
		self.fig = plt.figure() if fig is None else fig
		self.dims = force_dims
		self.title = title

		self.x = x
		self.y = y
		self.z = z
		self.xlabel = xlabel
		self.ylabel = ylabel
		self.zlabel = zlabel
		self.xlim = xlim
		self.ylim = ylim
		self.zlim = zlim

		self.X = []
		self.Y = []
		self.Z = []

	def parse_tsv(self, infile):
		autolabel = True if (self.xlabel is None and self.ylabel is None and self.zlabel is None) else False
		autocol = True if (self.x is None and self.y is None and self.z is None) else False

		if self.dims is None:
			autodims = True if (self.x is None and self.y is None and self.z is None) else False
		else: autodims = False


		for l in infile:
			if not l.strip(): continue
			elif l.lstrip().startswith('#'):
				if autolabel:
					if '\t' in l: sl = l[1:].split('\t')
					#elif re.split(', *', l)
					elif re.findall(', *', l): sl = l[1:].split(', *')

					else: sl = l[1:].split('\t')
					try: self.xlabel = sl.pop(0)
					except IndexError: pass
					try: self.ylabel = sl.pop(0)
					except IndexError: pass
					try: self.zlabel = sl.pop(0)
					except IndexError: pass
			else:
				if autocol:
					sl = l.split('\t')
					done = [False, False, False]
					for rawval in sl:
						val = numeric(rawval)
						if val is None: continue

						if not done[0]:
							self.X.append(val)
							done[0] = True
						elif not done[1]:
							self.Y.append(val)
							done[1] = True
						elif not done[2]:
							self.Z.append(val)
							done[2] = True
				else:
					sl = l.split('\t')
					try:
						if self.x is not None: 
							val = numeric(sl[self.x])
							#TODO: do some validation or something
							self.X.append(val)
						if self.y is not None:
							val = numeric(sl[self.y])
							#TODO: do some validation or something
							self.Y.append(val)
						if self.z is not None: 
							val = numeric(sl[self.z])
							#TODO: do some validation or something
							self.Z.append(val)
					except IndexError: pass #maybe this shouldn't happen
		if self.xlabel is None: self.xlabel = ''
		if self.ylabel is None: self.ylabel = ''
		if self.zlabel is None: self.zlabel = ''

	def plot(self):
		if self.title is not None: self.add_title()
		truedims = sum([bool(self.X), bool(self.Y), bool(self.Z)])
		if self.dims is None: self.dims = truedims
		print(self.dims, truedims)

		if self.dims < truedims:
			if self.dims == 1:
				if truedims == 2: stuff = zip([211, 212], [self.X, self.Y], [self.xlim, self.ylim], [self.xlabel, self.ylabel])
				else: stuff = zip([221, 222, 223], [self.X, self.Y, self.Z], [self.xlim, self.ylim, self.zlim], [self.xlabel, self.ylabel, self.zlabel])
				#TODO: generalize to 4+ series
				for subplotspec, series, lim, label in stuff:
					ax = self.fig.add_subplot(subplotspec)
					hist = ax.hist(series, bins='fd')
					if lim is not None: ax.set_xlim(lim)
					if label is not None: ax.set_xlabel(label)

				ax = self.fig.add_subplot(224, projection='3d')
				scatter = ax.plot(self.X, self.Y, self.Z, lw=0, marker='o')
				if self.xlim is not None: ax.set_xlim(self.xlim)
				if self.ylim is not None: ax.set_ylim(self.ylim)
				if self.zlim is not None: ax.set_zlim(self.zlim)

				if self.xlabel is not None: ax.set_xlabel(self.xlabel)
				if self.ylabel is not None: ax.set_ylabel(self.ylabel)
				if self.zlabel is not None: ax.set_zlabel(self.zlabel)

			elif self.dims == 2:
				stuff = zip([221, 222, 223],
					[[self.X, self.Z], [self.Y, self.Z], [self.X, self.Y]],
					[[self.xlim, self.zlim], [self.ylim, self.zlim], [self.xlim, self.ylim]],
					[[self.xlabel, self.zlabel], [self.ylabel, self.zlabel], [self.xlabel, self.ylabel]],
				)
				for subplotspec, series, lims, labels in stuff:
					ax = self.fig.add_subplot(subplotspec)
					scatter = ax.plot(series[0], series[1], lw=0, marker='o')
					if lims[0] is not None: ax.set_xlim(lims[0])
					if lims[1] is not None: ax.set_ylim(lims[1])
					if labels[0] is not None: ax.set_xlabel(labels[0])
					if labels[1] is not None: ax.set_ylabel(labels[1])
				ax = self.fig.add_subplot(224, projection='3d')
				scatter = ax.plot(self.X, self.Y, self.Z, lw=0, marker='o')
				if self.xlim is not None: ax.set_xlim(self.xlim)
				if self.ylim is not None: ax.set_ylim(self.ylim)
				if self.zlim is not None: ax.set_zlim(self.zlim)

				if self.xlabel is not None: ax.set_xlabel(self.xlabel)
				if self.ylabel is not None: ax.set_ylabel(self.ylabel)
				if self.zlabel is not None: ax.set_zlabel(self.zlabel)


		#elif self.dims == truedims:
		else:
			if truedims == 1:
				ax = self.fig.add_subplot(111)
				hist = ax.hist(self.X, bins='fd')
				if self.xlim is not None: ax.set_xlim(self.xlim)
			elif truedims == 2:
				ax = self.fig.add_subplot(111)
				scatter = ax.plot(self.X, self.Y, marker='o', lw=0)
				if self.xlim is not None: ax.set_xlim(self.xlim)
				if self.ylim is not None: ax.set_ylim(self.ylim)
			elif truedims == 3:
				ax = self.fig.add_subplot(111, projection='3d')
				scatter = ax.plot(self.X, self.Y, self.Z, marker='o', lw=0)
				if self.xlim is not None: ax.set_xlim(self.xlim)
				if self.ylim is not None: ax.set_ylim(self.ylim)
				if self.zlim is not None: ax.set_zlim(self.zlim)

				if self.xlabel is not None: ax.set_xlabel(self.xlabel)
				if self.ylabel is not None: ax.set_ylabel(self.ylabel)
				if self.zlabel is not None: ax.set_zlabel(self.zlabel)
	def add_title(self, title=None):
		if title is None: title = self.title
		ax = self.fig.add_subplot(111)
		#ax.set_xticks([])
		ax.tick_params(labelcolor=(0,0,0,0), color=(0,0,0,0))
		ax.axes.get_yaxis().set_visible(False)
		ax.set_frame_on(False)
		ax.set_title(title)
	def show(self):
		self.fig.set_tight_layout(True)
		plt.show()

	def save(self, outfile, dpi=300):
		self.fig.savefig(outfile, dpi=dpi)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('-o', '--outfile', help='Where to save the figure')
	parser.add_argument('infile', help='TSV file to read data from')
	parser.add_argument('--force-dims', default=None, type=int, help='Force the resulting plot to be this many dimensions. Can be less than or equal to the number of data series.')
	parser.add_argument('--title', default=None, help='Figure title')

	parser.add_argument('-x', default=None, type=int, help='(0-indexed) column with x values')
	parser.add_argument('-y', default=None, type=int, help='(0-indexed) column with y values')
	parser.add_argument('-z', default=None, type=int, help='(0-indexed) column with z values')
	parser.add_argument('--xlabel', default=None, help='Label for x axis')
	parser.add_argument('--ylabel', default=None, help='Label for y axis')
	parser.add_argument('--zlabel', default=None, help='Label for z axis')

	parser.add_argument('--xlim', nargs=2, type=float, default=None, help='Limits for x axis')
	parser.add_argument('--ylim', nargs=2, type=float, default=None, help='Limits for y axis')
	parser.add_argument('--zlim', nargs=2, type=float, default=None, help='Limits for z axis')


	args = parser.parse_args()

	ahfig = AdHocFigure(force_dims=args.force_dims, x=args.x, y=args.y, z=args.z, xlabel=args.xlabel, ylabel=args.ylabel, zlabel=args.zlabel, xlim=args.xlim, ylim=args.ylim, zlim=args.zlim, title=args.title)
	with open(args.infile) as f: ahfig.parse_tsv(f)
	ahfig.plot()
	ahfig.show()
	if args.outfile: ahfig.save(args.outfile)

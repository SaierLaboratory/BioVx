#!/usr/bin/env python2

from __future__ import print_function
import os, sys, subprocess, re 
import gzip

import warnings
warnings.filterwarnings("ignore")

import time, hashlib
import tempfile

import numpy as np

import quod, tcblast

#import Bio.Entrez
import Bio.pairwise2, Bio.SubsMat.MatrixInfo

DEBUG = 0
VERBOSITY = 1

def warn(*msgs):
	''' prints warnings so that debug prints are more greppable '''
	for l in msgs: print('[WARNING]', l, file=sys.stderr)
def error(*msgs):
	''' prints error messages and exits with return code 1 '''
	for l in msgs: print('[ERROR]', l, file=sys.stderr)
	exit(1)
def info(*msgs):
	''' prints info messages '''
	for l in msgs: print('[INFO]', l, file=sys.stderr)

def run_pfam(indir, outdir, pfamdb):
	if not os.path.isdir(outdir): os.mkdir(outdir)
	for fn in os.listdir(indir):
		if not fn.endswith('.fa'): continue
		outfn = '{}/{}.pfam'.format(outdir, os.path.basename(os.path.splitext(fn)[0]))
		if VERBOSITY: redirect = '/dev/stderr'
		else: redirect = '/dev/null'
		cmd = ['hmmscan', '--cpu', '2', '--noali', '--cut_ga', '-o', redirect, '--domtblout', outfn, pfamdb, '{}/{}'.format(indir, fn)]
		subprocess.call(cmd)
		s = ''
		#for arg in cmd: s += arg + ' '
		#info(s)

def parse_pfam(infile, styledict, color=None, y=-2.75, size=8):
	entities = []
	spans = []
	with open(infile) as f:
		for l in f:
			if l.startswith('#'): continue
			elif not l.strip(): continue
			else:
				sl = l.strip().split()
				#label = l[181:].strip()
				#start = int(l[152:157].strip())
				#end = int(l[158:163].strip())
				#label = sl[-1]
				domain = sl[1]
				#TODO: pass clandict around to speed up clan searches
				clan = get_clan(domain)

				style = 3
				if clan.strip():
					if clan in styledict: style = styledict[clan]
					else: 
						if styledict: style = max(styledict.values())+1
						style = 3
						styledict[clan] = style
				elif domain in styledict:
					style = styledict[domain]
				else: 
					if styledict: style = max(styledict.values())+1
					else: style = 3
					styledict[domain] = style
				#print(styledict)

				start = int(sl[19])
				end = int(sl[20])

				dy = 0
				ystep = 0.7
				for span in spans:
					if (span[0] <= start <= span[1]) or (span[0] <= end <= span[1]): dy += ystep
				entities.append(quod.Region([[start, end]], [y-0.2+dy, 0.2], domain, style=style, size=size))
				#entities.append(quod.Region([[start, end]], [y-0.2+dy+ystep, 0.2], domain, style=style, size=size))
				#entities.append(quod.Region([[start, end]], [y-0.2+dy+ystep*2, 0.2], domain, style=style, size=size))
				#entities.append(quod.Region([[start, end]], [y-0.2+dy+ystep*3, 0.2], domain, style=style, size=size))
				#entities.append(quod.Region([[start, end]], [y-0.2+dy+ystep*4, 0.2], domain, style=style, size=size))
				spans.append([start, end])
	return entities, styledict

def parse_blast_tbl(arr):
	hitdata = {}
	hitdata['qid'] = arr[0]
	hitdata['sacc'] = arr[1]
	hitdata['sid'] = arr[2]
	hitdata['bitscore'] = float(arr[3])
	hitdata['expect'] = float(arr[4])
	hitdata['ident'] = float(arr[5])/100
	hitdata['qstart'] = int(arr[6])
	hitdata['qend'] = int(arr[7])
	hitdata['qlen'] = int(arr[8])
	hitdata['sstart'] = int(arr[9])
	hitdata['send'] = int(arr[10])
	hitdata['slen'] = int(arr[11])
	hitdata['ssciname'] = arr[12]
	hitdata['qseq'] = arr[13]
	hitdata['sseq'] = arr[14]

	return hitdata

def fetch(accessions, email=None, db='protein'):
	''' grabs PDBs from locally installed TCDB BLAST databases, I'm pretty sure '''
	if not accessions: return ''
	if db == 'tcdb':
		out = ''
		for acc in accessions:
			try: 
				if DEBUG: info('Running blastdbcmd')
				fa = subprocess.check_output(['blastdbcmd', '-db', 'tcdb', '-target_only', '-entry', acc])
				out += fa + '\n'
			except ValueError: raise ValueError
		return out
	else:
		if DEBUG: info('Preparing to fetch non-TCDB sequences')
		acclist = ''
		for x in accessions: acclist += ',' + x
		acclist = acclist[1:]

		try:
			if DEBUG: info('Running blastdbcmd')
			p = subprocess.Popen(['blastdbcmd', '-db', 'nr', '-entry', acclist], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			#out = re.sub('>', '\n', out) + '\n'

			if err.startswith('BLAST Database error'): raise subprocess.CalledProcessError('Database error', '1')

			remotes = ''
			for l in err.split('\n'):
				if l.strip(): 
					if 'Entry not found' in l: remotes += '%s,' % l.split()[-1]
			remotes = remotes[:-1]
			#out += subprocess.check_output(['curl', '-d', 'db=%s&id=%s&rettype=fasta&retmode=text' % (db, acclist), 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'])
			if remotes: 
				if DEBUG: info('Could not fetch some sequences locally; fetching from remote')
				out += subprocess.check_output(['curl', '-d', 'db=%s&id=%s&rettype=fasta&retmode=text' % (db, remotes), 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'])
			#out += subprocess.check_output(['curl', '-d', 'db=protein&id=Q9RBJ2&rettype=fasta&retmode=text', 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'])

			return out
		except subprocess.CalledProcessError:
			info('Could not find nr, falling back to Entrez efetch')

			if not email:
				if 'ENTREZ_EMAIL' in os.environ: email = os.environ['ENTREZ_EMAIL']
				else: 
					raise TypeError('Missing argument email')

			if DEBUG: info('Fetching from remote')
			out += subprocess.check_output(['curl', '-d', 'db=%s&id=%s&rettype=fasta&retmode=text' % (db, acclist), 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'])
			return out
	if DEBUG: info('Done fetching a batch from %s' % db)

def parse_p2report(p2report, minz=15, maxz=None, musthave=None, thispair=None):
	''' parses Protocol2 TSV reports '''

	if musthave and thispair: error('Arguments musthave and thispair are not mutually compatible')

	line = 0

	if minz == None: minz = -2**16
	if maxz == None: maxz = 2**16

	bcs = []
	alnregs = {}
	stats = {}
	for l in p2report.split('\n'):
		line += 1
		if line == 1:
			fams = [re.split('[, :]+', l)[2], re.split('[, :]+', l)[4]]
			bcs = {fams[0]:[], fams[1]:[]}
		elif line == 2: pass
		else: 
			if not l.strip(): continue
			ls = l.split('\t')
			z = float(ls[3])
			if minz <= z <= maxz: 
				#bcs.append(ls[:2])
				if musthave and ls[0] not in musthave and ls[1] not in musthave: continue
				found = 1
				if thispair:
					found = 0
					for pair in thispair:
						if ls[:2] != pair and ls[:2][::-1] != pair: continue
						else: 
							found = 1
							break
				if not found: continue
				bcs[fams[0]].append(ls[0])
				bcs[fams[1]].append(ls[1])
				try: alnregs[ls[0]][ls[1]] = (ls[6], ls[7])
				except KeyError: alnregs[ls[0]] = {ls[1]:(ls[6], ls[7])}
				try: stats[ls[0]][ls[1]] = ls[2:6]
				except KeyError: stats[ls[0]] = {ls[1]:ls[2:6]}
			
	return fams, bcs, alnregs, stats

def seek_initial(p1ds, bcs):
	''' Grabs detailed hit information '''
	hits = {}
	for fam in sorted(bcs):
		hits[fam] = {}
		for bc in sorted(bcs[fam]): hits[fam][bc] = []

	fs = {}
	fams = sorted(bcs.keys())
	for p1d in p1ds:
		if os.path.isfile(p1d):
				#fs[bc] = p1d

			#this is indeed an xor/xnor case, but that may be wrong for some directory naming schemes
			if fams[0] in p1d and fams[1] in p1d: 
				for bc in fams: fs[bc] = p1d
			elif fams[0] in p1d: fs[fams[0]] = p1d
			elif fams[1] in p1d: fs[fams[1]] = p1d
			else: 
				for bc in fams: fs[bc] = p1d
		elif os.path.isdir(p1d):
			for fam in sorted(bcs):
				if os.path.isfile('%s/%s.tbl' % (p1d, fam)): fs[fam] = '%s/%s.tbl' % (p1d, fam)

				elif os.path.isfile('%s/%s/psiblast.tbl' % (p1d, fam)): fs[fam] = '%s/%s/psiblast.tbl' % (p1d, fam)

				else: error('Could not find famXpander results table in %s' % p1d)

		else: error('Could not find p1d %s' % p1d)

	for bc in sorted(bcs):
		with open(fs[bc]) as f:
			for l in f:
				if not l.strip(): continue
				if l.lstrip().startswith('#'): continue
				if '\t' not in l: continue
				if '#' in l: continue

				ls = l.split('\t')

				hitdata = parse_blast_tbl(ls)

				try: 
					#hits[bc][ls[1]].append((float(ls[4]), ls[0], (int(ls[6]), int(ls[7])), (int(ls[9]), int(ls[10]))))
					hits[bc][ls[1]].append(hitdata)
				except KeyError: 
					#hits[bc][ls[1]] = [(float(ls[4]), ls[0], (int(ls[6]), int(ls[7])), (int(ls[9]), int(ls[10])))]
					hits[bc][ls[1]] = [hitdata]

	for fam in sorted(bcs):
		for bc in sorted(hits[fam]): 
			try: #hits[fam][bc] = sorted(hits[fam][bc])[0]
				#sortme = [(hits[fam][bc]
				besteval = None
				besthit = None
				for hit in hits[fam][bc]:
					if besteval is None: 
						besthit = hit
						besteval = besthit['expect']
					elif hit['expect'] < besteval: 
						besthit = hit
						besteval = besthit['expect']
				hits[fam][bc] = besthit
			except IndexError: error('Could not find any hits for {}/{}: Was psiblast.tbl deleted?'.format(fam, bc))
	return hits

def clean_fetch(accs, outdir, force=False, email=None):
	''' also fetches sequences but different? '''
	if DEBUG: info('Fetching %s' % accs)
	if not force:
		removeme = []
		for acc in accs:
			if os.path.isfile(outdir + '/%s.fa' % acc): removeme.append(acc)
		for acc in removeme: accs.remove(acc)

	if not os.path.isdir(outdir): os.mkdir(outdir)

	dlme = []
	tcdlme = []
	for acc in accs:
		if os.path.isfile(outdir + '/%s.fa' % acc): continue
		else: 
			if re.match('[0-9]\.[A-Z]\.[0-9]+\.', acc): tcdlme.append(acc)
			else: dlme.append(acc)

	allfaa = ''
	if dlme: 
		if VERBOSITY: info('Downloading %d sequence(s)' % len(dlme))
		allfaa += fetch(dlme, email=email)

	if tcdlme:
		if VERBOSITY: info('Loading %d TCDB sequence(s)' % len(tcdlme))
		allfaa += fetch(tcdlme, db='tcdb', email=email)

		if VERBOSITY: info('Done loading %d TCDB sequence(s)' % len(tcdlme))

	with open('%s/allseqs.faa' % outdir, 'w') as f: f.write(allfaa)

	with open('%s/allseqs.faa' % outdir) as f: 
		faa = Bio.SeqIO.parse(f, format='fasta')
		for record in faa: 
			for desc in record.description.split('>'):
				name = desc.split()[0]
				if name.count('.') < 4: name = name[:name.find('.')]
				if name.count('|') == 1: name = name.split('|')[1]

				if DEBUG > 1: info('Saving %s' % name)
				with open('%s/%s.fa' % (outdir, name), 'w') as f: f.write('>%s\n%s' % (desc, record.seq))

	fastas = {}

	#for fa in allfaa.split('\n\n'):
	#	if not fa.strip(): continue
	#	for acc in accs:
	#		if acc in fastas: pass
	#		if fa.startswith('>' + acc):
	#			fastas[acc] = fa

	#for x in sorted(fastas): 
	#	if DEBUG: info('Saving %s' % x)
	#	f = open(outdir + '/%s.fa' % x, 'w')
	#	f.write(fastas[x])
	#	f.close()

def quod_fragquod(frags, fulls, title, outfile, dpi=300, width=30, height=5.5, kernel=None):

	fig = quod.plt.figure()
	ax = fig.add_subplot(111)
	plot = quod.Plot(fig=fig, ax=ax)
	plot.add(quod.FragmentWhat(frags[0], fulls[0], style=0, kernel=kernel))
	plot.add(quod.FragmentWhat(frags[1], fulls[1], style=1, kernel=kernel))

	correl, cov = plot[0].entdict['hydro'].get_correlation(plot[1].entdict['hydro'])


	plot.width = width
	plot.height = height
	plot.render()
	plot.ax.set_title(title)
	plot.fig.savefig(outfile, dpi=dpi)


	return correl, cov

def quod_set(seqids, sequences, indir, outdir, dpi=300, force=False, bars=[], prefix='', suffix='', silent=False, pars=[], kernel=None):
	''' generates QUOD plots for batches of sequences '''
	if not os.path.isdir(outdir): os.mkdir(outdir)

	#wedges = [[[x, 2 * (0.5 - (i % 2))] for i, x in enumerate(span)] for span in bars]

	ove = lambda x: int(2 * (0.5 - (x % 2)))

	wedges = []
	for i, span in enumerate(bars):
		wedges.append([])
		if 1 <= i <= 2: y = -2
		else: y = -2
		wedges[-1].append(quod.Wall(spans=[span], y=y, ylim=[0,0.5]))

	medges = []
	for i, span in enumerate(pars):
		medges.append([])
		y = 2
		medges[-1].append(quod.Wall(spans=[span], y=y, ylim=[0.5,1]))

	domains = []
	styledict = {}
	for i, seqid in enumerate(seqids):
		#if i < 2: color = 'red'
		#else: color = 'blue'
		#domains.append(parse_pfam('{}/../pfam/{}.pfam'.format(indir, seqid), color=color))
		entities, styledict = parse_pfam('{}/../pfam/{}.pfam'.format(indir, seqid), styledict=styledict)
		domains.append(entities)

	#Draw A: barred by B
	#quod.what([sequences[seqids[0]]], force_seq=True, title=seqids[0], imgfmt='png', outdir=outdir, outfile=(seqids[0] + '_' + seqids[1] + '.png'), dpi=dpi, hide=1, entities=wedges[0]+domains[0], silent=True, width=15, height=3)

	halfwidth = 7.5
	halfheight = 2
	fig_a = quod.plt.figure()
	ax_a = fig_a.add_subplot(111)
	plot_a = quod.Plot(fig=fig_a, ax=ax_a)
	plot_a.add(quod.What(sequences[seqids[0]], style=0, kernel=kernel))
	for e in wedges[0]+domains[0]:
		plot_a.add(e)
	plot_a.width = halfwidth
	plot_a.height = halfheight
	plot_a.render()
	plot_a.ax.set_title(seqids[0])
	plot_a.fig.savefig('{}/{}_{}.png'.format(outdir, seqids[0], seqids[1]), dpi=dpi)

	#Draw B: barred by C
	#quod.what([sequences[seqids[1]]], force_seq=True, title=seqids[1], imgfmt='png', outdir=outdir, outfile=(seqids[1] + '_' + seqids[2] + '.png'), dpi=dpi, hide=1, entities=wedges[1]+medges[0]+domains[1], silent=True, width=15, height=3)
	fig_b = quod.plt.figure()
	ax_b = fig_b.add_subplot(111)
	plot_b = quod.Plot(fig=fig_b, ax=ax_b)
	plot_b.add(quod.What(sequences[seqids[1]], style=0, kernel=kernel))
	for e in wedges[1]+medges[0]+domains[1]:
		plot_b.add(e)
	plot_b.width = halfwidth
	plot_b.height = halfheight
	plot_b.render()
	plot_b.ax.set_title(seqids[1])
	plot_b.fig.savefig('{}/{}_{}.png'.format(outdir, seqids[1], seqids[2]), dpi=dpi)


	#Draw C: barred by B
	#quod.what([sequences[seqids[2]]], force_seq=True, title=seqids[2], imgfmt='png', outdir=outdir, outfile=(seqids[2] + '_' + seqids[1] + '.png'), dpi=dpi, hide=1, color=1, entities=wedges[2]+medges[1]+domains[2], silent=True, width=15, height=3)
	fig_c = quod.plt.figure()
	ax_c = fig_c.add_subplot(111)
	plot_c = quod.Plot(fig=fig_c, ax=ax_c)
	plot_c.add(quod.What(sequences[seqids[2]], style=1, kernel=kernel))
	for e in wedges[2]+medges[1]+domains[2]:
		plot_c.add(e)
	plot_c.width = halfwidth
	plot_c.height = halfheight
	plot_c.render()
	plot_c.ax.set_title(seqids[2])
	plot_c.fig.savefig('{}/{}_{}.png'.format(outdir, seqids[2], seqids[1]), dpi=dpi)

	#Draw D: barred by C
	#quod.what([sequences[seqids[3]]], force_seq=True, title=seqids[3], imgfmt='png', outdir=outdir, outfile=(seqids[3] + '_' + seqids[2] + '.png'), dpi=dpi, hide=1, color=1, entities=wedges[3]+domains[3], silent=True, width=15, height=3)
	fig_d = quod.plt.figure()
	ax_d = fig_d.add_subplot(111)
	plot_d = quod.Plot(fig=fig_d, ax=ax_d)
	plot_d.add(quod.What(sequences[seqids[3]], style=1, kernel=kernel))
	for e in wedges[3]+domains[3]:
		plot_d.add(e)
	plot_d.width = halfwidth
	plot_d.height = halfheight
	plot_d.render()
	plot_d.ax.set_title(seqids[3])
	plot_d.fig.savefig('{}/{}_{}.png'.format(outdir, seqids[3], seqids[2]), dpi=dpi)




def get_pfam(bc, prefix):
	print(prefix)
	domaindefs = []
	for acc in bc[:4]:
		with open('{}/pfam/{}.pfam'.format(prefix, acc)) as f:
			for l in f:
				if not l.strip(): continue
				elif l.startswith('#'): continue
				domaindefs.append(l.strip())
	return domaindefs

def build_html(bc, indir, blasts, outdir='hvordan_out/html', filename='test.html', lastpair=None, nextpair=None, pfam=None, alnstats=None):
	''' build an HTML report '''

	if not os.path.isdir(outdir): os.mkdir(outdir)
	if not os.path.isdir(outdir + '/assets'): os.mkdir(outdir + '/assets')

	if pfam is None: pfam = get_pfam(bc, prefix=indir)
		

	if not os.path.isfile(outdir + '/assets/openclose.js'):
		f = open(outdir + '/assets/openclose.js', 'w')
		f.write('function toggle_section(sectionid, selfid) {\n\tvar section = document.getElementById(sectionid);\n\tvar me = document.getElementById(selfid);\n\t//console.log([section, section.style.display]);\n\tif (section.style.display == \'none\') {\n\t\tsection.style.display = \'block\';\n\t\tme.innerHTML = \'-\';\n\t} else { \n\t\tsection.style.display = \'none\'; \n\t\tme.innerHTML = \'+\';\n\t}\n}')
		f.close()
	if not os.path.isfile(outdir + '/assets/nice.css'):
		f = open(outdir + '/assets/nice.css', 'w')
		#f.write('body {\n\tfont-family: sans-serif;\n\theight: 100%;\n}\ndiv {\n\tdisplay: block;\n}\ndiv.tcblast {\n\tmax-width: 1500px;\n}\ndiv.fullblast {\n\twidth: 50%;\n\tfloat: left;\n}\ndiv.tabular1 {\n\twidth: 49%;\n\tfloat: left;\n\theight: 100%;\n}\ndiv.tabular2 {\n\twidth: 49%;\n\tfloat: right;\n\theight: 100%;\n}\nimg.bluebarplot {\n\tmax-width: 100%;\n\theight: auto;\n}\n.clear { clear: both; }\n.scrollable {\n\toverflow-y: scroll;\n}\n.resizeable {\n\tresize: vertical;\n\toverflow: auto;\n\tborder: 1px solid gray;\n\tdisplay: block;\n\tpadding-bottom: 1ex;\n}\n.bluebars {\n\theight: 25vh;\n}\n.pairwise {\n\theight: 50vh;\n}\n.whatall {\n\theight: 50vh;\n}\n.whataln {\n\twidth: 100%;\n}\n#seqs {\n\tdisplay: none;\n}\n\n\n\n.summtbl {\n\tfont-family: monospace, courier;\n\tfont-size: 75%;\n}\n.oddrow {\n\tbackground-color: #d8d8d8;\n}\ntd {\n\tpadding-right: 1em;\n}\n.red {\n\tcolor: red;\n}\nimg {\n\tborder: 1pt solid black;\n}\n.monospace {\n\tfont-family: monospace;\n}')
		f.write('''body {

	font-family: sans-serif;
	height: 100%;
}
div {
	display: block;
}
div.tcblast {
	max-width: 1500px;
}
div.fullblast {
	width: 50%;
	float: left;
}
div.tabular1 {
	width: 49%;
	float: left;
	height: auto;
}
div.tabular2 {
	width: 49%;
	float: right;
	height: auto; 
}
img.bluebarplot {
	max-width: 99%;
	height: auto;
}
.clear { clear: both; }
.scrollable {
	overflow-y: scroll;
}
.resizeable {
	resize: vertical;
	overflow: auto;
	border: 1px solid gray;
	display: block;
	padding-bottom: 1ex;
}
.bluebars {
	height: 25vh;
}
.pairwise {
	height: 50vh;
}
.whatall {
	margin: 8pt;
	/*height: 50vh;*/
}
.whataln {
	width: 100%;
}
#seqs {
	display: none;
}



.summtbl {
	font-family: monospace, courier;
	font-size: 75%;
}
.oddrow {
	background-color: #d8d8d8;
}
td {
	padding-right: 1em;
}
.red {
	color: red;
}
img {
	border: 1pt solid black;
}
.monospace {
	font-family: monospace;
}
.miscinfo {
	font-size: 4pt;
}
		''')
		f.close()
	#bc := [WP_1234567890, AP_1234567890]
	title = 'HVORDAN summary: %s vs %s' % tuple(bc[1:3])

	out = '<html><head><title>%s</title>' % title
	out += '\n<link rel="stylesheet" type="text/css" href="assets/nice.css"/>'
	out += '\n<script src="assets/openclose.js"></script>'
	out += '\n</head><body>'

	out += '\n<h1>%s</h1>' % title
	if lastpair or nextpair:
		out += '\n'
		if lastpair: out += '<a href="%s_vs_%s.html">&#9664; %s vs %s</a> ' % (lastpair[1], lastpair[2], lastpair[1], lastpair[2])
		if nextpair: out += '<a href="%s_vs_%s.html">%s vs %s &#9654;</a> ' % (nextpair[1], nextpair[2], nextpair[1], nextpair[2])
		out += '<br/>'
	out += '\n<h1><button class="showhide" id="tocsh" onclick="toggle_section(\'toc\', \'tocsh\')">-</button>'
	out += '\nTable of contents</h1>'

	out += '\n<div class="toc" id="toc"> <ol> <li><a href="#summary">Summary</a></li> <li><a href="#tcsummary">TCBLAST Summary</a></li> <li><a href="#pairwise">Pairwise</a></li> <li><a href="#abcd">ABCD hydropathy plots</a></li> <li><a href="#bc">BC hydropathy plot</a></li> <li><a href="sequences">Sequences</a></li> <li><a href="domains">Domains</a></li>  </ol> </div>'

	#stats
	out += '\n<h2><button class="showhide" id="summarysh" onclick="toggle_section(\'summary\', \'summarysh\')">-</button>'
	out += '\nSummary</h2>'

	out += '\n<div class="whataln" id="summary">'
	out += '\nSS Z-score: %s<br/>' % bc[8]
	out += '\nGSAT Z-score: %s<br/>' % bc[9]
	out += '\nSubject align-length: %s<br/>' % bc[10]
	out += '\nTarget align-length: %s<br/>' % bc[11]
	out += '\n</div>'


	#bluebars
	out += '\n<h2><button class="showhide" id="tcblastsh" onclick="toggle_section(\'tcblast\', \'tcblastsh\')">-</button>'
	out += '\nTCBLAST</h2>'

	out += '\n<div class="tcblast" id="tcblast"><a name="tcsummary"><h3>TCBLAST Summary</h3></a>'
	out += '\n<div class="resizeable bluebars"><div class="scrollable tabular1">'
	out += '\n<img class="bluebarplot" src="../graphs/TCBLAST_%s.png"/>' % bc[1]
	out += '\n</div><div class="scrollable tabular2">'
	out += '\n<img class="bluebarplot" src="../graphs/TCBLAST_%s.png"/>' % bc[2]
	out += '\n</div></div>'

	#pairwise
	out += '\n<div class="clear"></div><a name="pairwise">'
	out += '\n<h3>Pairwise</h3></a><div class="resizeable pairwise"><div class="scrollable tabular1">'
	out += '\n%s' % blasts[0][1]
	out += '</div><div class="scrollable tabular2">'
	out += '\n%s' % blasts[1][1]
	out += '\n</div></div></div>'


	#abcd bc
	out += '\n<div class="clear"></div><a name="abcd">'
	out += '\n<h3><button class="showhide" id="abcdsh" onclick="toggle_section(\'abcd\', \'abcdsh\')">-</button>'
	out += '\nABCD Hydropathy plots</h3></a>'

	out += '\n<div class="whatall" id="abcd">'
	out += '\n<div class="tabular1">'
	out += '\nA<br/><img class="bluebarplot" id="plota" src="../graphs/%s_%s.png"/><br/>' % (bc[0], bc[1])
	out += '\nA-B<br/><img class="bluebarplot" id="plotab" src="../graphs/%s_vs_%s.png"/><br/>' % (bc[0], bc[1])
	out += '\nB<br/><img class="bluebarplot" id="plotb" src="../graphs/%s_%s.png"/><br/>' % (bc[1], bc[2])
	if alnstats and 'ab' in alnstats:
		out += '\nA-B stats:<br/>'
		out += '\nCorrelation: R = {:0.2f}<br/>'.format(alnstats['ab']['pearson'])
		out += '\nCoverage: {:0.0%}<br/>'.format(alnstats['ab']['coverage'])

	out += '\n</div><div class="tabular2">'

	out += '\nD<br/><img class="bluebarplot" id="plotd" src="../graphs/%s_%s.png"/><br/>' % (bc[3], bc[2])
	out += '\nD-C<br/><img class="bluebarplot" id="plotcd" src="../graphs/%s_vs_%s.png"/><br/>' % (bc[3], bc[2])
	out += '\nC<br/><img class="bluebarplot" id="plotc" src="../graphs/%s_%s.png"/><br/>' % (bc[2], bc[1])
	if alnstats and 'dc' in alnstats:
		out += '\nD-C stats:<br/>'
		out += '\nCorrelation: R = {:0.2f}<br/>'.format(alnstats['dc']['pearson'])
		out += '\nCoverage: {:0.0%}<br/>'.format(alnstats['dc']['coverage'])

	out += '\n</div></div>'

	out += '\n<div class="clear"></div><br/>'
	out += '\n<h3><button class="showhide" id="bcsh" onclick="toggle_section(\'bc\', \'bcsh\')">-</button>'
	out += '\n<a name="bc">BC hydropathy plot</a></h3>'
	out += '\n<div class="resizeable whataln" id="bc"><div class="scrollable">'
	out += '<img class="bluebarplot" id="plotbc" src="../graphs/%s_vs_%s.png"/><br/>' % (bc[1], bc[2])

	if alnstats and 'bc' in alnstats:
		out += '\nB-C stats:<br/>'
		out += '\nCorrelation: R = {:0.2f}<br/>'.format(alnstats['bc']['pearson'])
		out += '\nCoverage: {:0.0%}<br/>'.format(alnstats['bc']['coverage'])

	out += '\n</div></div>'

	#out += '\n<button class="showhide" id="tcblastsh" onclick="toggle_section(\'tcblast\', \'tcblastsh\')">Hide</button>'

	#sequences
	out += '\n<div class="clear"></div><br/><a name="sequences">'
	out += '\n<h3><button class="showhide" id="seqsh" onclick="toggle_section(\'sequences\', \'seqsh\')">-</button>'
	out += '\nSequences</h3></a>'
	out += '\n<div class="resizeable whataln monospace" id="sequences"><div class="scrollable">'
	out += ('\n%s\n%s\n%s\n%s' % tuple(bc[4:8])).replace('\n', '<br/>\n')
	out += '\n</div></div>'

	#pfam
	out += '\n<div class="clear"></div><br/><a name="domains">'
	out += '\n<h3><button class="showhide" id="domsh" onclick="toggle_section(\'domains\', \'domsh\')">-</button>'
	out += '\nDomains</h3></a>'
	out += '\n<div class="resizeable whataln monospace" id="domains"><div class="scrollable"><pre>'
	for domainstr in pfam: out += '\n{}<br/>'.format(fmt_pfam(domainstr))
	out += '\n</pre></div></div>'

	out += '\n<hr/><div class="miscinfo">'
	out += time.strftime('Generated %Y-%m-%dT%H:%M:%SZ', time.gmtime())

	out += '\n</body></html>'

	f = open(outdir + '/' + filename, 'w')
	f.write(out)
	f.close()

def get_clan(domain, clandict=None):
	if clandict:
		if domain in clandict: return clandict[domain]

	if '.' in domain: domain = domain[:domain.find('.')]

	with gzip.open(os.environ['PFAMCLANSDB']) as f:
		for l in f:
			if domain in l:
				return l.split('\t')[1]#print(l.split('\t'))
	return ''

def fmt_pfam(domainstr, clandict=None):
	domain = re.findall('PF[0-9]{5}', domainstr)[0]
	clan = get_clan(domain, clandict=clandict)
	if clan.strip():
		domainstr = re.sub('(PF[0-9]+)', r'<a href="https://pfam.xfam.org/clan/{}">{}</a>/<a href="https://pfam.xfam.org/family/\1">\1</a>'.format(clan, clan), domainstr)
	else:
		domainstr = re.sub('(PF[0-9]+)', r'NOCLAN/<a href="https://pfam.xfam.org/family/\1">\1</a>'.format(clan, clan), domainstr)
	domainstr = re.sub('([0-9]\.[A-Z]\.[0-9]+\.[0-9]+\.[0-9]+)', r'<a href="http://tcdb.org/search/result.php?tc=\1">\1</a>', domainstr)
	domainstr = re.sub('([WX]P_[0-9]+)', r'<a href="https://ncbi.nlm.nih.gov/protein/\1">\1</a>', domainstr)
	domainstr = re.sub('([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})', r'<a href="https://www.uniprot.org/uniprot/\1">\1</a>', domainstr)
	domainstr = re.sub('([A-Z]{3}[0-9]{6,})', r'<a href="https://ncbi.nlm.nih.gov/protein/\1">\1</a>', domainstr)
	domainstr = re.sub(' ([0-9][A-Z][0-9A-Z]{2})', r' <a href="https://ncbi.nlm.nih.gov/protein/\1">\1</a>', domainstr)
	return domainstr

def get_fulltrans(fams, bcs, abcd):
	''' collect A, B, C, and D into one convenient data structure '''

	pairs = zip(bcs[fams[0]], bcs[fams[1]])
	origs = [abcd[fams[0]], abcd[fams[1]]]

	fulltrans = []
	for p in pairs:
		fulltrans.append(tuple([origs[0][p[0]]['qid'], p[0], p[1], origs[1][p[1]]['qid']]))
	return fulltrans

def blastem(acc, indir, outdir, dpi=300, force=False, seqbank={}, tmcount={}, maxhits=50):
	''' generates TCBLAST plots '''
	f = open(indir + '/sequences/' + acc + '.fa')
	seq= f.read()
	f.close()

	return tcblast.til_warum(seq, outfile='%s/graphs/TCBLAST_%s.png' % (outdir, acc), title=acc, dpi=dpi, outdir='%s/blasts' % outdir, clobber=force, seqbank=seqbank, tmcount=tmcount, silent=True, maxhits=maxhits)

	#fn = outdir + '/' + filename + '.png'

	#blasts = tcblast.til_warum(seq, fn, dpi=dpi)
	#blasts = [tcblast.til_warum(l[0], args.o + '/images/' + accs[0] + '.png', dpi=args.r, html=2, outdir=args.o + '/hmmtop'), tcblast.til_warum(l[1], args.o + '/images/' + accs[1] + '.png', dpi=args.r, html=2, outdir=args.o + '/hmmtop')]

def identifind(seq1, seq2):
	''' obtains qstart, qend, sstart, send '''
	#Seq1 = Bio.Seq.Seq(seq1, Bio.Alphabet.ProteinAlphabet())
	if seq1.startswith('>'): seq1 = seq1[seq1.find('\n')+1:]
	if seq2.startswith('>'): seq2 = seq2[seq2.find('\n')+1:]
	seq1 = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', seq1.upper())
	seq2 = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', seq2.upper())

	if DEBUG: info('Starting an alignment')
	#alns = Bio.pairwise2.align.localds(seq1, seq2, Bio.SubsMat.MatrixInfo.ident, -10, -0.5)
	#out = subprocess.check_output(['ggsearch36'])
	aln = ggsearch(seq1, seq2)

	if DEBUG: info('Finished an alignment')

	subjstart = 0

	#sngap = re.findall('^-+', aln[0])
	#if sngap: sngap = len(sngap[0])
	#else: sngap = 0

	#scgap = re.findall('-+$', aln[0])
	#if scgap: scgap = len(aln[0]) - len(scgap[0]) - 1
	#else: scgap = len(aln[0])-1

	#tngap = re.findall('^-+', aln[1])
	#if tngap: tngap = len(tngap[0])
	#else: tngap = 0

	#tcgap = re.findall('-+$', aln[1])
	#if tcgap: tcgap = len(aln[1]) - len(tcgap[0]) - 1
	#else: tcgap = len(aln[1])-1

	#if sngap: 
	#	sstart = 0
	#	tstart = sngap
	#else: 
	#	sstart = tngap
	#	tstart = 0

	igap1 = re.findall('^-+', aln[0])
	igap2 = re.findall('^-+', aln[1])
	tgap1 = re.findall('-+$', aln[0])
	tgap2 = re.findall('-+$', aln[1])

	#print(seq1)
	#print(seq2)
	#print(aln[0])
	#print(aln[1])
	if igap1:
		#1 -----CYFQNCPRG
		#2 CYFQNCPRGCYFQN
		qstart = 0
		sstart = len(igap1[0])
	elif igap2:
		#1 CYFQNCPRGCYFQN
		#2 -----CYFQNCPRG
		qstart = len(igap2[0])
		sstart = 0
	else:
		#1 CYFQNCPRGCYFQN
		#2 CYFQNCPRG-----
		qstart = 0
		sstart = 0
	if tgap1:
		#1 CYFQNCPRG-----
		#2 CYFQNCPRGCYFQN
		qend = len(seq1)-1
		send = len(seq2)-1-len(tgap1[0])
	elif tgap2:
		#1 CYFQNCPRGCYFQN
		#2 CYFQNCPRG-----
		qend = len(seq1)-1-len(tgap2[1])
		send = len(seq2)-1
	else:
		#1 CYFQNCPRGCYFQN
		#2 -----CYFQNCPRG
		qend = len(seq1)-1
		send = len(seq2)-1

	return qstart+1, qend+1, sstart+1, send+1

		#I prefer 0-indexing, but pretty much everyone 1-indexes (at least for protein sequences)

def ggsearch(seq1, seq2):
	''' runs ssearch '''
	if not seq1.startswith('>'): seq1 = '>seq1\n' + seq1
	if not seq2.startswith('>'): seq2 = '>seq2\n' + seq2

	try:
		f1 = tempfile.NamedTemporaryFile(delete=False)
		f1.write(seq1)
		f1.close()

		f2 = tempfile.NamedTemporaryFile(delete=False)
		f2.write(seq2)
		f2.close()

		cmd = ['ssearch36', '-a', '-m', '3', f1.name, f2.name]
		out = subprocess.check_output(cmd).replace(' ', '-')

	finally:
		os.remove(f1.name)
		os.remove(f2.name)
	
	seqi = 0
	alns = []
	for l in out.split('\n'):
		if l.startswith('>'): seqi += 1
		if seqi:
			if not l.strip(): seqi = 0
			#elif l.startswith('>'): alns.append(l + '\n')
			#else: alns[-1] += l + '\n'
			elif l.startswith('>'): alns.append('')
			else: alns[-1] += l

	return alns


def summarize(p1d, p2d, outdir, minz=15, maxz=None, dpi=100, force=False, email=None, musthave=None, thispair=None, fams=None, maxhits=50, pfamdb='./Pfam-A.hmm', kernel=None):
	''' summarize stuff '''
	if thispair is not None:
		if len(thispair) % 2: error('Unpaired sequence found')
		else:
			truepairs = [thispair[i:i+2] for i in range(0, len(thispair), 2)]
	else: truepairs = None

	if not os.path.isdir(outdir): os.mkdir(outdir)

	if VERBOSITY: info('Reading Protocol2 report')
	try: f = open(p2d + '/report.tbl')
	except IOError:
		if os.path.isfile(p2d):
			f = open(p2d)
			warn('Opening %s as a Protocol2 results table' % p2d)
		else:
			try:
				famvfam = '%s_vs_%s' % tuple(fams)
				try: 
					f = open('%s/%s/report.tbl' % (p2d, famvfam))
					info('Could not find report.tbl in %s, falling back on family vs family subdirectory' % p2d)
				except IOError:
					try: f = open('%s/%s/%s/report.tbl' % (p2d, famvfam, famvfam))
					except IOError: error('Could not find a Protocol2 directory for %s and %s' % tuple(fams))
			except TypeError: error('Specify families if using Protocol2 root directories')
	p2report = f.read()
	f.close()

	fams, bcs, alnregs, stats = parse_p2report(p2report, minz, maxz, musthave=musthave, thispair=truepairs)

	if VERBOSITY: info('Selecting best A-B C-D pairs')
	abcd = seek_initial(p1d, bcs)
	#for k in abcd: 
	#	for j in abcd[k]: 
	#		print(k, j, abcd[k][j])
	fulltrans = get_fulltrans(fams, bcs, abcd)

	fetchme = set()
	pairstats = {}
	#for fam in abcd:
	#	for bc in abcd[fam]:
	#		fetchme.add(bc) # B|C
	#		fetchme.add(abcd[fam][bc][1]) #A|D
	#		try: pairstats[bc][abcd[fam][bc][1]] = abcd[fam][bc]
	#		except KeyError: pairstats[bc] = {abcd[fam][bc][1]:abcd[fam][bc]}
	for fam in abcd:
		for bc in abcd[fam]:
			try: pairstats[bc][abcd[fam][bc][1]] = abcd[fam][bc]
			except KeyError: pairstats[bc] = {abcd[fam][bc]['qid']:abcd[fam][bc]}
			#if 'WP_051443908' in pairstats:
			#	print('#'*80)
			#	print('WP_051443908', pairstats['WP_051443908'])

	for pair in fulltrans:
		for acc in pair: fetchme.add(acc)

	#grab all relevant sequences and store them
	if VERBOSITY: info('Retrieving %d sequence(s)' % len(fetchme))

	clean_fetch(fetchme, outdir + '/sequences', force=force, email=email)
	run_pfam(indir=(outdir + '/sequences'), outdir='{}/pfam'.format(outdir), pfamdb=pfamdb)

	if VERBOSITY: info('Done retrieving %d sequences' % len(fetchme))

	#prepare correspondences for identifind (marks B, C)
	allseqs = []
	bars = []
	seqs = {}
	pars = []
	if VERBOSITY: info('Aligning subsequences to sequences (x%d)' % len(fulltrans))

	paths = {}

	for i, pair in enumerate(fulltrans):
		[allseqs.append(x) for x in pair]

		if pair[0] not in paths: 
			paths[pair[0]] = {}
			if pair[1] not in paths[pair[0]]: paths[pair[0]] = {}
			paths[pair[0]][pair[1]] = pair[2]
		if pair[3] not in paths:
			paths[pair[3]] = {}
			if pair[2] not in paths[pair[3]]: paths[pair[3]] = {}
			paths[pair[3]][pair[2]] = pair[1]

		#bar A
		#bars.append(pairstats[pair[1]][pair[0]][3])
		#pars.append(pairstats[pair[1]][pair[0]][2])
		bars.append([pairstats[pair[1]][pair[0]]['qstart'], pairstats[pair[1]][pair[0]]['qend']])
		pars.append([pairstats[pair[1]][pair[0]]['sstart'], pairstats[pair[1]][pair[0]]['send']])
		#pars.append(pairstats[pair[1]][pair[0]][3])
		#bar B, C

		try: seqb = seqs[pair[1]]
		except KeyError:
			with open('%s/sequences/%s.fa' % (outdir, pair[1])) as f: seqb = seqs[pair[1]] = f.read()

		try: seqc = seqs[pair[2]]
		except KeyError:
			with open('%s/sequences/%s.fa' % (outdir, pair[2])) as f: seqc = seqs[pair[2]] = f.read()

		
		if DEBUG: info('Performing 2 subsequence-sequence alignments')
		bars.append(identifind(alnregs[pair[1]][pair[2]][0], seqb)[2:4])
		bars.append(identifind(alnregs[pair[1]][pair[2]][1], seqc)[2:4])

		#bar D
		#bars.append(pairstats[pair[2]][pair[3]][3])
		#pars.append(pairstats[pair[2]][pair[3]][2])
		#bars.append(pairstats[pair[2]][pair[3]][2])
		#pars.append(pairstats[pair[2]][pair[3]][3])
		bars.append([pairstats[pair[2]][pair[3]]['qstart'], pairstats[pair[2]][pair[3]]['qend']])
		pars.append([pairstats[pair[2]][pair[3]]['sstart'], pairstats[pair[2]][pair[3]]['send']])

		try: subseqs = alnregs[pair[1]][pair[2]]
		except KeyError: subseqs = alnregs[pair[2]][pair[1]]


	#make graphs for all individual full-lengthers
	if VERBOSITY: info('Generating QUOD plots')

	for x in allseqs:
		try: seqs[x]
		except KeyError: 
			with open('%s/sequences/%s.fa' % (outdir, x)) as f: seqs[x] = f.read()

	for i in range(0, len(allseqs), 4):
		quod_set(tuple(allseqs[i:i+4]), seqs, outdir + '/sequences', outdir + '/graphs/', dpi=dpi, force=force, bars=bars[i:i+4], silent=not i, pars=pars[i//2:i//2+2], kernel=kernel)

	#make graphs for all pairs of sequences
	halfwidth = 7.5
	halfheight = 2
	i = 0


	alnstats = {}

	for s1 in alnregs: 
		for s2 in alnregs[s1]: 
			#quod.what(alnregs[s1][s2], force_seq=True, labels=[s1,s2], title='%s (red) vs %s (blue)' % (s1,s2), imgfmt='png', outdir=outdir+'/graphs', outfile='%s_vs_%s.png' % (s1,s2), dpi=dpi, hide=1, width=30, height=3)
			acc0 = allseqs[i]
			acc1 = allseqs[i+1]
			acc2 = allseqs[i+2]
			acc3 = allseqs[i+3]

			seq1 = open('{}/sequences/{}.fa'.format(outdir, acc1)).read()
			seq2 = open('{}/sequences/{}.fa'.format(outdir, acc2)).read()
			correl, cov = quod_fragquod(alnregs[s1][s2], (seq1, seq2), title='{} (red) vs. {} (blue)'.format(s1, s2), outfile='{}/graphs/{}_vs_{}.png'.format(outdir, s1, s2), dpi=dpi, width=halfwidth*2, height=halfheight*2, kernel=kernel)

			alnstats['bc'] = {'pearson': correl, 'coverage': cov}
			i += 4

			seq0 = open('{}/sequences/{}.fa'.format(outdir, acc0)).read()
			frag0 = pairstats[acc1][acc0]['qseq']
			frag1 = pairstats[acc1][acc0]['sseq']
			frag2 = pairstats[acc2][acc3]['sseq']
			frag3 = pairstats[acc2][acc3]['qseq']
			seq3 = open('{}/sequences/{}.fa'.format(outdir, acc3)).read()
			#A-B
			correl, cov = quod_fragquod([frag0, frag1], [seq0, seq1], title='{} (red) vs. {} (blue)'.format(acc0, acc1), outfile='{}/graphs/{}_vs_{}.png'.format(outdir, acc0, acc1), dpi=dpi, width=halfwidth, height=halfheight, kernel=kernel)
			alnstats['ab'] = {'pearson': correl, 'coverage': cov}
			#C-D
			correl, cov = quod_fragquod([frag3, frag2], [seq3, seq2], title='{} (red) vs. {} (blue)'.format(acc3, acc2), outfile='{}/graphs/{}_vs_{}.png'.format(outdir, acc3, acc2), dpi=dpi, width=halfwidth, height=halfheight, kernel=kernel)
			alnstats['dc'] = {'pearson': correl, 'coverage': cov}

	if VERBOSITY: info('Generating TCBLAST plots')
	blasts = {}
	tmcount = {}
	seqbank = {}
	for pair in fulltrans:
		#blasts[tuple(pair)] = [blastem(pair[1], indir=outdir, outdir=outdir, dpi=dpi), blastem(pair[2], indir=outdir, outdir=outdir, dpi=dpi, force=force, seqbank=seqbank, tmcount=tmcount, maxhits=maxhits)]
		blasts[tuple(pair)] = [blastem(pair[i+1], indir=outdir, outdir=outdir, dpi=dpi, maxhits=maxhits) for i in range(2)]

	if fulltrans:
		if VERBOSITY: info('Generating %d HTML reports' % len(fulltrans))
		for i, pair in enumerate(fulltrans):
			pairseqs = []
			for seq in pair: 
				try: pairseqs.append(seqs[seq])
				except KeyError: 
					with open('%s/sequences/%s.fa' % (outdir, seq)) as f: pairseqs.append(f.read())

			if i > 0: lastpair = fulltrans[i-1]
			else: lastpair = None
			if i < (len(fulltrans)-1): nextpair = fulltrans[i+1]
			else: nextpair = None
			build_html(pair + tuple(pairseqs) + tuple(stats[pair[1]][pair[2]]), indir=outdir, blasts=blasts[tuple(pair)], outdir=(outdir + '/html'), filename='%s_vs_%s.html' % tuple(pair[1:3]), lastpair=lastpair, nextpair=nextpair, alnstats=alnstats)
	else:

		if minz is None: zmin = '-inf' 
		else: zmin = '%0.1f' % minz
		if maxz is None: zmax = '+inf'
		else: zmax = '%0.1f' % maxz

		info('Generated 0 HTML reports: No significant Protocol2 hits found with Z-scores between %s and %s' % (zmin, zmax))

if __name__ == '__main__':

	import argparse

	parser = argparse.ArgumentParser(description='HTML Visualization of Reasonable, Decent Alignment Networks')

	parser.add_argument('--p1d', metavar='PATH', default=['.'], nargs='+', help='famXpander directories or table(s) (generally psiblast.tbl). Note: Running "cut -f1-12" on psiblast.tbl will greatly improve performance, but compatibility with famXpander/9.X.99/psiblast.tbl directory structures is implemented. Directory traversal is not implemented yet.')
	parser.add_argument('--p2d', metavar='PATH', default='.', help='Protocol2 directory or results table (generally results.tbl). If using on root Protocol2 directories, -f is required.')

	parser.add_argument('-o', '--outdir', metavar='DIR', default='hvordan_out', help='output directory {default:hvordan_out}')

	parser.add_argument('-f', '--fams', metavar='FAMILY', default=None, nargs=2, help='families to inspect. Required if using --p2d on root Protocol2 directories')

	parser.add_argument('-z', '--z-min', default=15, type=int, help='minimum Z score {default:15}')
	parser.add_argument('-Z', '--z-max', default=None, type=int, help='maximum Z score {default:none}')

	parser.add_argument('-c', '--clobber', action='store_true', help='force redownloads/regenerates where applicable')
	parser.add_argument('-r', '--dpi', type=int, default=100, help='resolution of graphs {default:100}')
	parser.add_argument('-m', '--max-hits', type=int, default=10, help='how many TCBLAST hits to BLAST for. Contributes significantly to execution time for small famXpander results. {default:10}')

	if 'ENTREZ_EMAIL' in os.environ:
		parser.add_argument('-e', '--email', default=None, help='Working email in case too many requests get sent and the NCBI needs to initiate contact. Defaults to checking $ENTREZ_EMAIL if set. {current value: $ENTREZ_EMAIL == %s}' % os.environ['ENTREZ_EMAIL'])
	else: parser.add_argument('-e', '--email', default=None, help='Working email in case too many requests get sent and the NCBI needs to initiate contact. Defaults to checking $ENTREZ_EMAIL if set. {unset}')

	if 'PFAMDB' in os.environ:
		parser.add_argument('-d', '--pfamdb', default=os.environ['PFAMDB'], help='Which PFAM database to use. Defaults to checking $PFAMDB if set. (default: {})'.format(os.environ['PFAMDB']))
	else: parser.add_argument('-d', '--pfamdb', default='/ResearchData/pfam/pfamdb/Pfam-A.hmm', help='Which PFAM database to use. Defaults to checking $PFAMDB if set. (default: {})'.format('/ResearchData/pfam/pfamdb/Pfam-A.hmm'))
	parser.add_argument('--kernel', help='Use one of the named kernels instead of the old-fashioned nearest-neighbor kernel')

	parser.add_argument('-i', metavar='ACC', nargs='+', help='Operate only on pairs containing these accessions')
	parser.add_argument('-p', metavar='ACC', nargs='+', help='Operate only on these specific pairs.')

	args = parser.parse_args()

	if args.p1d == '.' and args.p2d == '.': 
		parser.print_help()
		exit()

	summarize(args.p1d, args.p2d, args.outdir, minz=args.z_min, maxz=args.z_max, dpi=args.dpi, force=args.clobber, email=args.email, musthave=args.i, thispair=args.p, fams=args.fams, maxhits=args.max_hits, pfamdb=args.pfamdb, kernel=args.kernel)

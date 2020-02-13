#! /usr/bin/env python

import os,sys
import argparse
from decimal import Decimal


'''
A small Script to take tsv files containing results form genome analysis projects and creating an overall table
'''
 
def readLines(file):

	openFile = open(file,'r')

	lines = openFile.readlines()

	return lines[1:]

def mapTCID(lines,tcids,accessions,substrates):

	tcidMap = {}

	for line in lines:

		line = line.rstrip()

		fields = line.split('\t')

		tcid = fields[2]
		tcids.add(tcid)
	
		acc = fields[1]

		hitTMS = fields[7]
			
		query = [fields[6],'%.2E' % Decimal(fields[4]),fields[0]]

		substrateData = fields[10:]


		#Overall Data
		if tcid not in accessions:

			accessions[tcid] = []

		if (acc,hitTMS) not in accessions[tcid]:

			accessions[tcid].append((acc,hitTMS))


		if tcid not in substrates:

			substrates[tcid] = substrateData

		#Genome-specifc data
		if tcid not in tcidMap:

			tcidMap[tcid] = {}


		if acc in tcidMap[tcid]:

	
			if float(query[1]) < float(tcidMap[tcid][acc][1]):


				tcidMap[tcid][acc] = query

		else:


			tcidMap[tcid][acc] = query




	return tcidMap,tcids,accessions,substrates




def printTable(genomes,tcids,tcidMaps,accessions,substrates,output):



	outputFile = open(output,'w')

	outputFile.write('#TCID\tAcc\tCategory\tSubcategory\tSubstrate\thit_tms_no\t{}\tquery_tms_no\te_value\tquery_acc\tquery_tms_no\te_value\tquery_acc\tquery_tms_no\te_value\tquery_acc\tquery_tms_no\te_value\tquery_acc\n'.format('\t'.join(genomes)))


	for tcid in tcids:


		for acc,hitTMS in accessions[tcid]:
		
			hits = []
			pos = []

			for genome in genomes:


				if tcid in tcidMaps[genome]:

					if acc in tcidMaps[genome][tcid]:

						hits.append('\t'.join(tcidMaps[genome][tcid][acc]))
						pos.append('+')
		
					else:

						hits.append('none\tnone\tnone')
						pos.append('-')
				else:

					hits.append('none\tnone\tnone')
					pos.append('-')

			substrateData = substrates[tcid]


			outputFile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(tcid,acc,'\t'.join(substrateData),hitTMS,'\t'.join(pos),'\t'.join(hits)))



	outputFile.close()


'''
Parse the command line arguments
'''

def parseArgs():

	parser = argparse.ArgumentParser(description='Given several GBLAST files, this script generates a table with the similarities and differences between the files based on TCIDs.')

	parser.add_argument('files',nargs='+',help='A space seperated list of paths to the relevant GBLAST files.')
	parser.add_argument('-out','--output', action='store',dest='output',help='REQUIRED. The path to the outfile.')
	parser.add_argument('-id',action='store',dest='id',nargs='+',help='REQUIRED. The space separated list of names that will be used in the header to identify each genome.Should be in the same order as the filenames.')

	#Parse the arguments
	args = parser.parse_args()

	gFiles = args.files
	output = args.output
	genomes = args.id

	if len(genomes) != len(gFiles):

		print('Number of files does not match the number of IDs.')
		quit()


	genomeFiles = {}

	for genome,fileName in zip(genomes,gFiles):

		if not os.path.exists(fileName):

			print('File for {} ({}) does not exist'.format(genome,fileName))
			quit()
		else:

			genomeFiles[genome] = fileName

	return genomes,genomeFiles,output	

if __name__ == "__main__":



	#Initialize tcids
	tcids = set()

	#initialize maps
	tcidMaps = {}		


	substrates = {}

	accessions = {}


	#Parse arguments
	genomes,genomeFiles,output = parseArgs()


	for genome in genomes:

		lines = readLines(genomeFiles[genome])

		tcidMaps[genome],tcids,accessions,substrates = mapTCID(lines,tcids,accessions,substrates)

	tcids = sorted(list(tcids),key=lambda x: x.split('.'))
		
	printTable(genomes,tcids,tcidMaps,accessions,substrates,output)


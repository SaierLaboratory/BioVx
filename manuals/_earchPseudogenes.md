# Documentation for script: _searchPseudogenes.py_


## Summary
Sometimes, after running [_getMultiComSystems_](https://github.com/SaierLaboratory/TCDBtools/blob/master/manuals/getMultCompSystems.md) 
there are components still missing in a transport system.  When this happens, we need to search the DNA 
sequence of the (meta)genome for the footprint of the missing component. The component could have not 
been detected due to sequencing/annotation errors, or because it became a pseudogene. This program uses 
BLASTX to compare the protein sequence of the missing components against the DNA sequence of the query (meta)genome.


## Contributors  
Yichi Zhang and Arturo Medrano-Soto


## Dependencies
The following Python module needs to be available to the script: 

1. **Python 2.X**  
You can download the Python 2 from the [official website](https://www.python.org/). The following
modules are required:  
  a) Biopython 1.75  
  b) pandas 0.24.2  
  c) requests 2.24.0  
  d) formatQuodDomain (included in the BioVx distribution)  
  e) hmmscanParser (included in the BioVx distribution)  

2. **_blast+ 2.6.0 to 2.10.0_**  
Other versions of blast may require minor adaptations. Visit the
[download site](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).  

3. **_extractFamily.pl_**  
This program is part of the [TCDBtools distribution](https://github.com/SaierLaboratory/TCDBtools/blob/master/manuals/extractFamily.md).  


## Command line options
The following options are available:

    usage: searchPseudogenes.py [-h] [-g <string input_file>] [-p <query string>]
                            [-e <evalue>] [-c] [-l] [-out <directory name>]
                            [-f <string format>] [-o]

    optional arguments:
  -h, --help            show this help message and exit
  -g <string input_file>, --genome <string input_file>
                        MANDATORY. Name or location of a complete DNA sequence
                        file in FASTA format, bzip2 format, or gzip format
  -p <query string>, --protein <query string>
                        MANDATORY. Provide one of following:
                        proteinsequence.faa - file containing protein sequence
                        of interest in FASTA format. Single or list of TCID-
                        Accessions seperated by commas - eg:
                        1.A.1.1.1-P0A334,1.A.1.10.1-P08104. Name of a file
                        containing TCID-Accessions seperated by new lines.
  -e <evalue>, --evalue <evalue>
                        the e-value threashold for blast. Default value is 1.0
  -c, --comp-based-stats
                        use composition-based statistics. Disabled by default
  -l, --low-comp-filter
                        filter low complexity regions. Disabled by default
  -out <directory name>
                        the name of output directory which contains all output
                        files. Location with name of directory is also
                        accepted. Default name is "pseudogenesResult"
                        generated in the current directory
  -f <string format>, --format <string format>
                        table - generate the table of results;alignment -
                        generate paired alignments;both - generate both table
                        and alignment output file.Default is table
  -o, --overwrite       overwrite all output files. Disabled by default

# Documentation for script: _showsubnet.py_

## Summary
This program uses the fourr network files created by script [_YutanpaNet.py_](YutanpaNet.md) 
to generate an interactive graphical layout of any subnetwork involving user-specified 
proteins and/or TC systems.


## Contributors  
Yichi Zhang and Arturo Medrano-Soto  


## Dependencies
The following Python module needs to be available to the script: 

1. **Python 3.X**  
You can download the latest version of Python 3 from the [official website](https://www.python.org/).
YutanpaNet was tested with Python 3.7.4 in MacOS Catalina 10.5.6. In addition, this program requires
access to the functions in module YutanpaNet, which is part of the distribtion.  

2. **Files generated by script _YutanpaNet.py_**  
Run first [_YutanpaNet.py_](YutanpaNet.md) and then pass the generated files to this script.  

3. **HTML template: _network_gui_templete.html.bkp_**  
This file is part of the BioVx distribution.


## Command line options
The following options are available:

    -h, --help          show this help message and exit
    
    -ft <genome feature table>, --featuretable <genome feature table>
                        MANDATORY. Path to the feature table file as
                        downloaded from NCBI genome assemblies. This file is
                        used to extract the genomic context of genes and it
                        must be compressed in gzip format.
                        
    -lt, --locus_tag    Flag. If set, the program will use locus_tag
                        accessions to identify proteins in the query genome.
                        Users are responsible for deciding whether RefSeq or
                        locus_tag accessions will be used in the analysis.
                        Default value is RefSeq.
                        
    -ad <directory path>, --address <directory path>
                        Path to the main output directory generated by program
                        getMultCompSystems.pl. This directory will be used to
                        generate links to hydropathy plots for all matches
                        between the query genome and multicomponent systems in
                        TCDB. If not specified, the edges of subnetworks will
                        not be linked to the corresponding hydropathy plots.
                        
    -l <accessions list>, --list <accessions list>
                        Mandatory. List of TC systems or protein accessions or
                        both, in the genome seperated by comma (e.g.,
                        "1.A.30.2.1,dmul_13829").
                        
    -w, --whole         Flag. If set, visualize the whole processed network.
                        This action supersedes any accessions provided with
                        the -l option.
                        
    -o <string>         The name of the subnetwork that will be generated.
                        Default value: isolated_systems.
                        
    -c, --circular      Flag. If set, the replicon structure of the query
                        genome is regarded as circular. This information is
                        used to calculate the distance between pairs of genes.
                        The default setting is: linear.
                        
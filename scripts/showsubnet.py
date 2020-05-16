#!/usr/bin/env python
# coding: utf-8

# In[1]:


import YutanpaNet as ytn
import argparse
import sys


# In[5]:


if __name__ == '__main__':
    # parse arguments
    parser = argparse.ArgumentParser()
    # create a parser to handle options and arguments. Show help info if no args
    parser.add_argument( type = str, dest = 'rn', metavar = '<raw network nodes>')
    parser.add_argument( type = str, dest = 're', metavar = '<raw network edges>')
    parser.add_argument( type = str, dest = 'pn', metavar = '<processed network nodes>')
    parser.add_argument( type = str, dest = 'pe', metavar = '<processed network edges>')
    parser.add_argument( '-ft', '--featuretable', type = str, dest = 'ft', required = True, metavar = '<featuretable.txt.gz>', help = 'MANDATORY. the file name or file name with the absolute address of the genomic feature table in gzip, .gz, compressed format.' )
    parser.add_argument( '-ad', '--address', dest = 'ad', type = str, default = '', help = 'The absolute address of a directory containing all plots generated from multi-component system anlysis. If not specified, edges of subnetworks will not direct to corresponding hydropathy plots.' )
    parser.add_argument( '-l', '--list', dest = 'list', type = str, required = True, help = 'Mandatory. list of TC systems or protein candidates in the genome seperated by comma in a string. eg ."1.A.30.2.1,1.A.30.2.5,dmul_13829"' )
    parser.add_argument( '-w', '--whole', action = 'store_true', dest = 'whole', default = False, help = 'if set, visualize the whole processed network. This action overwrites the list input from the user.' )
    parser.add_argument( '-o', dest = 'out', default = 'isolated_systems.html', help = 'the name of the subnetwork shown.' )
    parser.add_argument( '-linear', action = 'store_false', dest = 'linear', default = True, help = 'Flag. if set, replicon structure of the genome in this analysis is circular. Otherwise, use the default setting, linear.' )
    #parser.add_argument( '-p', '--physics', dest = 'physics', action = 'store_true', default = False, help = 'if set, the sub-network HTML page will show a user interface of the physics engine that users can turn on/off or change other physics engine settings.')
    
    args = parser.parse_args()
    if len(sys.argv) < 6:
        parser.print_help()
        sys.exit(1)
    G = ytn.load_network(args.rn,args.re)
    P = ytn.load_network(args.pn,args.pe)
    if args.linear == True:
        status = 'linear'
    else:
        status = 'non-linear'
    if args.whole == True:
        if 'rest' in args.pn:
            ytn.network_visualization_v2(P,args.ft,args.ad,status,args.out,True)
        else:
            ytn.network_visualization_v2(P,args.ft,args.ad,status,args.out,False)
    else:
        node_list = [str(item) for item in args.list.split(',')]
        ytn.show_subnetwork(P,G, node_list,args.ft,args.ad,status,args.out,False)


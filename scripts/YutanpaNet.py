#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
from pyvis.network import Network
import networkx as nx
import math
import os
from hmmscanParser import hmmParser
import argparse
from toolz import unique
from Bio import SeqIO
import subprocess
import more_itertools as mit
import operator
import json
import argparse
import sys
import os.path
header = ["tname", "tacc", "tlen", "qname","qacc", "qlen", "E-value", "score", "bias", "#", "of","c-Evalue", "i-Evalue", "score", "bias", "hfrom", "hto","afrom", "ato", "efrom", "eto", "acc", "description of target"]


# In[2]:


def get_rawdata(raw_table_file):
    df = pd.read_csv(raw_table_file,sep='\t',encoding = "ISO-8859-1")
# reorder cols
    cols = ['subject_accession','tcid','query_accession','subject_tms','query_tms','status','query_length','subject_length','evalue','perc_idenity','alignment_length','query_coverage','subject_coverage','neighbors']
    df = df[cols]
    df.sort_values(by=['tcid','subject_accession','evalue','neighbors'],ascending=[True,True,True,False],inplace=True)
    return df


# In[3]:




def get_mc(tcdb_faa):
    '''
    input a tcdb fasta file and return two lists. One for all single-component systems, one for multi-component systems
    '''
    systems = dict()
    with open(tcdb_faa, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            tcid = record.id.split('-')[0]
            if tcid not in systems:
                systems[tcid] = [record.id.split('-')[1]]
            else:
                systems[tcid].append(record.id.split('-')[1])

    mc = {k:systems[k] for k in systems if len(systems[k])>1}

    return mc

def get_info(qacc_list, sacc_list, address='/ResearchData/Users/amedrano/RalfRabus/MultiomponentSystems/Desulfococcus_multivorans',clan_info_file='/ResearchData/pfam/download/Pfam-A.clans.tsv.gz' ):
    # eg: add = 'plots/Q92TN4_vs_Dmul_22110/ssearch_Dmul_22110_vs_Q92TN4/files/'
    pfam_dic = dict()
    clan_df = pd.read_csv(clan_info_file,sep='\t')
    clan_df.columns = ['dacc','cacc','name','des1','des2']
    clan_dic = dict(zip(clan_df['dacc'].values, clan_df['cacc'].values))
    pairs = zip(qacc_list,sacc_list)
    for p in pairs:
        if str(p[1])=='nan':
            continue
        hmmscan = os.path.join(address,'plots/{}_vs_{}/ssearch_{}_vs_{}/files/hmmscan.out'.format(p[0],p[1],p[1],p[0]))
        hmm = hmmParser(hmmscan)
        df=pd.DataFrame(hmm.matrix,columns=header)
        for index, row in df.iterrows():
            qname = row['qname'].split('-')[-1]
            if qname not in pfam_dic:
                pfam_dic[qname] = {}
            pacc = row['tacc'].split('.')[0]
            if pacc not in pfam_dic[qname]:
                pfam_dic[qname][pacc] = {} 
                pfam_dic[qname][pacc]['cood'] = [[int(row['efrom']),int(row['eto'])]]
                pfam_dic[qname][pacc]['clan'] = clan_dic[pacc]
            else:
                pfam_dic[qname][pacc]['cood'].append([int(row['efrom']),int(row['eto'])])
            pfam_dic[qname][pacc]['cood'] = list(map(list, unique(map(tuple, pfam_dic[qname][pacc]['cood']))))
            pfam_dic[qname][pacc]['cood'].sort()
            
    return pfam_dic


def get_cood(qacc, sacc, address):
    ssearch = os.path.join(address,'plots/{}_vs_{}/ssearch_{}_vs_{}/files/ssearch36.out'.format(qacc,sacc,sacc,qacc))
    goal = ''
    with open(ssearch,"r") as f:
        for line in f:
            if line.startswith("Smith-Waterman"):
                goal = line
                break
    coods = goal.split('overlap')[1]
    s_coods = coods.split(':')[0].split('(')[1]
    q_coods = coods.split(':')[1].split(')')[0]
    s_cood = [int(s_coods.split('-')[0]),int(s_coods.split('-')[1])]
    q_cood = [int(q_coods.split('-')[0]),int(q_coods.split('-')[1])]           
    return q_cood, s_cood

def alignment_has_domain(qacc,sacc,pfam_dic,address ):
    q_cood,s_cood = get_cood(qacc, sacc,address) # get the alignment coordinates
    # check if one of the protein does not has any domain predicted
    if not qacc in list(pfam_dic.keys()) or not sacc in list(pfam_dic.keys()):
        return [False,False]

    q_domains = list(pfam_dic[qacc].keys())
    s_domains = list(pfam_dic[sacc].keys())
    # check if domain exists in both alignments and if there is a common domain
    qdwa = []
    sdwa = []
    for qdacc in q_domains:
        for loc in pfam_dic[qacc][qdacc]['cood']:
            if (loc[0] <= q_cood[1] and loc[0] >= q_cood[0]) or (loc[1] <= q_cood[1] and loc[1] >= q_cood[0]) or (loc[1] > q_cood[1] and loc[0] < q_cood[0]):
                qdwa.append(qdacc)
    for sdacc in s_domains:
        for loc in pfam_dic[sacc][sdacc]['cood']:
            if (loc[0] <= s_cood[1] and loc[0] >= s_cood[0]) or (loc[1] <= s_cood[1] and loc[1] >= s_cood[0]) or (loc[1] > s_cood[1] and loc[0] < s_cood[0]):
                sdwa.append(sdacc)
      
    if len(qdwa) == 0 or len(sdwa) == 0:
        return [False,False] # at least one protein has no domian info within the alignment
    # check intersection
    intersact = list(set(qdwa) & set(sdwa))
    if len(intersact) == 0:
        return [True,False] # alignment has no common domain
    else:
        return [True,True] # alignment has common domains
            
def domain_extension(qacc,sacc,pfam_dic,address):
    q_cood,s_cood = get_cood(qacc, sacc,address) # get the alignment coordinates
    qde,sde = get_cood(qacc, sacc,address)
    q_domains = list(pfam_dic[qacc].keys())
    s_domains = list(pfam_dic[sacc].keys())
    # check if domain exists in both alignments and if there is a common domain
    qdwa = []
    sdwa = []
    for qdacc in q_domains:
        for loc in pfam_dic[qacc][qdacc]['cood']:
            if (loc[0] <= q_cood[1] and loc[0] >= q_cood[0]) or (loc[1] <= q_cood[1] and loc[1] >= q_cood[0]) or (loc[1] > q_cood[1] and loc[0] < q_cood[0]):
                qdwa.append(qdacc)
    for sdacc in s_domains:
        for loc in pfam_dic[sacc][sdacc]['cood']:
            if (loc[0] <= s_cood[1] and loc[0] >= s_cood[0]) or (loc[1] <= s_cood[1] and loc[1] >= s_cood[0]) or (loc[1] > s_cood[1] and loc[0] < s_cood[0]):
                sdwa.append(sdacc)
    intersact = list(set(qdwa) & set(sdwa))
    for i in intersact:
        for coods in pfam_dic[qacc][i]['cood']:
            if coods[1] < q_cood[0]:
                continue
            if coods[0] > q_cood[1]:
                continue
            if coods[0] < q_cood[0] and coods[1] < q_cood[1] and coods[0] < qde[0]:
                qde[0] = coods[0]
            if coods[0] < q_cood[1] and coods[1] > q_cood[1] and coods[1] > qde[1]:
                qde[1] = coods[1]
        for coods in pfam_dic[sacc][i]['cood']:
            if coods[1] < s_cood[0]:
                continue
            if coods[0] > s_cood[1]:
                continue
            if coods[0] < s_cood[0] and coods[1] < s_cood[1] and coods[0] < sde[0]:
                sde[0] = coods[0]
            if coods[0] < s_cood[1] and coods[1] > s_cood[1] and coods[1] > sde[1]:
                sde[1] = coods[1]
    return qde,sde

def count_complete_systems(G,nodes): # G is the target network for checking completeness
    # number of complete systems
    systems = [x for x,y in G.nodes(data=True) if y['attr_dic']['tp'] == 'system']
    num = 0
    for s in systems:
        situation = check_complete(nodes,G,s)
        if situation[0] == True:
            num = num + 1
    return num

def get_gene_feature(gene_feature_file):
    stdoutdata = subprocess.getoutput("zgrep CDS {} | cut -f 8,9,10,17".format(gene_feature_file))
    #a = subprocess.Popen("zgrep gene {} | cut -f 8,9,10,17".format(gene_feature_file), stdout=subprocess.PIPE)
    data = stdoutdata.split()
    #print(data)
    gene_feature = []
    for i in range(int(len(data)/4)):
        index = i*4
        gene_feature.append([data[index+3],data[index+2],int(data[index]),int(data[index+1])])
        # gene_accession/locus_tag, strand, start, end

    gene_feature.sort(key=lambda x: x[2]) # sort by the left coordinate
    
    # return a dic
    gene_dic = {}
    for gene in gene_feature:
        gene_dic[gene[0]] = {'strand':gene[1],'start':gene[2],'end':gene[3]}
    return gene_dic

def get_strand(gene_dic, gacc):
    return gene_dic[gacc]['strand']

def get_genetic_distance(gene_dic, gacc_list, condition): # returns genetic distance between genes according to different situations based on biological fact of that genome
    # genetic distance is represented as the number of genes apart from any two genes
    # len(gacc_list) must be greater than 1
    if len(gacc_list) < 2:
        raise ValueError('Need more gene accessions!')
    gene_list = list(gene_dic.keys())
    gene_distance = {}
    for gacc in gacc_list:
        gene_distance[gacc] = {}
        for g in gacc_list:
            if condition == 'linear':
                gene_distance[gacc][g] = abs(gene_list.index(gacc)-gene_list.index(g))
            else:
                dis = abs(gene_list.index(gacc)-gene_list.index(g))
                gene_distance[gacc][g] = min(dis, len(gene_list)-dis)
    return gene_distance
#eg: et_genetic_distance(gene_dic, gacc_list=['Dmul_00010','Dmul_02460','Dmul_04350'], condition='linear')
  # return: {'Dmul_00010': {'Dmul_00010': 0, 'Dmul_02460': 252, 'Dmul_04350': 440},
 #          'Dmul_02460': {'Dmul_00010': 252, 'Dmul_02460': 0, 'Dmul_04350': 188},
 #.         'Dmul_04350': {'Dmul_00010': 440, 'Dmul_02460': 188, 'Dmul_04350': 0}}
                


# In[4]:


def generate_weight(qacc,sacc,evalue,qcov,scov,pfam_dic,address):
    # domin_info are dictionaries that key is domain_clan, value are coordinates,
    # domain orders should be sorted based on the first coordinate
    # normalized evalue is the -log(eval,10) and assigned with tiers. 
    weight = 0 # possibility that a protein genome travels to a component system
    if evalue == 0:
        log_eval = 100
    else:
        log_eval = -1 * math.log(evalue,10)
    normalized_eval = 0
    if log_eval >= 30 :
        normalized_eval = 1
    elif log_eval >= 11:
        normalized_eval = 0.8
    elif log_eval >= 7:
        normalized_eval = 0.6
    else:
        normalized_eval = 0.4
    situations = alignment_has_domain(qacc,sacc,pfam_dic,address)
    if situations[0] == False:
        weight = 0.5*normalized_eval + 0.5* max(scov,qcov)/100 
        # I expected one of them has a high coverage in two types of fusion
        # number of difference of tms within the alignment range will be considered, but in a limited affect
    else:
        if situations[1] == True:
            weight = 0.25*normalized_eval + 0.25* max(scov,qcov)/100 +0.5
            # image a protein like this: subject plot:  --  --  --  --  -- 
            #------------------------------------------| | | | | | | | | |, 5 tms in total
            # query in genome hit the left part without any tms:
            #     -------------------------------------
            # coverage of query: ~90%, coverage of subject: ~50%, totoal TMS difference: 5, normalized-eval: 1
            # domain overlap ? YES
            # according to the equation stated, the weight/possibility/confidence will be really close to 1 
            # since only the larger percent coverage will be used 
            # during the graph search step, this protein still has a high chance to travel to this component, 
            # and we expect a second protein hits a membrane part
            
        else:
            print('punish_case: '+qacc + ', '+sacc)
            weight = 0.25*normalized_eval + 0.25* max(scov,qcov)/100 
    return int((1-weight)*100)


# In[28]:


def add_system(nodes,edges,T_com,tcid,S): #returns a modified subnetwork by adding a new system into the netowrk #nodes/edges is a dic of all nodes generated by the raw_network
    if tcid not in S:
        S.add_node(tcid, attr_dic=nodes[tcid])
    for cpt in list(T_com.predecessors(tcid)):
        if cpt not in S:
            S.add_node(cpt,attr_dic=nodes[cpt])
            S.add_edge(cpt,tcid,attr_dic=edges[(cpt,tcid)])
        for cand in list(T_com.predecessors(cpt)):
            if cand not in S:
                S.add_node(cand, attr_dic=nodes[cand])
            S.add_edge(cand,cpt,attr_dic=edges[(cand,cpt)])
def add_system_v2(nodes,edges,T_com,tcid,S): #returns a modified subnetwork by adding a new system into the netowrk #nodes/edges is a dic of all nodes generated by the raw_network
    if tcid not in S:
        S.add_node(tcid, attr_dic=nodes[tcid]['attr_dic'])
    for cpt in list(T_com.predecessors(tcid)):
        if cpt not in S:
            S.add_node(cpt,attr_dic=nodes[cpt]['attr_dic'])
            S.add_edge(cpt,tcid,attr_dic=edges[(cpt,tcid)])
        for cand in list(T_com.predecessors(cpt)):
            if cand not in S:
                S.add_node(cand, attr_dic=nodes[cand]['attr_dic'])
            S.add_edge(cand,cpt,attr_dic=edges[(cand,cpt)])
    if len(list(T_com.predecessors(tcid))) != len(nodes[tcid]['attr_dic']['components']):
            for mcpt in nodes[tcid]['attr_dic']['components']:
                if mcpt not in S:
                    S.add_node(mcpt,attr_dic={'tp':'nohit'})
                    S.add_edge(tcid,mcpt,attr_dic={'tp':'nohit'})


# In[6]:


def show_subnetwork(T,raw_tree, node_list,gene_feature_file,condition='linear',name='isolated_system.html',raw_network=False): # T is a refined tree that certain systems has been selected
    H = nx.DiGraph()                                                 # the original tree generated from the table must be provided
    nodes = dict(raw_tree.nodes(data=True))
    edges = dict(((u,v),e) for u,v,e in raw_tree.edges(data=True))
    candidates = [x for x,y in T.nodes(data=True) if y['attr_dic']['tp'] == 'subject']
    systems = [x for x,y in T.nodes(data=True) if y['attr_dic']['tp'] == 'system']
    for nacc in node_list:
        if nacc in systems and nacc not in H:
            H.add_node(nacc,attr_dic=nodes[nacc]['attr_dic'])
            for cpt in list(T.predecessors(nacc)):
                H.add_node(cpt,attr_dic=nodes[cpt]['attr_dic'])
                H.add_edge(cpt,nacc,attr_dic=edges[cpt,nacc])
                for cand in list(T.predecessors(cpt)):
                    if cand not in H:
                        H.add_node(cand,attr_dic=nodes[cand]['attr_dic'])
                    H.add_edge(cand,cpt,attr_dic=edges[cand,cpt])
            if len(list(T.successors(nacc))) > 0: # consider incomplete trees
                for cpt in list(T.successors(nacc)):
                    H.add_node(cpt,attr_dic=nodes[cpt]['attr_dic'])
                    H.add_edge(nacc,cpt,attr_dic=edges[nacc,cpt])
            if len(list(T.predecessors(nacc))) != len(nodes[nacc]['attr_dic']['components']):
                for mcpt in nodes[nacc]['attr_dic']['components']:
                    if mcpt not in H:
                        H.add_node(mcpt,attr_dic={'tp':'nohit'})
                        H.add_edge(nacc,mcpt,attr_dic={'tp':'nohit'})
        if nacc in candidates and nacc not in H:
            first_cpt = list(T.successors(nacc))[0]
            first_tcid = list(T.successors(first_cpt))[0]
            linked_sys = [first_tcid]
            for ls in linked_sys:
                if ls in T:
                    for cpt in list(T.predecessors(ls)):
                        for cand in list(T.predecessors(cpt)):
                            for cp in list(T.successors(cand)):
                                p_sys = list(T.successors(cp))[0]
                                if p_sys not in linked_sys:
                                    linked_sys.append(p_sys)
            for ls in linked_sys:
                add_system_v2(nodes,edges,T,ls,H)
       
    network_visualization_v2(H,gene_feature_file,condition,name,raw_network)


# In[7]:


def get_tcdic(tvt_out,percent):
    tvt = open(tvt_out,'r')
    tdic = {}
    for line in tvt.readlines():
        if '#' in line:
            continue

        li = line.split()
        q = li[0].split('-')[1]
        s = li[1].split('-')[1]
        if q == s:
            continue
        qstart = int(li[2])
        qend = int(li[3])
        qlen = int(li[4])
        sstart = int(li[5])
        send = int(li[6])
        slen = int(li[7])
        qmaxunalign = max(qstart-1,qlen-qend)/qlen
        smaxunalign = max(sstart-1,slen-send)/slen
        if qmaxunalign <= percent and smaxunalign <= percent:
            if q not in tdic:
                tdic[q] = [s]
            else:
                tdic[q].append(s)
    return tdic


# In[8]:


def CrossValidation(q, s, tcdic, G):#query is tc component, subjuect is candidate protein
    for component in list(G.successors(s)):
        if component == q:
            continue
        if G[s][component]['fusion'] == False:
            if component in tcdic:
                if q in tcdic[component]:
                    return True
    return False


# In[9]:


def initialize_raw_network(df,pdic,tcdb,address,pog=0.24):
    if not os.path.isfile('tcdb_vs_tcdb.out'):
        tc_components = list(set((df['tcid'] + '-'+df['query_accession']).tolist()))
        out = open('tcdb_entries.faa', 'wb')
        for tc in tc_components:
            p = subprocess.Popen('blastdbcmd -db tcdb -entry {} -target_only'.format(tc), stdout=subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
            for line in p.stdout:
                out.write(line)
        os.system("blastp -query tcdb_entries.faa -subject tcdb_entries.faa -use_sw_tback -out tcdb_vs_tcdb.out -evalue 1e-6 -outfmt '7 qacc sacc qstart qend qlen sstart send slen length evalue pident'")
    tcdic = get_tcdic('tcdb_vs_tcdb.out',pog)
    G = nx.DiGraph()
    mc = get_mc(tcdb)
    fusion_candidates = []
    for index,row in df.iterrows():
        if not math.isnan(float(row['evalue'])):
            if row['tcid'] not in G:
                G.add_node(row['tcid'],tp='system',components=mc[row['tcid']]) 
            if row['subject_accession'] not in G:
                G.add_node(row['subject_accession'],tp='subject', length = row['subject_length'],tms = row['subject_tms'])
            if row['query_accession'] not in G:
                G.add_node(row['query_accession'],tp='component', length = row['query_length'],tms = row['query_tms'])
        # assign potential fusions
        # definition: largest unaligned potion/lengh > 23%
            q_cood,s_cood = get_cood(row['query_accession'], row['subject_accession'], address)
            if q_cood[0]-0 > row['query_length']-q_cood[1]:
                q_max_unalign = [0,q_cood[0]]
            else:
                q_max_unalign = [q_cood[1],row['query_length']]
            if s_cood[0]-0 > row['subject_length']-s_cood[1]:
                s_max_unalign = [0,s_cood[0]]
            else:
                s_max_unalign = [s_cood[1],row['subject_length']]
            q_portion = (q_max_unalign[1]-q_max_unalign[0])/row['query_length']
            s_portion = (s_max_unalign[1]-s_max_unalign[0])/row['subject_length']
            ds = alignment_has_domain(row['query_accession'], row['subject_accession'],pdic,address)
            if (q_portion <= 0.23 and s_portion <= 0.23):
                G.add_edge(row['subject_accession'],row['query_accession'],tp = 'q_vs_s',evalue=row['evalue'],qcov = row['query_coverage'],scov =row['subject_coverage'],fusion=False,weight=generate_weight(row['query_accession'],row['subject_accession'],row['evalue'],
                row['query_coverage'],row['subject_coverage'],pdic,address))
            elif ds[0] == True and ds[1] == True:
                qde,sde = domain_extension(row['query_accession'], row['subject_accession'],pdic,address)
                if (qde[1]-qde[0])/row['query_length'] >= 0.77 and (qde[1]-qde[0])/(q_cood[1]-q_cood[0]) <=2 and (sde[1]-sde[0])/row['subject_length'] >= 0.77 and (sde[1]-sde[0])/(s_cood[1]-s_cood[0]) <= 2:
                    G.add_edge(row['subject_accession'],row['query_accession'],tp = 'q_vs_s',evalue=row['evalue'],qcov = row['query_coverage'],scov =row['subject_coverage'],fusion=False,weight=generate_weight(row['query_accession'],row['subject_accession'],row['evalue'],
                    row['query_coverage'],row['subject_coverage'],pdic,address))
                else:
                    G.add_edge(row['subject_accession'],row['query_accession'],tp = 'q_vs_s',evalue=row['evalue'],qcov = row['query_coverage'],scov =row['subject_coverage'],fusion=True,weight=generate_weight(row['query_accession'],row['subject_accession'],row['evalue'],
                    row['query_coverage'],row['subject_coverage'],pdic,address))
                    fusion_candidates.append([row['subject_accession'],row['query_accession']])
            else:
                G.add_edge(row['subject_accession'],row['query_accession'],tp = 'q_vs_s',evalue=row['evalue'],qcov = row['query_coverage'],scov =row['subject_coverage'],fusion=True,weight=generate_weight(row['query_accession'],row['subject_accession'],row['evalue'],
                row['query_coverage'],row['subject_coverage'],pdic,address))
                fusion_candidates.append([row['subject_accession'],row['query_accession']])
            G.add_edge(row['query_accession'],row['tcid'],tp = 'hit',title='',weight=1)
        else:
            if row['tcid'] not in G:
                G.add_node(row['tcid'],tp='system',components=mc[row['tcid']])
            G.add_node(row['query_accession'],tp='component', length = row['query_length'],tms = row['query_tms'])
            G.add_edge(row['tcid'],row['query_accession'],tp='nohit',title='Nohit',weight=1)

    
    for fc in fusion_candidates:
        if CrossValidation(fc[1],fc[0],tcdic,G) == True:
            G[fc[0]][fc[1]]['fusion'] = False
    
    return G


# In[10]:


def check_complete(nodes,G, s): # s for a system node in graph, nodes is a dic of all node 
    complete = True
    missing = []
    if len(nodes[s]['components']) != len(list(G.predecessors(s))):
        missing = [x for x in nodes[s]['components'] if x not in list(G.predecessors(s))]
        return [False,missing]
    else:
        for cpt in list(G.predecessors(s)):
            if len(list(G.predecessors(cpt))) == 0:
                complete = False
                missing.append(cpt)
        return [complete,missing]

# recommend system , returns two new trees that attributes are different than the original tree. added attr_dic = {all original attributes}
def raw_recom(H,gene_dic,degree_of_tolerance=2,condition='linear'): 
    nodes = dict(H.nodes(data=True))
    edges = dict(((u,v),e) for u,v,e in H.edges(data=True))
    candidates = [x for x,y in H.nodes(data=True) if y['tp'] == 'subject']
    systems = [x for x,y in H.nodes(data=True) if y['tp'] == 'system']
    cpts = [x for x,y in H.nodes(data=True) if y['tp'] == 'component']
    raw_com = nx.DiGraph() # build up a network that contains only complete systems found from the raw_network
    T = nx.DiGraph() # inlcudes every node, but structures/edges are re-assinged/discarded after the recommand process

    for s in systems:
        situation = check_complete(nodes,H,s)
        if situation[0] == True:
            add_system(nodes,edges,H,s,raw_com)
 
    raw_com_candidates = [x for x,y in raw_com.nodes(data=True) if y['attr_dic']['tp'] == 'subject']
    for c in raw_com_candidates:
        length,path = nx.single_source_dijkstra(raw_com,c)
        #print (length)
        # get the shortest lengh & path of a system + degree_of_tolerance
        node = list(length.keys())
        
        cpts_len = {k:length[k] for k in node if k in cpts}
        #sys_path = {k:path[k] for k in node if k in systems}
        # check if multiple len recorded in the system and choose the shortest lengths
        key_min = min(cpts_len.keys(), key=(lambda k: cpts_len[k]))
        adj_len = cpts_len[key_min] + degree_of_tolerance
        #print(adj_len)
        cpts_len_con = {k:cpts_len[k] for k in cpts_len if cpts_len[k]<=adj_len}
        for nd in cpts_len_con:
            for n in path[nd]:
                if n not in T:
                    T.add_node(n,attr_dic=nodes[n])

            # check if edge has already been added
            if not T.has_edge(path[nd][0],path[nd][1]):
                T.add_edge(path[nd][0],path[nd][1],attr_dic=edges[path[nd][0],path[nd][1]])
            tcnode = list(H.successors(nd))[0]
            if tcnode not in T:
                T.add_node(tcnode,attr_dic=nodes[tcnode])

            if not T.has_edge(nd,tcnode):
                T.add_edge(nd,tcnode,attr_dic=edges[nd,tcnode])

            

    #systems = [x for x,y in H.nodes(data=True) if y['tp'] == 'system']
    raw_com_systems = [x for x,y in raw_com.nodes(data=True) if y['attr_dic']['tp'] == 'system']
    for s in raw_com_systems:
        situation = check_complete(nodes,T,s)
        if situation[0] == True:
            print(s + ' is complete: ')
            for cpt in list(T.predecessors(s)):
                print('component '+cpt+' has candidates: ')
                for cand in list(T.predecessors(cpt)):
                    if edges[(cand,cpt)]['fusion'] == True:
                        print(cand+'(fusion)')
                    else:
                        print(cand)
            print('\n')
        else:
            for cpt in situation[1]:
                if cpt in nodes:
                    T.add_node(cpt,attr_dic=nodes[cpt]) #TCDB version error, data not coherence !!!!!!! Q746T6
                else:
                    T.add_node(cpt,attr_dic={'tp':'component', 'length' : 9999,'tms' : 9999})
                T.add_edge(s,cpt,attr_dic={'tp':'nohit','title':'Nohit','weight':1})
        
    print('The total number of complete systems is: '+str(count_complete_systems(T,nodes)))

    return T


# In[11]:


def network_visualization_v2(G,gene_feature_file,condition = 'linear',name='ms_test.html',raw_network=False): # works on every tree except the original tree generated from the table data
    ms_net = Network(height = '100%', width = '100%', directed = True, notebook = False, bgcolor = '#222222', font_color = 'white')
    ms_net.barnes_hut(gravity=-100, central_gravity=0.005, spring_length=250, spring_strength=0.05, damping=0.7, overlap=0)
    # add edges
    #ms_net.prep_notebook()
    # nodes grouped by its genetic context
    nodes = dict(G.nodes(data=True))
    edges = dict(((u,v),e) for u,v,e in G.edges(data=True))
    groups = {}
    systems = [x for x,y in G.nodes(data=True) if y['attr_dic']['tp'] == 'system']
    cands = [x for x,y in G.nodes(data=True) if y['attr_dic']['tp'] == 'subject']
    gene_dic = get_gene_feature(gene_feature_file)
    groups = {c:'#1ad8f6' for c in cands}
    for s in systems:
        all_cands = []
        for cpt in list(G.predecessors(s)):
            for cand in list(G.predecessors(cpt)):
                if cand not in all_cands:
                    all_cands.append(cand)
        if len(all_cands) > 1:
            gene_distance = get_genetic_distance(gene_dic, all_cands,condition)
            for gene in gene_distance:
                for g in gene_distance[gene]:
                    if g == gene:
                        continue
                    if gene_distance[gene][g] <= 15:
                        groups[gene] = '#f6f31a'
                        groups[g] = '#f6f31a'
                                
    
    #print(edges)
    if raw_network == False:
        blue_edges = [(u,v) for u,v,e in G.edges(data=True) if e['attr_dic']['tp'] == 'q_vs_s']
        white_edges = [(u,v) for u,v,e in G.edges(data=True) if e['attr_dic']['tp'] == 'hit']
        red_edges = [(u,v) for u,v,e in G.edges(data=True) if e['attr_dic']['tp'] == 'nohit']
    else:
        blue_edges = [(u,v) for u,v,e in G.edges(data=True) if e['tp'] == 'q_vs_s']
        white_edges = [(u,v) for u,v,e in G.edges(data=True) if e['tp'] == 'hit']
        red_edges = [(u,v) for u,v,e in G.edges(data=True) if e['tp'] == 'nohit']
        
    for b in blue_edges:
        ms_net.add_node(b[0],b[0],title='# of TMS: {}<br>'.format(nodes[b[0]]['attr_dic']['tms']),shape='dot',physics=True,color = groups[b[0]],mass = 20)
        ms_net.add_node(b[1],b[1],shape='triangle',color='green',title='# of TMS: {}<br>'.format(nodes[b[1]]['attr_dic']['tms']),physics=True,mass = 10)
        if raw_network == False:
            evalue = edges[(b[0],b[1])]['attr_dic']['evalue']
        else:
            evalue = edges[(b[0],b[1])]['evalue']
        if evalue == 0:
            n_e = 1000
        else:
            n_e = -1*math.log(evalue)
        if n_e >= 50:
            wt = 4
        elif n_e >= 25:
            wt = 3
        elif n_e >= 15:
            wt = 2
        else:
            wt = 1
        if raw_network == False:   
            if edges[(b[0],b[1])]['attr_dic']['fusion'] == True:
                ms_net.add_edge(b[0],b[1],width = wt,physics=True,color='purple' )
            else:
                ms_net.add_edge(b[0],b[1],width = wt,physics=True,color = '#05c1f8')
        else:
            if edges[(b[0],b[1])]['fusion'] == True:
                ms_net.add_edge(b[0],b[1],width = wt,physics=True,color='purple' )
            else:
                ms_net.add_edge(b[0],b[1],width = wt,physics=True, color = '#05c1f8')
    
    for w in white_edges:
        if raw_network == False:
            ms_net.add_node(w[1],w[1],shape='square',color='red',title='# of components: {}, missing '.format(len(nodes[w[1]]['attr_dic']['components'])),physics=False,mass = 50)
        else:
            ms_net.add_node(w[1],w[1],shape='square',color='red',title='# of components: {}, missing '.format(len(nodes[w[1]]['attr_dic']['components'])),physics=False,mass = 50)
        ms_net.add_edge(w[0],w[1],width = 1,color = 'white',physics=True)

    if len(red_edges) > 0: # to visualize in-complete systems
        for r in red_edges:
            if r[0] not in ms_net.nodes:
                ms_net.add_node(r[0],r[0],shape='square', color='red',title='# of components: {}, missing '.format(len(nodes[w[1]]['attr_dic']['components'])),physics=False,mass = 50)
            ms_net.add_node(r[1],r[1],shape='triangle',hidden = False, color='green',title='Nohit. # of TMS: {}<br>'.format(nodes[b[1]]['attr_dic']['tms']),physics=True,mass = 10)
            ms_net.add_edge(r[0],r[1],width = 1,physics=True)       
        
    neighbor_map = ms_net.get_adj_list()
    for node in ms_net.nodes:
        node['title'] += "<br>".join(neighbor_map[node['id']])
        #print(node)
    ms_net.set_edge_smooth('discrete')
    ms_net.show_buttons(['physics'])
    ms_net.save_graph(name)
    # open the saved network and add navigation buttons
    networkfile = open(name,'r')
    contents = networkfile.readlines()
    networkfile.close()
    index = 0
    style = 0
    button = 0
    nb = 0
    button_function = 0
    config = 0
    for line in contents:
        if '        "dragNodes": true,' in line:
            nb = index
        elif '<style type="text/css">' in line:
            style = index
        elif '<div id = "mynetwork"></div>' in line:
            button = index
        elif 'var options, data;' in line:
            button_function = index
        elif '#config {' in line:
            config = index
        index += 1
    contents.insert(config+1,'background-color: #ffffff;position: fixed;margin-top: 60px;display: none;')
    contents.insert(nb+1,'        "navigationButtons": true,')
    contents.insert(style+1,'.button {position:absolute;top: 10px;left: 10px;background-color: #4CAF50;border: none;color: white;padding: 15px 32px;text-align: center;text-decoration: none;display: inline-block;font-size: 16px;margin: 4px 2px;cursor: pointer;}')
    contents.insert(button+3,'<button class="button" onclick="showphysics()">Animation Settings</button>')
    contents.insert(button_function+4,'    function showphysics() {var x = document.getElementById("config");if (x.style.display === "none") {x.style.display = "block";} else {x.style.display = "none";}}')
    contents.remove('            width: 400px;\n')
    contents.remove('            height: 600px;\n')
    networkfile = open(name, "w")
    contents = "".join(contents)
    networkfile.write(contents)
    networkfile.close()
    


# In[12]:


def global_selection(G,T_com,gene_dic,pfam_dic,df,address,green_name,yellow_name,condition='linear',max_cycle=10): # df is raw table data
    systems = [x for x,y in T_com.nodes(data=True) if y['attr_dic']['tp'] == 'system']
    nodes = dict(G.nodes(data=True))
    block = []
    total = 0
    green_cases = [] #list of tuples (subject, query)
    yellow_cases = [] #list of tuples
    discarded = [] # list of unassigned hits
    for s in systems:
        if s not in block:
            count = 0
            temp = selection(G,T_com,s,gene_dic,pfam_dic,block,count,green_cases,yellow_cases,discarded,address,condition,max_cycle)
            total = total + temp
    print('The total number of potentially complete systems is: '+str(count_complete_systems(T_com,nodes))) 
    print('The total number of complete systems is (rescued included): '+str(total))
    green_index = []
    yellow_index = []
    for index,row in df.iterrows():
        if (row['subject_accession'],row['query_accession']) in green_cases:
            green_index.append(index)
        if (row['subject_accession'],row['query_accession']) in yellow_cases:
            yellow_index.append(index)
    df_green = df.loc[green_index]
    df_green.to_csv(green_name, sep='\t')
    if len(yellow_index) > 0:
        df_yellow = df.loc[yellow_index]
        df_yellow.to_csv(yellow_name, sep='\t')
    return green_cases, yellow_cases, discarded
    
    
def verify_fusion(S,address):
    cands = [x for x,y in S.nodes(data=True) if y['attr_dic']['tp'] == 'subject']
    for cand in cands:
        outgoing_edges = S.out_edges(cand,data=True)
        fusion_edges = [ie for ie in outgoing_edges if ie[2]['attr_dic']['fusion']==True]
        if len(fusion_edges) == 0:
            continue
        candil = [fe for fe in fusion_edges if fe[2]['attr_dic']['qcov'] > fe[2]['attr_dic']['scov']] # candidate is larger
        if len(candil) == 0:
            continue
        else: # determine if the candidate hit multiple different components in the same system
            # remove edges that this candiate hits only once
            hit_dic = {}
            for cdil in candil:
                hit_sys = list(S.successors(cdil[1]))[0]
                if hit_sys not in hit_dic:
                    hit_dic[hit_sys] = {}
                    hit_dic[hit_sys]['cpts'] = [cdil[1]]
                    hit_dic[hit_sys]['counts'] = 1
                else:
                    hit_dic[hit_sys]['cpts'].append(cdil[1])
                    hit_dic[hit_sys]['counts'] = hit_dic[hit_sys]['counts'] + 1
            
            for hit in hit_dic:
                if hit_dic[hit]['counts'] == 1:
                    if S.has_edge(cand,hit_dic[hit]['cpts'][0]):
                        S.remove_edge(cand,hit_dic[hit]['cpts'][0])
                        #rescue.remove_edge(cand,hit_dic[hit]['cpts'][0])
                else: #verify true fusion
                    complemant = False
                    s_cood_mark = []
                    for cpt in hit_dic[hit_sys]['cpts']:
                        if complemant == True:
                            continue
                        q_cood,s_cood = get_cood(cpt, cand,address)
                        if len(s_cood_mark) == 0:
                            s_cood_mark = s_cood
                        else:
                            if s_cood[0] > s_cood_mark[1] or s_cood_mark[0] > s_cood[1]: #no overlap
                                complemant = True
                            elif (s_cood[0] < s_cood_mark[0] and s_cood[1] > s_cood_mark[1]) or (s_cood[0] > s_cood_mark[0] and s_cood[1] < s_cood_mark[1]):
                                complemant = False
                            else: # overlap
                                overlap_cood = s_cood + s_cood_mark
                            # check the overlap of coordinates
                                overlap_cood.sort()
                                if (overlap_cood[2]-overlap_cood[1])/(overlap_cood[3]-overlap_cood[0]) <= 0.5:
                                    complemant = True
                    if complemant == False: # remove all edges
                        for cpt in hit_dic[hit_sys]['cpts']:
                            if S.has_edge(cand,cpt):
                                S.remove_edge(cand,cpt)
                                #rescue.remove_edge(cand,cpt)
    # deal with situation that component is larger
    cpts = [x for x,y in S.nodes(data=True) if y['attr_dic']['tp'] == 'component' and len(list(S.successors(x)))==1]
    for comp in cpts:
        in_edges = S.in_edges(comp,data=True)
        fusion_edges = [ie for ie in in_edges if ie[2]['attr_dic']['fusion']==True]
        if len(fusion_edges) == 0:
            continue
        compil = [fe for fe in fusion_edges if fe[2]['attr_dic']['qcov'] < fe[2]['attr_dic']['scov']] # component is larger
        if len(compil) == 0:
            continue
        elif len(compil) == 1: # remove
            if S.has_edge(compil[0][0],compil[0][1]):
                S.remove_edge(compil[0][0],compil[0][1])
                #rescue.remove_edge(compil[0][0],compil[0][1])
        else: # determine overlap
            complemant = False
            q_cood_mark = []
            for cpil in compil:
                if complemant == True:
                    continue
                q_cood,s_cood = get_cood(comp, cpil[0],address)
                if len(q_cood_mark) == 0:
                    q_cood_mark = q_cood
                else:
                    if q_cood[0] > q_cood_mark[1] or q_cood_mark[0] > q_cood[1]: #no overlap
                        complemant = True
                    elif (q_cood[0] < q_cood_mark[0] and q_cood[1] > q_cood_mark[1]) or (q_cood[0] > q_cood_mark[0] and q_cood[1] < q_cood_mark[1]):
                            complemant = False
                    else:
                        overlap_cood = q_cood + q_cood_mark
                        overlap_cood.sort()
                        if (overlap_cood[2]-overlap_cood[1])/(overlap_cood[3]-overlap_cood[0]) <= 0.5:
                            complemant = True
            if complemant == False: # remove all edges
                for cpil in compil:
                    if S.has_edge(cpil[0],comp):
                        S.remove_edge(cpil[0],comp)
                        #rescue.remove_edge(cpil[0],comp)
                        
                        
# define the selection function. It will assgin the detailed candidates to a potential complete system
def selection(G,T_com,tcid,gene_dic,pfam_dic,visited_sys,count,green_cases,yellow_cases,discarded,address,condition,max_cycle,cycle=1): # T is a recommended network that only contain all potential complete systems. G is the raw_network, to ensure data extraction easily
    # first, check all candidates in the input system. if a candidate linked to other different systems, bring them in. until all connected systems are analyzed together
    nodes = dict(G.nodes(data=True))
    edges = dict(((u,v),e) for u,v,e in G.edges(data=True))
    #T_edges = dict(((u,v),e) for u,v,e in T_com.edges(data=True)) # Important, dont change
    #print(edges)
    # construct a sub network to include all linked systems and analyze them together
    #S = nx.DiGraph()
    rescue = nx.DiGraph()
    if cycle == 1:
        S = nx.DiGraph()
        linked_sys = [tcid]
        for ls in linked_sys:
            if ls in T_com:
                for cpt in list(T_com.predecessors(ls)):
                    for cand in list(T_com.predecessors(cpt)):
                        for cp in list(T_com.successors(cand)):
                            p_sys = list(T_com.successors(cp))[0]
                            if p_sys not in linked_sys:
                                linked_sys.append(p_sys)
        for ls in linked_sys:
            visited_sys.append(ls)
            add_system(nodes,edges,T_com,ls,S)
            add_system(nodes,edges,T_com,ls,rescue)
    else:
        S = T_com
        T_com_sys = [x for x,y in T_com.nodes(data=True) if y['attr_dic']['tp'] == 'system']
        for tcs in T_com_sys:
            #add_system(nodes,edges,T_com,tcs,S)
            add_system(nodes,edges,T_com,tcs,rescue)
    print('The total number of potentially complete systems in cycle {} is: '.format(cycle)+str(count_complete_systems(S,nodes)))
    #print(S.nodes(data=True))
    cands = [x for x,y in S.nodes(data=True) if y['attr_dic']['tp'] == 'subject']
    verify_fusion(S,address) # removed unqualifed fusions

                        
    # for the complete subnetwork S, perform pedigree function. pedigree will generate a dic showing the intersection of membrane protein candidates 
    # inside the network. if several system share a membrane candidate, or several membrane candidates, they will be analyzed together by paternity_test
    if cycle == 1:
        pedigree(S, linked_sys,gene_dic,pfam_dic,address,condition)
    else:
        pedigree(S, T_com_sys,gene_dic,pfam_dic,address,condition)
    
    # check / print complete systems
    # give graph
    verify_fusion(S,address) # removed unassigned fusions
    
    systems = [x for x,y in S.nodes(data=True) if y['attr_dic']['tp'] == 'system']
    for s in systems:
        situation = check_complete(nodes,S,s)
        if situation[0] == True:
            print('(cycle {})'.format(cycle) + s + ' is complete: ')
            rescue.remove_node(s)
            for cpt in list(S.predecessors(s)):
                print('component '+cpt+' has candidates: ')
                rescue.remove_node(cpt)
                for cand in list(S.predecessors(cpt)):
                    if S[cand][cpt]['attr_dic']['fusion'] == True:
                        print(cand + ' (fusion)')
                    else:
                        print(cand)
                    if cand in rescue:
                        rescue.remove_node(cand)
                    # add this to the table outputs
                    if cycle == 1:
                        green_cases.append((cand,cpt))
                    else:
                        yellow_cases.append((cand,cpt))
            print('\n')
        else:
            print('(cycle {})'.format(cycle) + s + ' is incomplete: ')
            for cpt in list(S.predecessors(s)):
                print('component '+cpt+' has candidates: ')
                for cand in list(S.predecessors(cpt)):
                    print(cand)
            for ocpt in list(S.successors(s)):
                print('component '+ocpt+' has no candidates!')
            print('\n')
            
    r_systems = [x for x,y in rescue.nodes(data=True) if y['attr_dic']['tp'] == 'system']
    for rs in r_systems:
        situation = check_complete(nodes,rescue,rs)
        if situation[0] == False:
            rescue.remove_node(rs)
            for rcpt in nodes[rs]['components']:
                for rcand in list(rescue.predecessors(rcpt)):
                    rescue.remove_edge(rcand,rcpt)
                rescue.remove_node(rcpt)
    print('The total number of complete systems in cycle {} is: '.format(cycle)+str(count_complete_systems(S,nodes)))
    print('\n')
    count = count + count_complete_systems(S,nodes)
    if cycle == max_cycle:
        print('Maximum cycle reached, no convergence observed, human intervention required for above cases!!!!!')
        print('\n')
        return count
    print('Candidates that are discarded/Human intervention required in cycle {}:'.format(cycle))
    for cand in cands:
        if len(list(S.successors(cand))) == 0: # make a bet function here. try to see if it can be remaining in the resuce network
            #print('{} is removed from the initial assignment'.format(cand,cycle))
            if cand in rescue:
                # logic of bet: if a hit is discarded, it is either too good to be extra, or just being not qualified fusion(some part has not hit)
                # if a canidate, 1-portion of largest unaligned part is > 0.5, change all its edges to blue
                if len(list(rescue.successors(cand))) == 0:
                    rescue.remove_node(cand)
                    if cycle == 1:
                        print('{} is discarded in cycle {} as a good extra hit.'.format(cand,cycle) )
                        discarded.append(cand)
                    else:
                        print('{} is discarded in cycle {} as a human intervention case.'.format(cand,cycle) )
                    continue
                for succ in list(rescue.successors(cand)):
                    if abs(int(nodes[succ]['tms'])-int(nodes[cand]['tms'])) < 3: # start bet function
                        ds = alignment_has_domain(succ,cand,pfam_dic,address)
                        if ds[0] == True and ds[1] == True:
                            rescue[cand][succ]['attr_dic'] = {'tp' : 'q_vs_s','evalue':edges[(cand,succ)]['evalue'],'qcov':edges[(cand,succ)]['qcov'],'scov' :edges[(cand,succ)]['scov'],'fusion':False,'weight':edges[(cand,succ)]['weight']}
                            
                        else:
                            if rescue.has_edge(cand,succ):
                                rescue.remove_edge(cand,succ)
                    else:
                        if rescue.has_edge(cand,succ):
                            rescue.remove_edge(cand,succ)
                if len(list(rescue.successors(cand))) == 0:
                    rescue.remove_node(cand)
                    print('{} is discarded in cycle {} as a human intervention case.'.format(cand,cycle) )
                    discarded.append(cand)
                else:
                    print('{} is rescued and put it into the next cycle.'.format(cand))
                
                
    
    print('\n')
    if count_complete_systems(rescue,nodes) == 0:
        return count
    else:
        return selection(G,rescue,tcid,gene_dic,pfam_dic,visited_sys,count,green_cases,yellow_cases,discarded,address,condition,max_cycle=max_cycle,cycle=cycle+1)
    #return S,count_complete_systems(S),count_complete_systems(R),systems
    
def pedigree(S, linked_sys,gene_dic,pfam_dic,address,condition):
    # generate a dic showing intersection of shared candidates
    nodes = dict(S.nodes(data=True))
    edges = dict(((u,v),e) for u,v,e in S.edges(data=True))
    paternity = {}
    candidates = [x for x,y in S.nodes(data=True) if y['attr_dic']['tp'] == 'subject']
    for cand in candidates:
        cpts_involved = list(S.successors(cand))
        cpts_involved.sort() # to eliminate the effect of insertion order
        if len(cpts_involved) == 1:
            continue # this is the case that candidate is not shared
        # if shared
        sys_involved = []
        for cpt_involved in cpts_involved:
            sys_involved.append(list(S.successors(cpt_involved))[0]) # doesnt matter if system in this list duplicates. It means a candidate hits diffrent components in the same system
        sys_involved = list(set(sys_involved)) # remove duplicates
        # sort and make as a tuple to be the key of paternity dic
        sys_involved.sort()
        key = tuple(sys_involved)
        if key not in paternity:
            paternity[key] = {}
            paternity[key][tuple(cpts_involved)] = [cand]
        else:
            # check if cpts involved are same or not
            if not tuple(cpts_involved) in paternity[key]: # a new set of cpts that from different sys sharing a cand
                paternity[key][tuple(cpts_involved)] = [cand]
            else: # a same set of cpts from same set of sys share another new cand
                paternity[key][tuple(cpts_involved)].append(cand)

                    
    #print(paternity)
    # prepare for paternity_test
    # when performing paternity_test and looping the dictionary: dic_key will be sorted by length of keys, which means linked TCID keys will be examed first
    test_groups = [x for x in list(paternity.keys()) if len(x) > 1]
    test_groups_complexity = {}
    for tg in test_groups:
        if tg not in test_groups_complexity:
            test_groups_complexity[tg] = 0
        temp = 0
        for pairs in paternity[tg]:
            local = len(paternity[tg][pairs])
            if local > temp:
                temp = local
        test_groups_complexity[tg] = temp
    sorted_tgc = sorted(test_groups_complexity.items(), key=operator.itemgetter(1),reverse=True) # this is the order to keep edges
    tgc = dict(sorted_tgc)
    test_sequence = list(tgc.keys())
    #test_sequence = sorted(test_groups,key=len,reverse=True)
    #print(test_sequence)
    for t in test_sequence: # convert back from tuple to list later
        for rc in paternity[t]: # rc for related cpts
            paternity_test(S,nodes,edges,list(rc), paternity[t][rc],gene_dic,pfam_dic,address,condition) # a dictionary contain all assignments
            #print (complete_sys)
    cpts = [x for x,y in S.nodes(data=True) if y['attr_dic']['tp'] == 'component']
    for cpt in cpts:
        if len(list(S.predecessors(cpt))) == 0: #incomplete sys
            tcid = list(S.successors(cpt))[0]
            S.remove_edge(cpt,tcid)
            S.add_edge(tcid,cpt,attr_dic={'tp':'nohit','title':'Nohit','weight':1})
            
    
    
    

def paternity_test(S,nodes,edges,children, parents,gene_dic,pfam_dic,address,condition): #children: cpts_involved, parents: candidates_involved
    # paternity test will remove inrelavet edges in S if a candidate is assigned
    # paternity test will not determine if a system is complete or not. It only solves conflicts among shared candidates
    # paternity test should also solve fusion in both directions
    # if a candidate has severl outgoing edges, some of them are purple(fusion), blue edges have the higher priority to be assinged
    #print(children) 
    combination = {}
    fusion = {} # solve shared fusions

    for p in parents:
        for c in children:
            if S.has_edge(p,c) and edges[(p,c)]['attr_dic']['fusion'] == False:
                combination[(p,c)] = test_kit( S,nodes, edges, p, c,gene_dic,pfam_dic,address,condition )
                    
            else:
                combination[(p,c)] = 0
                if p not in fusion:
                    fusion[p] = [c]
                else:
                    fusion[p].append(c)

    sorted_combo = sorted(combination.items(), key=operator.itemgetter(1),reverse=True) # this is the order to keep edges
    combo = dict(sorted_combo)
    ps = []
    ss = [] # different systems involved for only blue edges assignment
    for key in combo:
        if combo[key] > 0:
            if key[0] not in ps:
                ps.append(key[0])
            if list(S.successors(key[1]))[0] not in ss:
                ss.append(list(S.successors(key[1]))[0])
    max_assignment = min(len(ps),len(ss)) # max number of different system can be assgined for blue edges
    # choose the first n assignments for only blue edges
    assigned_cand = []
    assigned_cpts = []
    assigned_sys = []
    #print(combo)
    for pc_pair in combo: # dimer situation considered
        if combo[pc_pair] > 0: #and pc_pair[0] not in fusion: # if a candidate involves in both fusion and complete cases, it favors non-fusion first
            if pc_pair[0] not in assigned_cand:
                if pc_pair[1] not in assigned_cpts:
                    if len(assigned_sys) < max_assignment:
                        assigned_cand.append(pc_pair[0])
                        assigned_cpts.append(pc_pair[1])
                        assigned_sys=assigned_sys+list(S.successors(pc_pair[1]))
                    else:
                        if S.has_edge(pc_pair[0],pc_pair[1]):
                            S.remove_edge(pc_pair[0],pc_pair[1])
                    
                else:
                    if S.has_edge(pc_pair[0],pc_pair[1]):
                        S.remove_edge(pc_pair[0],pc_pair[1])# remove extra candidate and free them for more potential cases
            
            else:
                # check dimer cases
                if pc_pair[1] not in assigned_cpts and list(S.successors(pc_pair[1]))[0] in assigned_sys:
                    assigned_cpts.append(pc_pair[1])
                elif S.has_edge(pc_pair[0],pc_pair[1]):
                    S.remove_edge(pc_pair[0],pc_pair[1])
        else:
            if pc_pair[0] in assigned_cand: # which means the candidate has both valid fusion edges and blue edges, and the blue edges are assigned in previous steps
                if S.has_edge(pc_pair[0],pc_pair[1]):
                    S.remove_edge(pc_pair[0],pc_pair[1])
    
    for f in fusion: # solve for shared fusion
        if S.has_edge(f,fusion[f][0]): # which means this candidate is not assgined for any blue edges in previous step
            if len(fusion[f]) == 1:
                continue # directly assign
            else:
                # determine how many different systems involved
                conflict = {}
                for fcpt in fusion[f]:
                    f_sys = list(S.successors(fcpt))[0]
                    if fcpt in assigned_cpts:
                        S.remove_edge(f, fcpt)
                    else:
                        if f_sys not in conflict:
                            conflict[f_sys] = [fcpt]
                        else:
                            conflict[f_sys].append(fcpt)
                if len(list(conflict.keys())) == 1: #multiple components fuse to the candidate, directly assgin
                    continue
                else: # choose the best pair by normalized progressive score, remove the rest edges
                    best_score = 0
                    best_sys = ''
                    
                    for cf in conflict:
                        local = 0
                        for cfcpt in conflict[cf]:
                            temp_score = test_kit(S,nodes,edges,f,cfcpt,gene_dic,pfam_dic,address,condition)
                            if temp_score > local:
                                local = temp_score
                        if local > best_score:
                            best_score = local
                            if len(best_sys) == 0:
                                best_sys = cf
                            else:
                                for cf_cpt in conflict[best_sys]:
                                    S.remove_edge(f,cf_cpt)
                                best_sys = cf
                                
                        else:
                            for cf_cpt in conflict[cf]:
                                S.remove_edge(f,cf_cpt)

                    win = True
                    new_best = best_score
                    best_blue = ''
                    #print(conflict)
                    for c in children:
                        if combo[(f,c)] > 0:
                            if combo[(f,c)] > new_best: # in these cases purple edges are not as good as blue edges
                                win = False
                                new_best = combo[(f,c)]
                                if len(best_blue) == 0:
                                    best_blue = c
                                else:
                                    S.remove_edge(f,best_blue)
                                    best_blue = c
                            else: # remove the correspoding blue edges
                                S.remove_edge(f,c)
                    if win == False: # remove all purple edges. If win == True, in previous loop all blue edges are gone                       
                        for sys in conflict:
                            for cpt in conflict[sys]:
                                if S.has_edge(f,cpt):
                                    S.remove_edge(f,cpt)
                        
        
    #print(max_sys)
    
    
def test_kit(S,nodes,edges,p,c,gene_dic,pfam_dic,address,condition): 
    # calculate progressive score to this pair of p-c
    score = progressive_score( S, nodes, edges, p, c,gene_dic,pfam_dic,address,condition)
    return score


# In[13]:


def calculate_individual_basic_score(S, nodes, edges, p, c,pfam_dic,address): # this is only based on single_edge_confidence
    evalue = edges[(p,c)]['attr_dic']['evalue']
    if evalue == 0:
        normalized_e = 100
    else:
        normalized_e = -1 * math.log(evalue,10)
    if normalized_e > 100:
        normalized_e = 100  
    else:
        normalized_e = 100 * (normalized_e/100)
    score = normalized_e + max(edges[(p,c)]['attr_dic']['qcov'],edges[(p,c)]['attr_dic']['scov'])
    situations = alignment_has_domain(c,p,pfam_dic,address)
    if situations[0] == True and situations[1] == True:
        score = score + 100
    return score
    


# In[14]:


def calculate_overall_basic_score(S, nodes, edges, p, c,pfam_dic,address): # this is based on the overall confidence.
    # for each component of that system, a best representative score will be chose.
    # They will add together and normalized by the number of different components that system has
    sys = list(S.successors(c))[0]
    other_cpts = list(S.predecessors(sys))
    total_score = calculate_individual_basic_score(S, nodes, edges, p, c,pfam_dic,address)
    for cpt in other_cpts:
        if cpt == c:
            continue
        best_local_score = 0
        for cand in list(S.predecessors(cpt)):
            temp_score = calculate_individual_basic_score(S, nodes, edges, cand, cpt,pfam_dic,address)
            if temp_score > best_local_score:
                best_local_score = temp_score # get the representative score for each component
        total_score = total_score + best_local_score
    normalized_score = total_score/len(nodes[sys]['attr_dic']['components'])
    return normalized_score


# In[15]:


def calculate_bonus_score(S, nodes, edges, p, c, p_vs_all,max_bonus=100):
    # max is 100
    # distance matrix should not include other candiates under the same component which is being investigated
    bonus = 0
     # check consecutivity of the rest of gene
    if p_vs_all[list(p_vs_all.keys())[1]] > 15:
        distances = list(p_vs_all.values())
        consec_dis = [list(group) for group in mit.consecutive_groups(distances)] # this consecutivity will never include the current candidate
        if len(consec_dis) < len(distances):
            bonus = bonus + 0.1*max_bonus
            return bonus
        else:
            return bonus
    else:
        bonus = bonus + 0.4*max_bonus
        distances = list(p_vs_all.values())
        consec_dis = [list(group) for group in mit.consecutive_groups(distances)] # this consecutivity will never include the current candidate
        if len(consec_dis) == len(distances): # nothing continuous
            return bonus
        else:
            # check if how good the consuctivity
            c_l = 0
            for d in consec_dis:
                if len(d) > 1:
                    c_l = c_l + len(d)
            bonus = bonus + 0.6*max_bonus * c_l/len(distances)   # the percentage that candidates involved in the consecutivity
    
            return bonus


# In[16]:


def progressive_score( S, nodes, edges, p, c ,gene_dic,pfam_dic,address,condition):
    sys = list(S.successors(c))[0]
    other_cpts = nodes[sys]['attr_dic']['components']
    other_cpts = [x for x in other_cpts if not x == c]
    distance_list = [p]
    for cpt in other_cpts:
        for oc in list(S.predecessors(cpt)):
            distance_list.append(oc)
    if len(list(set(distance_list))) == 1:
        # special case: during previous process , the system became incomplete
        return 0
    gene_distance = get_genetic_distance(gene_dic, distance_list,condition)
    x = gene_distance[p]
    sorted_x = sorted(x.items(), key=operator.itemgetter(1))        
    p_vs_all = dict(sorted_x)
    bonus_score = calculate_bonus_score(S, nodes, edges, p, c,p_vs_all)
    if p_vs_all[list(p_vs_all.keys())[1]] > 15:
        basic_score = calculate_individual_basic_score(S, nodes, edges, p, c,pfam_dic,address)
    else:
        basic_score = calculate_overall_basic_score(S, nodes, edges, p, c,pfam_dic,address)
    return basic_score + bonus_score
    


# tblock = []
# tcount=0
# tgreen_cases = [] #list of tuples (subject, query)
# tyellow_cases = [] #list of tuples
# tdiscarded = [] # list of unassigned hits
# temp = selection(G,T,'3.D.7.1.4',gene_dic,pdic,tblock,tcount,tgreen_cases,tyellow_cases,tdiscarded,'/ResearchData/Users/amedrano/RalfRabus/MultiomponentSystems/Desulfococcus_multivorans','linear',max_cycle=10)

# In[17]:


def get_df_red(df,green_cases,yellow_cases,discarded):
    hit_cand = []
    hit_cpts = []
    for g in green_cases:
        if g[0] not in hit_cand:
            hit_cand.append(g[0])
        if g[1] not in hit_cpts:
            hit_cpts.append(g[1])
    for y in yellow_cases:
        if y[0] not in hit_cand:
            hit_cand.append(y[0])
        if y[1] not in hit_cpts:
            hit_cpts.append(y[1])
    drop_list = []
    for index,row in df.iterrows():
        if row['subject_accession'] in hit_cand or row['query_accession'] in hit_cpts: #or row['subject_accession'] in discarded:
            drop_list.append(index)
    df_red = df.drop(index = drop_list)
    return df_red


# In[18]:


def save_network(G, node_filename, edge_filename):
    saved_nodes = {}
    saved_edges = nx.to_dict_of_dicts(G)
    for s_n in G.nodes(data=True):
        saved_nodes[s_n[0]] = s_n[1]
    json.dump(saved_edges, open(edge_filename,'w'))
    json.dump(saved_nodes, open(node_filename,'w'))


# In[19]:


def load_network(node_filename,edge_filename):
    load_nodes = json.load(open(node_filename))
    load_edges = json.load(open(edge_filename))
    recover_network = nx.from_dict_of_dicts(load_edges, create_using=nx.DiGraph)
    for ln in load_nodes:
        if 'attr_dic' in load_nodes[ln]:
            recover_network.add_node(ln,attr_dic=load_nodes[ln]['attr_dic'])
        else:
            recover_network.add_node(ln,attr_dic=load_nodes[ln])
    return recover_network


# In[20]:


def analyze_inconfident_data(df,green_cases,yellow_cases,discarded,tcdb,address, linear='linear',pog=0.24 ):
    outfile = open('unassigned_hits.txt','w')
    df_red = get_df_red(df,green_cases,yellow_cases,discarded)
    # check potential complete systems in the red cases
    rest = initialize_raw_network(df_red,pdic,tcdb,address,pog) #treat as a new G
    save_network(rest, 'rest_network_nodes.txt', 'rest_network_edges.txt')
    nodes = dict(rest.nodes(data=True))
    edges = dict(((u,v),e) for u,v,e in rest.edges(data=True))
    #print(rest.edges(data=True))
    rest_complete = nx.DiGraph()
    rsystems = [x for x,y in rest.nodes(data=True) if y['tp'] == 'system']
    for s in rsystems:
        situation = check_complete(nodes,rest,s)
        if situation[0] == True:
            add_system(nodes,edges,rest,s,rest_complete)
            print(s + ' is complete: ')
            for cpt in list(rest.predecessors(s)):
                print('component '+cpt+' has candidates: ')
                for cand in list(rest.predecessors(cpt)):
                    print(cand)
            print('\n')
    rgreen_cases,ryellow_cases,rdiscarded=global_selection(rest,rest_complete,gene_dic,pdic,df_red,address,'yellow_reds.tsv','should_not_exist.tsv',linear,10)
    for hit in rdiscarded:
        outfile.write(hit+'\n')
    outfile.close()
    df_nomore = get_df_red(df_red,rgreen_cases,ryellow_cases,rdiscarded)
    df_nomore.to_csv('reds.tsv',sep='\t')


# df = get_rawdata('reportMulticomponentSystems.tsv')

# pdic = get_info(df['query_accession'], df['subject_accession'])

# gene_dic = get_gene_feature('/ResearchData/Users/amedrano/RalfRabus/Genomes/GCA_001854245.1_ASM185424v1/GCA_001854245.1_ASM185424v1_feature_table.txt.gz')

# G = initialize_raw_network(df,pdic,'./tcdb.faa','/ResearchData/Users/amedrano/RalfRabus/MultiomponentSystems/Desulfococcus_multivorans')
# 

# save_network(G, 'raw_network_nodes.txt', 'raw_network_edges.txt')

# T = raw_recom(G,gene_dic)

# save_network(T, 'recom_network_nodes.txt', 'recom_network_edges.txt')

# green_cases,yellow_cases,discarded = global_selection(G,T,gene_dic,pdic,df,'/ResearchData/Users/amedrano/RalfRabus/MultiomponentSystems/Desulfococcus_multivorans','greens.tsv','yellow_greens.tsv')

# analyze_inconfident_data(df,green_cases,yellow_cases,discarded,'tcdb.faa','/ResearchData/Users/amedrano/RalfRabus/MultiomponentSystems/Desulfococcus_multivorans', linear='linear' )

# show_subnetwork(T,G,['1.A.30.2.10','1.A.30.2.7'],name='isolated_system_v3.html')

# In[240]:


if __name__ == '__main__':
    # parse arguments
    parser = argparse.ArgumentParser()
    # create a parser to handle options and arguments. Show help info if no args
    parser.add_argument( '-d', '--data', type = str, dest = 'data', metavar = '<raw data table>', help = 'MANDATORY. the file name or file name with the absolute address of the gblast table for multi-component transport systems in tsv format.' )
    parser.add_argument( '-db', '--database', type = str, dest = 'tcdb', metavar = '<tcdb.faa>', help = 'MANDATORY. the file name or file name with the absolute address of tcdb file in FASTA format.' )
    parser.add_argument( '-ra', '--rootaddress', type = str, dest = 'ra', metavar = '<address>', help = 'MANDATORY. the absolute address of the multi-component system analysis with all relative results generated by Gblast under that directory. eg. "/ResearchData/Users/amedrano/RalfRabus/MultiomponentSystems/Desulfococcus_multivorans" ' )
    parser.add_argument( '-ft', '--featuretable', type = str, dest = 'ft', metavar = '<featuretable.txt.gz>', help = 'MANDATORY. the file name or file name with the absolute address of the genomic feature table in gz compressed format.' )
    parser.add_argument( '-clan', '--clan_info_file', type = str, dest = 'clan', metavar = '<Pfam-A.clans.tsv.gz>', help = 'MANDATORY. the file name or file name with the absolute address of the Pfam table in gz compressed format.' )
    parser.add_argument( '-l', '--linear', action = 'store_false', dest = 'linear', default = True, help = 'Flag. if set, chromosomes of the genome in this analysis are non-linear.' )
    parser.add_argument( '-dot', '--degree_of_tolerance', dest = 'dot', type = int, default = 2, help = 'The degree of tolerance of the network. default is 2, means 2 percent.' )
    parser.add_argument( '-pog','--percent_of_homology',dest='pog', type = float, default=0.24, help = 'The maximum unaligned coverage that 2 TCDB proteins allowed to be recognized as homologous. Default is 0.24, that is, if maximum unaligned parts of 2 TCDB proteins are both smaller than 24 percent, they are homologous.')
    args = parser.parse_args()
    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)

    df = get_rawdata(args.data)
    pdic = get_info(df['query_accession'], df['subject_accession'],args.ra,args.clan)
    gene_dic = get_gene_feature(args.ft)
    G = initialize_raw_network(df,pdic,args.tcdb,args.ra,args.pog)
    # save the initial network
    save_network(G, 'raw_network_nodes.txt', 'raw_network_edges.txt')
    if args.linear == True:
        status = 'linear'
    else:
        status = 'non-linear'
    T = raw_recom(G,gene_dic,args.dot,status)
    # save the network of stage 2
    save_network(T, 'recom_network_nodes.txt', 'recom_network_edges.txt')
    green_cases,yellow_cases,discarded = global_selection(G,T,gene_dic,pdic,df,args.ra,'greens.tsv','yellow_greens.tsv',status)
    analyze_inconfident_data(df,green_cases,yellow_cases,discarded,args.tcdb,args.ra, status,args.pog )


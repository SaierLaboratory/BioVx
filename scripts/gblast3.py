#!/usr/bin/env python
import resource,copy
from Bio.Blast import NCBIXML
from Bio import SeqIO
from ProjectBio import ParseDefline, nr_dict, ParseTC
from Bio.Blast.Applications import NcbiblastpCommandline
from urllib import quote as urlencode
from urllib2 import urlopen
import hmmtop,shutil
from math import ceil
import sys,os,pickle
import subprocess
import matplotlib
import hmmgap
import tcdb
import cdd
from decimal import Decimal
from sys import exit
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab,re,tempfile
rsrc = resource.RLIMIT_DATA
soft, hard = resource.getrlimit(rsrc)
resource.setrlimit(rsrc, (1073741824, hard)) #limit to one gig, omg..

class Tools:

    def __init__(self):
        self.dbfile = '/db/tcdb'
        self.tcdb = os.environ['HOME']+self.dbfile
        self.query = False
        self.indir = False
        self.goodresults  = []
        self.bestresults  = []
        self.notmsresults = []
        self.substrates = {}
        self.debug=True
        self.tms = {}
        self.expect = 0.001
        self.minlen = 50
        self.mincov = 40.0
        self.myqueries = False
        self.mytcdb = False
        self.rows = []
        self.data = []
        self.globalcount = 0
        self.cdd_on = False
        self.abbreviations = {}
        self.ortho = False
        self.query_gis = False
        self.target_gis = False
        self.esort = False
        self.alabel = 'Query'
        self.blabel = 'TCDB-Hit'
        self.names  = tcdb.Names()
        self.loadSubstrates()

        #Sequence dictionaries
        self.queries = {}
        self.tcdbHits = {}        

        #self.substrates = tcdb.Substrates()
        tcdb.use_local()
        tcdb.use_local_betabarrel()

    def blast_all(self):
        try:
            os.makedirs(self.indir+"/xml")
        except:
            pass
        if self.ortho is not False:
            self.prep_orthologs()
        queries = SeqIO.parse(open(self.query),'fasta')
        for query in queries:
            query_file = tempfile.NamedTemporaryFile()
            SeqIO.write(query,query_file,'fasta')
            query_file.flush()
            blast_out = "'"+self.indir+'/xml/'+query.id+".xml"+"'"
            # removes the quotes in a rly cool way :)
            if os.path.exists(blast_out[1:-1]) is True:
                continue # Blast xml already exists...
            blastp = NcbiblastpCommandline(
                query=query_file.name,
                db=self.tcdb,
                evalue=self.expect,
                out=blast_out,
                outfmt=5,
                comp_based_stats='0'
            )
            blastp()
            print "Blasted :: %s" %query.id


    def prep_orthologs(self):
        queries = SeqIO.parse(open(self.query),'fasta')
        targets = SeqIO.parse(open(self.ortho),'fasta')
        # Load query Gis
        query_gis = open(self.query_gis,'r')
        target_gis = open(self.target_gis,'r')
        queries = SeqIO.to_dict(queries)
        targets = SeqIO.to_dict(targets)
        # Make ortho directory
        try:
            os.makedirs(self.indir+"/orthologs")
        except:
            pass
        # Write queries & targets
        myqueries = open(self.indir+'/orthologs/myqueries.faa','wb')
        mytargets = open(self.indir+'/orthologs/mytargets.faa','wb')
        for sgi in query_gis:
            SeqIO.write(queries[sgi.strip()],myqueries,'fasta')
        for tgi in target_gis:
            SeqIO.write(targets[tgi.strip()],mytargets,'fasta')
        myqueries.flush()
        mytargets.flush()
        #os.system('makeblastdb -dbtype prot -in '+mytargets.name)
        subprocess.call(('makeblastdb -dbtype prot -in '+mytargets.name).split())
        self.query = myqueries.name
        self.tcdb = mytargets.name


    def load_good(self):
        db = self.indir+"/goodresults.db"
        if os.path.exists(db) and self.debug:
            self.goodresults = pickle.load(open(db,'r'))
            return
        xml = os.listdir(self.indir+"/xml")
        rez = []
        for res in xml:
            if os.path.exists(self.indir+'/xml/'+res) is False:
                continue
            try:
                for blast in NCBIXML.parse(open(self.indir+'/xml/'+res)):
                    query =  blast.query
                    try:
                        descriptions = [i for i in blast.descriptions]
                        alignments = [i for i in blast.alignments]
                        results = zip(descriptions,alignments)
                        results.sort(key=lambda x:x[0].e,reverse=False)
                        record = results[0] # Contains top description and top alignment object in tuple.
                        hsps = [i for i in record[1].hsps]
                        hsps.sort(key=lambda x:x.expect)
                        hsp = hsps[0]


                        if not blast.query_length:
                                print('No Query Length')        
                                print(vars(blast))

                        if not record[1].length:

                                print('No hit length')
                                print(vars(record[1]))

                        q_len = blast.query_length
                        h_len = record[1].length

                        qcov = float('%.1f'%(((hsp.query_end-hsp.query_start)/float(q_len))*100))
                        hcov = float('%.1f'%(((hsp.sbjct_end-hsp.sbjct_start)/float(h_len))*100))

                        if qcov >= 100:
                
                                qcov = 100.0
                        
                        if hcov >= 100:

                                hcov = 100.0

                        
                        if record[0].e > self.expect or len(hsp.match) < self.minlen:
                                continue
                        

                        if qcov >= self.mincov or hcov >= self.mincov:

                                rez.append((query,record,hsp,q_len,h_len,qcov,hcov)) # (genome ID, hit record <e,title>, hit.hsp)
                    except:
                        pass
            except:
                continue
        self.goodresults = rez
        pickle.dump(self.goodresults,open(db,'wb'))
        return

    def write_fastas(self):
        mytcdb = []
        myquery = []
        tcdb = SeqIO.parse(self.tcdb,'fasta')
        self.tcdbHits = nr_dict(tcdb)
        queries= SeqIO.parse(self.query,'fasta')
        self.queries= SeqIO.to_dict(queries)
        ##### gabo's addition
        done = dict()
        for query,hit,hsp,q_len,h_len,qcov,hcov in self.goodresults:
            hit = hit[0]
            query=ParseDefline(query).id
            myquery.append(self.queries[str(query)])
            subject = ParseDefline(hit.title,True).id
            if subject not in done.keys():
                done[subject] = 1
                try:
                    mytcdb.append(self.tcdbHits[subject])
                except:
                    (family,tcid,acc) = ParseTC(hit.title)
                    #print tcid,acc
                    print ParseDefline(hit.title,True).id
                    #print hit.title
                    quit()

        query_file=open(self.indir+"/myqueries.faa",'wb')
        tcdb_file=open(self.indir+"/mytcdb.faa",'wb')
        SeqIO.write(myquery,query_file,'fasta')
        SeqIO.write(mytcdb,tcdb_file,'fasta')

    def hmmtop(self):
        db = self.indir+'/hmmtop.db'
        if os.path.exists(db) and self.debug:
            self.tms = pickle.load(open(db,'r'))
            return
        ht = hmmtop.tools()
        ht.add_library('queries',self.indir+"/myqueries.faa") # Genome
        ht.add_library('tcdb',self.indir+"/mytcdb.faa")
        ht.scan_libraries()
        pickle.dump(ht.results,open(db,'wb'))
        self.tms = ht.results
        return

    def calculate_tms_scores(self):
        for genome,tcdb,hsp,q_len,h_len,qcov,hcov in self.goodresults:
            genome=ParseDefline(genome).id
            tcdb=tcdb[0]
            delta  = hsp.sbjct_start-hsp.query_start
            tcdbid = ParseDefline(tcdb.title,True).id
            try:
                genome_tms = self.tms['queries'][genome].values()
                tcdb_tms = self.tms['tcdb'][tcdbid].values()

            except KeyError:
                # Genome or TCDB hit dont have a TMS.
                # These must be manually revised later!
                self.notmsresults.append((genome,tcdb,hsp,q_len,h_len,qcov,hcov,None))
                continue
            g_tms = [[i[0]+delta,i[1]+delta] for i in genome_tms]
            overlap = self.find_overlap(g_tms,tcdb_tms)
            row= (genome,tcdb,hsp,q_len,h_len,qcov,hcov,overlap)
            self.bestresults.append(row)

        if self.esort:

                self.bestresults.sort(key=lambda x:(x[1].e,x[7],self.tcsort(x[1].title)))
                self.notmsresults.sort(key=lambda x:(x[1].e,self.tcsort(x[1].title)))

        else:

                self.bestresults.sort(key=lambda x:(self.tcsort(x[1].title),x[1].e,-1.0*max(x[6],x[5])))
                self.notmsresults.sort(key=lambda x:(self.tcsort(x[1].title),x[1].e,-1.0*max(x[6],x[5])))        



    def find_overlap(self,query,target):
        queries = [range(i[0],i[1]+1) for i in query]
        targets = [range(i[0],i[1]+1) for i in target]
        overlap = []
        for sub in queries:
            for tar in targets:
                overlap.extend(set(sub)&set(tar))
        return float(len(set(overlap)))/20

    def title_extract(self,string): #returns: (acc,tcid)
        string = string.split(" ")
        gi = string[0].split('|')[-1]
        tcid = string[1]
        return (gi,tcid)

    def tcsort(self,title):
        if self.ortho is not False:
            return 1
        title=ParseDefline(title,True).description
        (family,tc,acc) = ParseTC(title)
        tc = tc.split('.')
        return ( int(tc[0]), str(tc[1]), int(tc[2]), int(tc[3]), int(tc[4]) )

    def build_view(self,data):
        (genome,tcdb,hsp,q_len,h_len,qcov,hcov,overlap) = data

        try:
            os.mkdir(self.indir+"/img")
        except:
            pass
        genome=ParseDefline(genome).id
        tid = ParseDefline(tcdb.title,True).id
        if os.path.exists(self.indir+"/img/"+genome+".png") is False:

            try:
                san = self.queryhmg(self.myqueries[genome],hsp.query)
                tan = self.tcdbhmg(self.mytcdb[tid],hsp.sbjct)
                self.what(hsp.query,hsp.sbjct,self.indir+"/img/"+genome+".png",[san,tan])
            except:
                print hsp.query
                print hsp.sbjct
                print genome
                print 'error, quit'
                quit()

        (family,tcid,acc) = ParseTC(ParseDefline(tcdb.title,True).description)
        try:
            query_tms = len(self.tms['queries'][genome])
        except:
            query_tms = 0
        try:
            hit_tms = len(self.tms['tcdb'][ParseDefline(tcdb.title,True).id])
        except:
            hit_tms = 0
        self.globalcount += 1
        if self.ortho is False:
            family = self.names.get_family_abr('.'.join(tcid.split('.')[0:3]))
        else:
            family = 'family_place_holder'


        '''
        Edits made by Vasu Pranav Sai Iddamsetty (VI)
        
        -Adding another file that is a tsv(tab seperated values) file so that the output of gblast
         may be parsed by another program easily.
         
         ****************************************************
         
                This is where the substrates are found.
                They populate the 'mysubstrate' variable.        
         
         ****************************************************
        '''
        ident = round((float(hsp.identities)/len(hsp.match))*100)

        '''
        substrate_info= self.substrates.get_tcid_substrates(tcid)
        mysubstrate_html = ''
        mysubstrate_tsv = ''
        
        if substrate_info is not None:
            
            category,gen_sub,spec_sub,chebi_id,status = substrate_info
            
            chebi_info = chebi_id.replace('"','').replace(' ','').split(',')
            chebi = ''

            for i in chebi_info:
            
                chebi += ' <a href="https://www.ebi.ac.uk/chebi/searchId.do;jsessionid=A9D16DCB24C6F74339FC28A4941EFBB6?chebiId=CHEBI:{}">{}</a>'.format(i,i)
        
            mysubstrate_html = '{},{},{},{},{}'.format(category,gen_sub,spec_sub,chebi,status)
        
            mysubstrate_tsv = ",".join(substrate_info)
        
        else:
        
            mysubstrate_html = 'None'
            mysubstrate_tsv = 'None'
        '''

        substrate = ''

        try:

            substrate = self.substrates[tcid]

        except:

            substrate = 'None'

        '''
        html_results (VI)
        '''

        glink = '<a href="content.html#%s">%s</a>'%(genome,genome)
        tclink = '<a href="http://tcdb.org/search/result.php?tc={}">{}</a>'.format(tcid,tcid)
        h_row = (glink,acc,tclink,tcdb.title,len(hsp.match),tcdb.e,ident,q_len,h_len,qcov,hcov,query_tms,\
        hit_tms,overlap,family,substrate,self.globalcount)
        h_row =['<td>'+str(i)+'</td>' for i in h_row]
        htmlrow = "<tr>\n%s\n</tr>"%("\n\t".join(h_row))


        '''
        text results (VI)
        '''

        row = (genome,acc,tcid,tcdb.title,len(hsp.match),tcdb.e,ident,q_len,h_len,qcov,hcov,query_tms,hit_tms,overlap,family,substrate,self.globalcount)
        row =[str(i) for i in row]
        txtrow = "%s\n"%("\t".join(row))

        '''
        content results (VI)
        '''

        self.plotHydro(genome,tcdb,acc,hsp)

        #self.rows.append(htmlrow)
        if self.cdd_on is True:
            mycdd = 'Query & TC-Hit Conserved Domains:<br><img src=\'cdd/%s.1.png\'><br><img src=\'cdd/%s.2.png\'><br>'%(genome,genome)
        else:
            mycdd = ''
        ol = float(overlap) if overlap is not None else float(0)
        content = '''<div class='result' id='%s'> <h3><a name='%s'>%s</a></h3>  <p>Hit Accession: %s<br>   Hit TCID: %s</p> <p>Hit Description: %s<br>
        <br>   Mach Len: %i<br>   e:%f</p> <p>Query TMS Count : %i<br>   Hit TMS Count: %i     <br>     TMS-Overlap Score: %f<br>
        Predicted Substrates:%s <br><br>     BLAST Alignment:<br>     <pre>     %s     </pre> <br>  <table><tr><th>Protein Hydropathy Plots:</th></tr> 
        <tr><td><img src='img/%s_hydro.png'></td> <td><img src='img/%s-%s_hydro.png'></td></tr><br>
        <tr><th><br> Pairwise Alignment-Hydropathy Plot:<br></th></tr>
        <tr><td colspan="2" style="text-align: center;"><img src='img/%s.png'></td></tr></table><br>%s </p> </div>\n''' %(genome,genome,genome,acc,tcid,tcdb.title,\
        len(hsp.match),tcdb.e,query_tms,hit_tms,ol,substrate,str(hsp),genome,acc,genome,urlencode(genome),mycdd)
        #self.data.append(content)
        return htmlrow,content,txtrow

    def write_results(self):
        bestresults = copy.deepcopy(self.bestresults)
        bestresults.extend(self.notmsresults)
        
        #create and open the necessary files for results (VI)
        html_results = open(self.indir+'/results.html','wb')
        content = open(self.indir+'/content.html','wb')
        tsv_results = open(self.indir+'/results.tsv', 'wb')
        
        '''
            
        We are unsure of why the local number variable( the one that is supposed to count which result it is) is the last number in the rows.  --Vasu Pranav Sai Iddamsetty
        '''
        
        
        #write to results.html and content.html (VI)
        html = '''<html><table width="100%%" border="1"> <tr> <td>Query ID</td> <td>Hit ID</td>
        <td>Hit TCID</td> <td>Hit Description</td> <td>Match Len</td> <td>e-Val</td> <td>% Identity</td> <td>Query Length</td> <td>Hit Length</td> <td>Query Coverage</td> <td>Hit Coverage</td> <td>Query TMS#</td>
        <td>Hit TMS#</td> <td>TM-Overlap Score</td> <td>Family Abrv.</td><td>Predicted Substrate</td> <td> #</td> </tr><br>\n\n'''
        html_results.write(html)
        content.write("<html>")
        
        
        #write the column descriptions to results.tsv (VI)
        columnDescription = '''#Query_id\tHit_xid\tHit_tcid\tHit_desc\tMatch_length\te-value\t%_identity\tQuery_Length\tHit_Length\tQuery_Coverage\tHit_Coverage\tQuery_n_TMS\tHit_n_TMS\tTM_Overlap_Score\tFamily_Abrv\tPredicted_Substrate\trow_number\n'''
        tsv_results.write(columnDescription)
        
        
        
        self.myqueries = SeqIO.parse(self.indir+'/myqueries.faa','fasta')
        self.mytcdb    = SeqIO.parse(self.indir+'/mytcdb.faa','fasta')
        self.myqueries = SeqIO.to_dict(self.myqueries)
        self.mytcdb    = SeqIO.to_dict(self.mytcdb)
        if self.cdd_on:
            self.cdd_extract()
        self.queryhmg = hmmgap.annotate()
        self.tcdbhmg  = hmmgap.annotate()
        self.queryhmg.hmmtop = self.tms['queries']
        self.tcdbhmg.hmmtop  = self.tms['tcdb']
        for res in bestresults:
            
            #retrieve relevant formaatted data and write it to the files (VI)
            (row,data,txt) = self.build_view(res)
            html_results.write(row)
            content.write(data)
            tsv_results.write(txt)
            
            print "Generated Results for :: %s" %ParseDefline(res[0]).id
        
        #end tags for the .html files
        html_results.write("</table></html>")
        content.write("</html>")

    def loadSubstrates(self):

        print('Loading Substrates')

        substrateData = urlopen('http://tcdb.org/cgi-bin/projectv/getSubstrates.py')


        for line in substrateData:

            data = line.replace('\n','').split('\t')
            #print(data)
            self.substrates[data[0]] = data[1].replace('|',', ')


        return
                         

    def cdd_extract(self):
        if os.path.exists(self.indir+'/cdd') is False:
            os.mkdir(self.indir+'/cdd')
        fastas  = dict(self.myqueries,**self.mytcdb)
        thisdir = self.indir+'/cdd/'
        for query,hit,hsp in self.goodresults:
            if os.path.exists(thisdir+query.title()+'.1.png') is False:
                cdd.fetch(str(hsp.query),thisdir+query.title()+'.1.png')
            if os.path.exists(thisdir+query.title()+'.2.png') is False:
                cdd.fetch(str(hsp.sbjct),thisdir+query.title()+'.2.png')

    '''
    Plots individual Protein hydropathies (VI)
    '''

    def plotHydro(self,genome,hit,acc,hsp):

        #File Paths

        queryPath = self.indir+'/img/'+genome+'_hydro.png'
        hitPath = self.indir+'/img/'+acc+'-'+genome+'_hydro.png'

        quod = ' -s -q -d {} --width 10 --height 3 --xticks 50'.format(self.indir+"/img/")


        #Query Hydropathy
        if not os.path.exists(queryPath):

            query=ParseDefline(genome).id
            querySeq = self.queries[str(query)].seq
            query = 'quod.py {} -o {} -c blue -w {}-{} -l {}'.format(querySeq,genome+'_hydro',hsp.query_start,hsp.query_end,genome) + quod

            #os.system(query)
            subprocess.call(query.split())

        #Hit Hydropathy
        if not os.path.exists(hitPath):
        
            hitID = ParseDefline(hit.title,True).id
            hitSeq = self.tcdbHits[str(hitID)].seq
            hit = 'quod.py {} -o {} -c red -w {}-{} -l {}'.format(hitSeq, acc+'-'+genome+'_hydro',hsp.sbjct_start,hsp.sbjct_end,acc) + quod

            #os.system(hit)
            subprocess.call(hit.split())

    '''
    Within the function hydro, we added a (0,0) value for the '*' symbols, so that 
    hydropathies can still be visualized for pseudogenes. - Vasu Pranav Sai Iddamsetty
    '''

    def hydro(self,gseq):
        seq=gseq.replace('-','')
        window = 19
        prev = 0
        index = {'G':(-0.400,0.48),
                 'I':(4.500,1.38),
                 'S':(-0.800,-0.18),
                 'Q':(-3.500,-0.85),
                 'E':(-3.500,-0.74),
                 'A':(1.800,0.62),
                 'M':(1.900,0.64),
                 'T':(-0.700,-0.05),
                 'Y':(-1.300,0.26),
                 'H':(-3.200,-0.4),
                 'V':(4.200,1.08),
                 'F':(2.800,1.19),
                 'C':(2.500,0.29),
                 'W':(-0.900,0.81),
                 'K':(-3.900,-1.5),
                 'L':(3.800,1.06),
                 'P':(-1.600,0.12),
                 'N':(-3.500,-0.78),
                 'D':(-3.500,-0.90),
                 'R':(-4.500,-2.53),
                 'X':(0,0),
                 'U':(0,0),
                 '*':(0,0)}
        midpt = (window+1)/2
        length = len(seq)
        hydro = []
        for i in range(length-window+1):
            total = 0
            for j in range(window):
                total +=index[seq[i+j]][0]
            total = total/window
            hydro.append(total)
        if len(seq) == len(gseq):
            return hydro
        replace = re.finditer('(-+)',gseq)
        inserts = {}
        for i in replace:
            inserts.setdefault(i.start(),i.end()-i.start())
        first = False
        newhydro = []
        for x, h in enumerate(hydro):
            if x in inserts.keys() and first is False:
                first = True
                for y in range(inserts[x]):
                    newhydro.append(0)
                newcount = x + inserts[x]
                continue
            if first is False:
                newhydro.append(h)
                continue
            if first is True and newcount in inserts.keys():
                for y in range(inserts[newcount]):
                    newhydro.append(0)
                newcount += inserts[newcount]
                continue
            else:
                newhydro.append(h)
                newcount +=1

        return newhydro


    def what(self,a,b,outfile,hmt=[]):
        #quit()
        ha = self.hydro(a)
        hb = self.hydro(b)
        omg = [len(ha),len(hb)]
        readbar = re.compile(r'\d+?[+_]+\+')
        omg.sort()
        ha =ha[0:omg[0]]
        hb =hb[0:omg[0]]
        x_data=range(0,len(ha))
        plt.figure()
        plt.axhline(y=0,color='black')
        plt.ylim(-3,3)
        plt.xlim(right=len(a))
        plt.plot(x_data,ha,linewidth=1,label=self.alabel,color='blue')
        plt.plot(x_data,hb,linewidth=1,label=self.blabel,color='red')
        plt.xlabel("Residue")
        plt.ylabel("Hydropathy (kcal/mol)")
        plt.legend(loc='lower right')
        # Draw TMS bars
        if len(hmt) == 2:
            sub = readbar.finditer(str(hmt[0]))
            tar = readbar.finditer(str(hmt[1]))
            for tms in sub:
                subtract = (len(hmt[0])-len(ha))/2
                if tms.end()-subtract>len(ha):
                    subtract = tms.end()-len(ha)
                plt.axvspan(tms.start()-subtract,tms.end()-subtract, facecolor="blue", alpha=0.3)
            for tms in tar:
                subtract = (len(hmt[1])-len(ha))/2
                if tms.end()-subtract>len(hb):
                    subtract = tms.end()-len(hb)
                plt.axvspan(tms.start()-subtract,tms.end()-subtract, facecolor="red", alpha=0.3)
        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(10,3)
        plt.savefig(outfile, dpi=100, format="png",bbox_inches='tight', pad_inches=0.003)
        plt.clf()
        plt.close(fig)



if __name__=="__main__":
    from optparse import OptionParser,OptionGroup
    desc = "Welcome to GBlast! Easily identify transporters in entire Genomes/Proteomes - By Vamsee Reddy"
    version = "GBlast V3"
    opts = OptionParser(description=desc,version=version)
    opts.add_option('-i',
                    action='store',
                    type='string',
                    dest='input',
                    help="Path to genome/proteome file"
    )
    opts.add_option('-o',
                    action='store',
                    type='string',
                    dest='output',default='genome.fsa',
                    help="Results output name"
    )
    opts.add_option('--evalue',
                    action='store',
                    type='float',
                    dest='evalue',
                    default=0.001,
                    help="Minimum e-Value [0.001]"
    )
    opts.add_option('--cov',
                    action='store',
                    type='float',
                    dest='mincov',
                    default=40.0,
                    help="Minimum Proten Coverage [40.0]"
    )
    opts.add_option('--esort',
                    action='store_true',
                    dest='esort',
                    default=False,
                    help="Use e-value as the preliminary criteria for sorting (otherwise sorted by TCID)"
    )
    '''    
    opts.add_option('--betabarrel',
                    action='store_true',
                    dest='bb',
                    default=False,
                    help="Find Beta Barrels instead of TMS"
    )
    opts.add_option('--cdd',
                    action='store_true',
                    default=False,
                    dest='cdd',
                    help="Include CDD analysis (Takes a while)"
    )
    opts.add_option('--orthologs',
                    action='store',
                    dest='ortho',
                    default=False,
                    help="Find orthologs with this fasta file (optional)"
    )
    opts.add_option('--query_gi',
                    action='store',
                    dest='subgi',
                    default=False,
                    help="Ortholog search using list of only these queries"
    )
    opts.add_option('--target_gi',
                    action='store',
                    dest='targi',
                    default=False,
                    help="Ortholog search using list of only these targets"
    )
    '''
    (cli,args)=opts.parse_args()
    if cli.input is not None and cli.output is not None:
        GB = Tools()
        #if(cli.bb):
            #GB.dbfile = '/db/betabarrel'
        GB.ortho  = False
        GB.indir  = cli.output
        GB.cdd_on = False
        GB.query  = cli.input
        GB.query_gis  = False
        GB.target_gis = False
        GB.expect = cli.evalue
        GB.mincov = cli.mincov
        GB.esort=cli.esort
        GB.blast_all()
        print "Loading BLAST results"
        GB.load_good()
        print "Writing FASTAS"
        GB.write_fastas()
        print "Running HMMTOP"
        GB.hmmtop()
        print "Calculating TMS Overlap"
        GB.calculate_tms_scores()
        print "Writing Results"
        GB.write_results()
    else:
        opts.print_help()


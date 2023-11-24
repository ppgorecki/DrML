#!/usr/bin/python




DRMLVER="1.00"

# (c) Copyright 2009-2012 by Pawel Gorecki & Oliver Eulenstein
#
# Permission is granted to copy and use this program provided no fee is
# charged for it and provided that this copyright notice is not removed.
#
# Please send bug reports, comments etc. to: gorecki@mimuw.edu.pl
#
# Last update: December 17, 2012
#

import string, types, math, random, sys, getopt, os, time
from functools import reduce, cmp_to_key

lmbd = 0.005

VERBOSEQUIET=-1

VERBOSEDLS=1
VERBOSEDS=2
VERBOSEGEN=3     # with kappas from DP

# DEBUGs
debug=0
DDP=1    # dp debug
DSP=2    # seq pair debug
DDS=4    # ds debug
DMLR=8   # mlr debug
DEX=16   # exhaustive ml debug
DSPX=32   # seq pair ext
DSPX2=64   # seq pair ext2 (unused)

verbose=VERBOSEDLS
verbosehi=VERBOSEQUIET       # verbose level for hard instances

BBSINGLETONS=1
BBPAIR=2
BBSCORES=4
BBRAND=8
BBRANDX=16

DPTYPE=1

TPREFIX=-1
TINFIX=0
TPOSTFIX=1

TLEAF=1
TINT=2
TALL=3

SOLVEML=1
SOLVEDS=2
SOLVEDLS=3
SOLVEHI=4

dupspecsettingcnt=0

basename="drml"
gsmatch=None

dlstypes=['+','-','~']

def usage():
    print("""
Drml version %s.

Drml generates at most two files:
drml.log - log file
drml.dls - output file that contains:
  gene tree (from the input)
  species tree (from the input)
  dlstree(s)
A DLS-tree is a tree, representing evolutionary scenario, written in the nested parenthesis notation 
with additional decoration: + denotes a duplication node, ~ denotes a speciation node, 
- denotes a loss node.

Both file names can be changed with -o or -p options. See below for more details.

Options:

1. Defining gene and species trees.

-s TREE - defines a rooted species tree with branch lengths
-g TREE - defines rooted gene tree

-P FILE - input file containing a pair of gene and species trees (# starts a comment):
    GENE TREE
    SPECIES TREE

-p FILE - like above, but the log file is FILE.log instead of drml.log

2. Random trees

-r NUM - generate NUM random pairs of trees
-G NUM - number of leaves in random gene tree
-S NUM - number of leaves in random species tree
-l NUM - define max length of species tree branches (for random generator)

3. General options

-v NUM - verbose level
    0 - summary only
    1 - summary and scenarios (DLS trees)
    2 - summary, DS-settings and DLS trees
    3 - summary, dynammic programming output, DS-settings and DLS trees

-V NUM - verbose level for hard instances processing

-o BASENAME
    Defines basename. Default is drml. If -p FILE[.txt] option is given then the basename is FILE.

-R print pairs of input trees (gene tree, species trees,...)

3. Maximum likelihood estimation

-m  - run ML estimation

-d LEV
    Estimation type (use with "-m"):
     m - compute ML (DP only), no hard instances processing (default)
     d - compute duplication-speciation settings, no hard instances processing
     ra - compute all optimal reconciliations (DLS trees), no hard instances processing
     r - compute one optimal reconciliation (DLS trees), no hard instances processing
     ha - full, including hard instances, all optimal reconciliations
     h - full, show only one optimal reconciliation

-b ALG
    Type of algorithm for solving hard instances (use with "-m -dh" or "-m -dha"):
     1 - all singletons first
     2 - pairs
     5 - singletons+priority queue
     6 - pairs+priority queue (default)

-L - duplication rate, lambda for Poisson distribution

4. Exhaustive estimation and comparing results.

-n - ML by exhaustive enumeration, for testing on small trees only

5. Gene-species tree labels matching.
   
-M type - how the labels of the leaves in the gene tree determine a species name 

This option should be used when the labels in gene trees are different from the species names.

default (no -M) - exact match 
pNUM[,LEN] - where NUM>0 - species names start at NUM position in each gene label, LEN (optional) is the length of species names
p-NUM[,LEN] - where NUM>0 - the position is determined from the end of gene labels
aDELIM[,LEN] - DELIM is a string, species names start after DELIM
bDELIM[,LEN] - species names are before DELIM

Example:

drml.py -s "(a,b)" -g "(a1,(b1,a2))" -m -Mp0,1

6. Processing dls files.

-Df - print frequencies of events (duplications, losses and speciations)
-Dm - compute ml values for each DLS tree
-Dd - print dls trees only 

Note that the species tree should be the same for all dls files.

Example:

drml.py -Df *.dls

""" % DRMLVER)

# Cached values of log factorial
logfact=[ math.log(1), math.log(1) ] # logfactorial

def computelogfact(n):
    while n>len(logfact):
        logfact.append(logfact[-1]+math.log(len(logfact)))
            
def logP(tau,d):
    """ Compute log of Poisson distribution """   
    return -lmbd*float(tau)+d*math.log(lmbd*float(tau))-logfact[d]    

def nextcombination(c, k, n):
    if not k: return 0
    i = k - 1
    c[i]+=1
    while i >= 0 and c[i] >= n-k+1+i:
        i-=1
        c[i]+=1

    if c[0] > n - k: return 0

    for i in range(i + 1, k): c[i] = c[i-1] + 1
    return c


class Node:
    def __init__(self):
        self.v=None
        self.extlab=''

    def isleaf(self):
        return 0

    def parse(cls,s):
        if not len(s): raise NameError('Parsing tree: empty label')
        if s[0]=='(':
            (t1,s)=Node.parse(s[1:])
            if s[0]!=',': raise NameError('Expected , in ' + s.__repr__())
            (t2,s)=Node.parse(s[1:])
            if s[0]!=')': raise NameError('Expected ) in ' + s.__repr__())            
            t=Int(t1,t2)
        else:
            t=Leaf(s[0])

        s=s[1:]
        if s and s[0]==':':
            t.setv(float(s[1]))
            s=s[2:]
        else: t.setv(1.0)

        if s and s[0]=='[':
            t.setextlab(s[1])
            s=s[3:]
        else: t.setextlab("")
        
        t.addlabels=[]
        while s and s[0] not in ",)":
            t.addlabels.append(s[0])
            s=s[1:]
        return (t,s)

    parse= classmethod(parse)

    def __repr__(self):
        return self.__str__()

    def setv(self,v):
        self.v=v

    def setextlab(self,l):
        self.extlab=l

    def __str__(self):
        if self.v and self.v!=1.0: return ':'+str(self.v)
        return ''

    def setparent(self,par):
        if not par: self.height=0
        else: self.height=par.height+1
        self.parent=par
        

    #(leaves) <--- self <---- x <--- (root)
    def isdescendantof(self,x):
        par=self
        while par:
            if par==x: return 1
            par=par.parent
        return 0

    def lca(self,x):
        if x==self: return self
        y=self
        if y.height>x.height: x,y=y,x
        while y.height<x.height: x=x.parent
        while x!=y:
            x=x.parent
            y=y.parent
        return x

    def printdr(self):
        return ':%s [&&DR dup=%d spec=%d kappa=%d phi=%d]' % (self.v,self.delta,
                                                              self.sigma,self.kappa,self.phi)

    def randombranchlenghts(self,a,b):
        self.v=random.randrange(a,b)

   


class Leaf(Node):
    def __init__(self,label):
        Node.__init__(self)
        self.label=label
        self._internalpostfix=[self]

    def drmlcopy(self):
        return Leaf(self.label)

    def isleaf(self):
        return 1

    def _nodes(self, dir, type):  ## Leaf class
        if type & TLEAF: return [ self ]
        return []

    def label(self):
        return self.label

    def __str__(self):
        if self.extlab: return self.extlab
        return str(self.label)+Node.__str__(self)

    def internal(self): return []

    def setpostfixnum(self,n):
        self.num=n
        self.minnum=n
        return (n,n+1)

    def setsize(self): self.size=1

    def printdr(self): return self.label+Node.printdr(self)

    def printsmp(self): return self.label

    def labs(self): return self.label

    def sortnodes(self): return self.label

    def setdls(self):
        if self.label[-1]=='-':
            self.label=self.label[0:-1]+" -"

        self.cluster=self.label.split()
        if '-' in self.cluster:
            self.cluster.remove('-')
            self.dlstype='-'
        else:
            self.dlstype='~'  # this is a special speciation-gene node
        self.cluster=set(gsmatchtransform(self.cluster))

#end of Leaf class

class Int(Node):
    def __init__(self,l,r):
        Node.__init__(self)
        self.l=l
        self.r=r
        self._internalpostfix=self._nodes(TPOSTFIX,TINT)  ## Int class

    def drmlcopy(self):
        return Int(self.l.drmlcopy(),self.r.drmlcopy())

    def _nodes(self, dir, type):
        if not (type & TINT): return self.l._nodes(dir,type) + self.r._nodes(dir,type)
        if dir==TPREFIX: return [ self ] + self.l._nodes(dir,type) + self.r._nodes(dir,type)
        if dir==TPOSTFIX: return self.l._nodes(dir,type) + self.r._nodes(dir,type) + [ self ]
        return self.l._nodes(dir,type) + [ self ] + self.r._nodes(dir,type)

    def __str__(self):
        if self.extlab: return self.extlab
        return '(' + self.l.__str__() + ',' + self.r.__str__() + ')'+Node.__str__(self)

    def setparent(self,par):
        Node.setparent(self, par)
        self.l.setparent(self)
        self.r.setparent(self)

    def setpostfixnum(self,n):
        (minnum,n)=self.l.setpostfixnum(n)
        (x,n)=self.r.setpostfixnum(n)
        self.minnum=minnum
        self.num=n
        return (minnum,n+1)

    def setsize(self):
        self.l.setsize()
        self.r.setsize()
        self.size=1+self.l.size+self.r.size

    def printdr(self):
        return '('+self.l.printdr()+','+self.r.printdr()+')'+Node.printdr(self)

    def printsmp(self):
        if self.extlab: return self.extlab
        return '('+self.l.printsmp()+','+self.r.printsmp()+')'

    def randombranchlenghts(self,a,b):
        self.l.randombranchlenghts(a,b)
        self.r.randombranchlenghts(a,b)
        Node.randombranchlenghts(self, a, b)

    def labs(self):
        return self.l.labs()+" "+self.r.labs()

    def sortnodes(self):
        a=self.l.sortnodes()
        b=self.r.sortnodes()
        if a>b:
            x=self.r
            self.r=self.l
            self.l=x
            return b
        return a

    def setdls(self):
        self.l.setdls()
        self.r.setdls()
        self.cluster=self.l.cluster.union(self.r.cluster)
        # duplication            
        if self.r.cluster==self.l.cluster:
            self.dlstype='+'
        elif not self.r.cluster.intersection(self.l.cluster):
            self.dlstype='~'
        else: 
            print("Unknown type of internal node: %s, %s" % (self.l,self.r))
            print("Speciation or duplication expected")
            sys.exit(-1)

#end of Int class

class DPstats:

    def __init__(self,logml):
        self.logml=logml
        self.dls=0
        self.ds=0
        self.ds0=0

    def __str__(self):
        return "%s %d %d %d" % (self.logml.__repr__(),self.ds,self.ds0,self.dls)


class Drml:


    def __init__(self,basename,solve,single,bbalg):
        self.logfile=None
        self.stats=DPstats(None)     # global stats
        self.dpstats=DPstats(None)   # dp stats
        self.solve=solve
        self.bbalg=bbalg
        self.hi=0            # hard instance, 0 - no (first run), 1 - DP in hi
        self.himode=0        # computing DP for hard instances if single==0, using historage
        self.bbsteps=0       # bbsteps in hard instance processing
        self.initbasename(basename)
        
        self.single=single
        self.historage=[]    # storage of copies of optimal gt,gt trees, used when himode==1
        self.histart=0       # executed hi processing
        self.starttime=None

    def hiprocessing(self):
        if not self.single:
            self.single=1
            self.himode=1
        self.histart=1

    def stophiprocessing(self):
        self.single=0
        self.himode=0
        self.dpstats=DPstats(None)
        self.stats=DPstats(None)

    def initbasename(self,basename):        
        self.basename=basename
        if self.logfile: 
            self.logfile.close()
        self.logfile=open(basename+".log","w")
        self.logfile.write("DrML ver. " + DRMLVER + "\n")

    def start(self,gt,st):
        self.gt=gt
        self.st=st
        self.starttime=time.time()
        self.starttimeinfo=time.strftime("Start at %a %H:%M %m/%d/%y\n", time.localtime())

        self.log("# Gene tree\n" + self.gt.__repr__() + "\n")
        self.log("# Species tree\n" + self.st.__repr__() + "\n")




    def __del__(self):
        if self.logfile:
            if self.starttime:
                self.log(self.summary())
            self.log("\nLog written to: %s" % self.logfile.name)                     
            self.logfile.close()

    def stop(self):
        self.totaltime=time.time()-self.starttime

    def initdpstats(self,logml):
        self.dpstats=DPstats(logml)

    def savegtst(self):
        self.curdlsfile=self.basename+".dls"
        
        f=open(self.curdlsfile,'w')
        f.write("# Gene tree\n" + self.gt.__repr__() + "\n")
        f.write("# Species tree\n" + self.st.__repr__() + "\n")
        f.close()

    def newscenario(self,dlstree,dupnum):
        self.dpstats.dls+=1
        if self.stats.logml is None or self.stats.logml < self.dpstats.logml:
            # initialize global, new better scenario
            self.savegtst()
            self.stats.dls=0
            self.stats.ds=self.dpstats.ds
            self.stats.ds0=self.dpstats.ds0
            self.stats.logml=self.dpstats.logml
            self.stats.dupnum=dupnum
            if self.himode==1:
                # initialize historage
                self.historage=[]

        if self.dpstats.logml==self.stats.logml:
            self.stats.dls+=1
            f=open(self.curdlsfile,'a')
            if self.stats.dls==1:
                f.write("""#
#For more detailed theory of DLS trees read paper
#P.Gorecki and J.Tiuryn, DLS-trees: a model of evolutionary scenarios, TCS, Vol. 359, 2006
#Here:
#  + denotes duplication node
#  ~ denotes speciation node
#  - denotes loss node
#
""")
            f.write("# DlsTree no %d\n# logml=%f\n#dupnum=%d\n%s\n" % (self.stats.dls,self.stats.logml,dupnum,dlstree))

            if self.himode==1:
                # store currrent curlock setting
                # copy DP results
                higt=self.gt.drmlcopy()
                hist=self.st.drmlcopy()
                hist.setpostfixnum()
                lcamapping(higt, hist)

                for (s,d) in zip(self.gt.prefall,higt.prefall):
                    if not d.isleaf(): d.spec = s.spec
                    d.lock = s.lock
                    d.ints = s.ints

                for (s,d) in zip(self.st.prefall,hist.prefall):
                    d.phi = s.phi
                    d.speclocked = s.speclocked
                    if s.isleaf():
                        d.ML = s.ML
                    else:
                        d.MLreconstr = s.MLreconstr
                        d.sumlockeddual = s.sumlockeddual
                        initializeInternalDP(higt, d)

                self.historage.append((higt,hist))
                print (self.historage)

            f.close()


    def newdssetting(self):
        self.dpstats.ds+=1
        if self.dpstats.logml==self.stats.logml:
            self.stats.ds+=1

    def newds0setting(self):
        self.dpstats.ds0+=1
        if self.dpstats.logml==self.stats.logml:
            self.stats.ds0+=1


    def summary(self):

        s="\nDrML ver. %s - execution summary\n" % DRMLVER

        lcaml,lcadup=lcamlvalue(self.gt,self.st)
        fatml,fatdup=fatmlvalue(self.gt,self.st)

        s+="logML for the optimal scenario: %f\n" % ( self.stats.logml if self.stats.logml else self.dpstats.logml )    
        s+="logL for the lca-scenario: %f\n" % lcaml
        s+="logL for the fat-scenario: %f\n" % fatml

        if hasattr(self.stats,'dupnum'):
            s+="Number of gene duplications in the first optimal scenario: %d\n" % self.stats.dupnum
        s+="Number of gene duplications in the lca-scenario: %d\n" % lcadup
        s+="Number of gene duplications in the fat-scenario: %d\n" % fatdup
        
        s+="lambda (duplication rate): %f\n" % ( lmbd )
        s+="DrML mode: "
        if self.solve==SOLVEML:
            s+="compute logml only (no hard instance detecting)"
        elif self.solve==SOLVEDS:
            s+="logml, DS settings (no hard instance detecting)"
        elif self.solve==SOLVEDLS:
            if self.single:
                s+="logml, DS settings, one optimal reconciliation (no hard instance detecting)"
            else:
                s+="logml, DS settings, all optimal reconciliations (no hard instance detecting)"
        elif self.solve==SOLVEHI:
            if self.single:
                s+="(full) logml, single optimal reconciliation, resolve hard instances"
            else:
                s+="(full) logml, all optimal reconciliations, resolve hard instances"
            s+="\nBB algorithm:"
            if self.bbalg & BBSINGLETONS: s+="singletons first"
            if self.bbalg & BBPAIR : s+="by pairs"
            if self.bbalg & BBRAND: s+="by random"
            if self.bbalg & BBSCORES : s+=" and use priorities"
            s+="\n"
        else:
            s+=" unknown option?"

        if self.stats.dls:
            s+="\nDS settings: %d" % self.stats.ds
            s+="\nDS settings without reconciliation: %d" % self.stats.ds0
            s+="\nReconciliations-DLS: %d" %self.stats.dls
            if self.solve>=SOLVEHI:
                if self.histart:
                    s+="\nHard instance: Yes"
                    s+="\nBB-steps: %d" % self.bbsteps
                    if not self.single:
                        s+="\nNumber of optimal raised/locked configurations: %d" % len(self.historage)
                else:
                    s+="\nHard instance: No"
            else:
                s+="\nHard instance: No"
        else:
            if self.solve>=SOLVEDS:
                s+="\nDS settings: %d" % self.dpstats.ds
                s+="\nDS settings without reconciliation: %d" % self.dpstats.ds0
            if self.solve>=SOLVEDLS:
                 s+="\nReconciliations-DLS: %d" %self.stats.dls

            if self.solve<=SOLVEDS:
                s+="\nHard instance: Unknown"
                s+="\nNote: run with -dr or higher option to check the correctness of MLE"
            else:
                s+="\nHard instance: Yes"
                s+="\nNote: this is hard instance - run with -dh or -dha option to compute a valid MLE value"

        s+="\nTime: %.2f" % self.totaltime
        if self.stats.dls>0:
            s+="\nScenarios saved to: %s\n" % self.curdlsfile

        return s

    def timelog(self):
        s="curlogML:"
        if self.stats.logml: s+="%f" % self.stats.logml
        else: s+="None"
        self.log("[File %s time %d %s]\n" % (self.basename,time.time()-self.starttime,s))

    def log(self,*par):
        for i in par: print (i,end='')        
        for i in par: print (i,end='',file=self.logfile)
        


class Tree:

    def __init__(self,s):
        if isinstance(s,Node): self.root=s
        if type(s) is str: self.root=Node.parse(self.tokenize(s))[0]
        self.root.setparent(None)
        self.allleaves=self.root._nodes(0,TLEAF)
        self.prefall=self.nodes()              
        self.postint=self.internal()
        self.postall=self.nodes(TPOSTFIX)

    def __str__(self):
        return self.root.__str__()

    def __repr__(self):
        return self.root.__repr__()

    def internal(self):
        return self.root._nodes(TPOSTFIX,TINT)

    def nodes(self, dir=TPREFIX,type=TALL):
        return self.root._nodes(dir,type)

    def setsize(self):
        self.root.setsize()

    # tree uniqueness!
    def lab2nodes(self):
        return dict(zip( (x.label for x in self.allleaves),self.allleaves ))

    def lca(self,a,b):
        return a.lca(b)

    def setpostfixnum(self):
        self.root.setpostfixnum(0)

    def printdr(self):
        return self.root.printdr()

    def randombranchlenghts(self,a,b):
        self.root.randombranchlenghts(a,b)

    def sortnodes(self):
        self.root.sortnodes()

    def tokenize(self,s):
        """ Tree parser utility """
        ls=[]
        cs=''
        for i in s:
            if i in ')(,:[]':
                ls.append(cs.lstrip().rstrip())
                ls.append(i)
                cs=''
            else:
                cs=cs+i
        ls.append(cs.lstrip().rstrip())
        return [ s for s in ls if len(s)>0 ]

    def drmlcopy(self):
        return Tree(self.root.drmlcopy())


class DlsTree(Tree):
    def __init__(self,s):
        Tree.__init__(self,s)

        # postprocessing
        # cluster reconstruction

        self.root.setdls()
        #for i in self.nodes():
        #    print i.dlstype,i,i.cluster,i.addlabels

# cluster is a tuple of strings 
def gsmatchtransform(cluster):
    global gsmatch
    

    if not gsmatch: return cluster
    tp=gsmatch[0]
    par=gsmatch[1:]
    ln=0
    if ',' in par:
        par,ln=par.rsplit(',',1)        
        ln=int(ln)
    if tp=='p':
        pos=int(par)
        s=( c[pos:] for c in cluster)

    elif tp=='a':
        s=( c.split(par,1)[1] for c in cluster)
    elif tp=='b':
        s=( c.split(par,1)[0] for c in cluster)

    if ln: s=( c[:ln] for c in s)
        
    return s


def randomtreefromlabels(lableaves):
    """ Return random Tree from leaves """
    n = [ Leaf(l) for l in lableaves]
    while len(n)>1:
        #print random.randint(0,len(n)), len(n)
        l=n.pop(random.randint(0,len(n)-1))
        l.extlab=''
        r=n.pop(random.randint(0,len(n)-1))
        r.extlab=''
        n.append(Int(l,r))

    n[0].extlab=''
    return Tree(n[0])


def randomtree(speclist,unique=0,numleaves=-1):
    """
    Generate random tree
    speclist -- list of species names
    unique -- generate a tree with unique leaf labels
    numleaves -- number of leaves in the output tree (opt)
    """
    n=[]
    if unique and numleaves<=0: numleaves=len(speclist)

    for i in range(numleaves):
        if unique:
            l=speclist.pop(random.randint(0,len(speclist)-1))
        else: l=speclist[random.randint(0,len(speclist)-1)]
        n.append(l)
    return randomtreefromlabels(n)

def lcamapping(gt,st):
    """ Compute lca mapping from gt to st """
    stdict=st.lab2nodes()
    for g in gt.postall:
        if g.isleaf(): g.lcamapping=stdict[ list(gsmatchtransform( (g.label,) ))[0] ]
        else: g.lcamapping=st.lca(g.l.lcamapping,g.r.lcamapping)

def setlcadupspec(gt):
    """ Set lcadup lcaspec attributes in the nodes of gt """
    for g in gt.allleaves:
        g.lcadup=0
        g.lcaspec=0

    for g in gt.postint:
        g.lcadup=g.lcamapping==g.r.lcamapping or g.lcamapping==g.l.lcamapping
        g.lcaspec=not g.lcadup


################################################################
## SEQ PAIR

class SeqPairStruct:

    def __init__(self,mseq):
        self.z=0
        self.alfa=0
        self.beta=0
        self.seq=[]
        self.t=[]
        self.t.append([])
        self.sumx=0
        self.sumy=0
        self.curmax=0
        for i in mseq:
            if i==(0,0):
                self.z+=1
            elif i==(1,0):
                self.alfa-=1
                self.z+=1
            elif i==(0,1):
                self.beta-=1
                self.z+=1
            else:
                self.sumx+=i[0]
                self.sumy+=i[1]
                self.seq.append(i)
        #self.seq.sort(cmp=lambda x,y: x[0]+x[1]-y[0]-y[1])

        self.seq = sorted(self.seq, 
            key=cmp_to_key(lambda x,y: x[0]+x[1]-y[0]-y[1]))
            


    def checksmp(self,alfa,beta):       
        if len(self.t)<=alfa:
            self.t.extend( [[]]*(1+alfa-len(self.t)))

        if len(self.t[alfa])<=beta:
            self.t[alfa].extend([None]*(1+beta-len(self.t[alfa])))

        return self.t[alfa][beta]

    def adjalfabeta(self,alfa,beta):
        al=alfa+self.alfa
        be=beta+self.beta
        z=self.z
        if al<0:
            z+=al
            al=0
        if be<0:
            z+=be
            be=0
        return (al,be,z)

    def get(self,alfa,beta):
        
        smp=self.checksmp(alfa, beta)
        if smp!=None: return smp
        self.t[alfa][beta]=self.getdet(alfa,beta)
        return self.t[alfa][beta]

    def getdet(self,alfa,beta):

        
        (al,be,z)=self.adjalfabeta(alfa,beta)

        seq = list(filter( lambda p: p[0]<=al and p[1]<=be, self.seq))

        if not seq:
            if debug & DSP: print ("SMP%d.%d.res=%d" % (al,be,z),end='')
            return z

        self.curmax=0
        self.callind=0
        res = z+self.seqpairm(al, be, seq, 0)

        if debug & DSP:
            print ("sq:%d.%d.res=%d.ci=%d" % (alfa,beta,res,self.callind),'z=',z)

        return res

    def seqpairm(self,alfa, beta, seq, inserted,coall=0):
        self.callind+=1
        first=seq.pop()
        z1=0
        z2=0
        if alfa-first[0]>=0 and beta-first[1]>=0:
            al=alfa-first[0]
            be=beta-first[1]
            seq2=list(filter( lambda p: p[0]<=al and p[1]<=be, seq))
            if not seq2:
                z1 = 1
            else:
                if coall or inserted+1+len(seq)>self.curmax:
                    z1 = 1+self.seqpairm(al,be,seq2,inserted+1)
            self.curmax=max(self.curmax,z1+inserted)

        if seq:
            if coall or inserted+len(seq)>self.curmax:
                z2=self.seqpairm(alfa,beta,seq,inserted)
                self.curmax=max(self.curmax,z2+inserted)

        return max(z1,z2)

    def getmaxdiag(self,alfa,beta):

        msi=-1
        for i in range(alfa+beta+1):
            si=self.get(i,alfa+beta-i)
            if msi<0 or msi<si:
                msi=si
        return msi

class SeqPairStructExt(SeqPairStruct):

    def __init__(self,mseq):
        SeqPairStruct.__init__(self,mseq)
        seqd=dict()
        for i in self.seq:
            if not i in seqd:
                seqd[i]=1
            else:
                seqd[i]=seqd[i]+1
        self.seq = list(map( lambda p: (p[0][0],p[0][1],p[1]), seqd.items()))

        #self.seq.sort(cmp=lambda x,y: x[0]+x[1]-y[0]-y[1])

        self.seq = sorted(self.seq, 
            key=cmp_to_key(lambda x,y: x[0]+x[1]-y[0]-y[1]))


        n = len(self.seq)
        self.sums=list(map( lambda p: [p[0]*p[2],p[1]*p[2]], self.seq))
        self.pairs=list(map( lambda p: p[2], self.seq))
        for i in range(n-2,-1,-1):
            self.sums[i][0]+=self.sums[i+1][0]
            self.sums[i][1]+=self.sums[i+1][1]
            self.pairs[i]+=self.pairs[i+1]

        self.greedylist=[]
        self.fromlist=[]
        for i in range(n):
            self.fromlist.append(self.seq[i:])
            self.greedylist.append([[],[]])
            for j in [0,1]:
                self.greedylist[i][j]=list(map( lambda p: (p[j],p[2]), self.fromlist[i]))

                #self.greedylist[i][j].sort(cmp=lambda x,y: x[0]-y[0])

                self.greedylist[i][j] = sorted(self.greedylist[i][j],
                    key=cmp_to_key(lambda x,y: x[0]-y[0]))

    def checksmp(self,alfa,beta):
        

        try:
            return self.t[alfa][beta]
        except IndexError:
            self.t.extend( [[]]*(1+alfa-len(self.t)))
            self.t[alfa].extend([None]*(1+beta-len(self.t[alfa])))
                
        (al,be,z)=self.adjalfabeta(alfa,beta)
        if len(self.seq)<=1:
            if not len(self.seq): return z
            return self.minz(al, be, self.seq[0])+z
        return None

    def minz(self,al,be,d):
        res=d[2]
        if d[0]>0: res=min(res,al//d[0])
        if d[1]>0: res=min(res,be//d[1])
        return res

    def greedy(self,frompos,lim,ind):
        res=0
        for i in self.greedylist[frompos][ind]:
            v=i[1]
            if i[0]:
                v=min(v,lim//i[0])
                lim-=v*i[0]
            res+=v
        return res


    def seqpairm(self, alfa, beta, seq, inserted,coall=0):

        def prn(pos,al,be,cur,curmax):
            print ("al=%d be=%d cur=%d curmax=%d :::" % (al,be,cur,curmax),end='')
            for i in range(len(pos)):
                print (self.seq[i],pos[i],";",end='')
            print()

        self.callind+=1
        pos = list(map( lambda p: 0, self.seq))
        lastpos=0
        al=alfa
        be=beta
        curmax=0
        cur=0
        n = len(self.seq)        
        for i in range(n):
            pos[i] = self.minz(al, be, self.seq[i])
            al-=pos[i]*self.seq[i][0]
            be-=pos[i]*self.seq[i][1]
            curmax+=pos[i]
            cur+=pos[i]
            if pos[i]>0: lastpos=i

        if debug & DSPX: prn(pos,al,be,cur,curmax)

        while 1:
            #decrement lastpos
            ok=0
            
            while lastpos>=0:                
                if pos[lastpos]>0:
                    pos[lastpos]-=1
                    al+=self.seq[lastpos][0]
                    be+=self.seq[lastpos][1]
                    cur-=1                    
                    ok=1
                    if lastpos<n-1:
                        frompos=lastpos+1
                        if al>=self.sums[frompos][0] and be>=self.sums[frompos][1]:
                            sol=self.pairs[frompos]
                            if debug & DSPX2: print ("Smp solution al=%d be=%d fp=%d sol=%d" % (al,be,frompos,sol))
                            curmax=max(curmax,cur+sol)
                            ok=0
                        elif al>=self.sums[frompos][0]:
                            sol=self.greedy(frompos,be,1)
                            if debug & DSPX2: print ("Greedy solution al=%d be=%d fp=%d sol=%d" % (al,be,frompos,sol))
                            curmax=max(curmax,cur+sol)
                            ok=0
                        elif be>=self.sums[frompos][1]:
                            sol=self.greedy(frompos,al,0)
                            if debug & DSPX2: print ("Greedy solution  al=%d be=%d fp=%d sol=%d" % (al,be,frompos,sol))
                            curmax=max(curmax,cur+sol)
                            ok=0
                else:
                    lastpos-=1
                if ok: break

            if not ok:
                
                break

            prev=lastpos
            
            if debug & DSPX2: print ("Inc",al,be)
            # try to increment next positions
            ok=0
            while lastpos<n-1:
                lastpos+=1
                d=self.seq[lastpos]
                pos[lastpos]=self.minz(al, be, d)
                if pos[lastpos]>0:
                    al-=pos[lastpos]*d[0]
                    be-=pos[lastpos]*d[1]
                    
                    cur+=pos[lastpos]
                    ok=1

                    
            if ok:
                curmax=max(curmax,cur)
            else:
                # find lastpos

                if debug & DSPX2:
                    print ("not found, recovering lastpos",end='')
                    print ("Rec",al,be)
                while prev>=0 and not pos[prev]:
                    prev-=1
                if prev==-1:
                    
                    break
                lastpos=prev
                
            if debug & DSPX: prn(pos,al,be,cur,curmax)
        
        return curmax




################################################################
# ML ESTIMATION

def initializeInternalDP(gt,s):
    def finddual(gt,s,locked):
        """ returns list of tuples (rootg1,rootg2) of non-locked speciations nodes """
        res=[]
        for g in gt.postint:
            if g.spec==1 and g.lcamapping==s and ((locked and g.lock==1) or ((not locked) and g.lock!=1)):
                if g.l.lcamapping.num <=s.l.num: res.append((g.l,g.r))
                else: res.append((g.r,g.l))
        return res

    def findnondual(gt,s):
        """ returns tuple of two lists of non dual nodes ( left-roots,right-roots) """
        res1=[]
        res2=[]
        for g in gt.prefall:
            if s.l.minnum <= g.lcamapping.num <= s.l.num:
                # parent is mapped either outside S(s)
                # or on s and is a duplication

                if not g.parent:
                # g is root of tg
                # add g to parent of g.lcamapping
                    if g.lcamapping.parent.isdescendantof(s):  res1.append(g)
                elif g.parent.lcamapping.num > s.num or g.parent.lcamapping.num < s.minnum or \
                    (g.parent.lcamapping==s and g.parent.spec==0):
                    res1.append(g)

            if s.r.minnum <= g.lcamapping.num <= s.r.num:
                # parent is mapped either outside S(s)
                # or on s and is a duplication

                if not g.parent:
                # g is root of tg
                # add g to parent of g.lcamapping
                    if g.lcamapping.parent.isdescendantof(s):  res2.append(g)
                elif g.parent.lcamapping.num > s.num or g.parent.lcamapping.num < s.minnum or \
                (g.parent.lcamapping==s and g.parent.spec==0):
                    res2.append(g)

        return (res1,res2)

    s.dual=finddual(gt,s,0)
    s.lockeddual=finddual(gt,s,1)
    s.nondual=findnondual(gt,s)

def mlrInternalDP(gt,s,kappas,drml):
    """ DP for internal nodes """

    # MAIN DP starts here
    # compute parameters for seqpair

    initializeInternalDP(gt,s)

    seq=list(map( lambda p : ( p[0].ints,p[1].ints ), s.dual ))
    lseq=list(map( lambda p : ( p[0].ints,p[1].ints ), s.lockeddual ))

    gamma0=[]
    sumdual=[]
    s.sumlockeddual=[]

    for j in [0,1]:
        gamma0.append( reduce(lambda x,y : x+y, [0]+list(map( lambda t : t.ints , s.nondual[j]))))
        sumdual.append( reduce(lambda x,y : x+y, [0] + [ p[j] for p in seq ]) )
        s.sumlockeddual.append( reduce(lambda x,y : x+y, [0] + [ p[j] for p in lseq ]) )


    if debug & DDP:
        print ('SeqPair:',end='')
        print ('gamma0=',gamma0, ", sumdual=",sumdual, ", seq=",seq, "sumlockeddual=",s.sumlockeddual,end='')
        print()

    seqpair=SeqPairStructExt(seq)

    # for the reconstruction of duplication/speciation setting
    s.MLreconstr=dict()



    for kappa in kappas:
        s.MLreconstr[kappa]=[]
        if debug & DDP:
            print ( 'INTERNAL NODE:', s.labs(), 'kappa=%d (%d)' % (kappa,max(kappas))) #, "l-dup=%d l-spec=%d",
            print ('de si k0 k1 acc p q k0-p ka1-q -> si* M:(l,r,s)')
            print ('SeqPair:',end='')
            print ('gamma0=',gamma0, ", sumdual=",sumdual, ", seq=",seq, "sumlockeddual=",s.sumlockeddual,end='')
            print()

        curml=-10e1000
        curdef=0
        if DPTYPE==1:
            for delta in range(s.duplocked,kappa-len(s.l.locked)-len(s.r.locked)-s.speclocked+1):
                for sigma in range(s.speclocked,s.spec+1):
                    for kappa0 in range(len(s.l.locked),kappa-sigma-delta+1):
                        relkappa0=kappa0-s.l.minkappa
                        kappa1=kappa-sigma-delta-kappa0
                        relkappa1=kappa1-s.r.minkappa

                        if relkappa1>=len(s.r.ML) or relkappa1<0: continue
                        if relkappa0>=len(s.l.ML) or relkappa0<0: continue

                        if s.l.ML[relkappa0]==None: continue
                        if s.r.ML[relkappa1]==None: continue

                        if debug & DDP:
                            print  ("%2d %2d %2d %2d " % (delta,sigma,kappa0,kappa1),end='')

                        k0=kappa0-s.sumlockeddual[0]
                        k1=kappa1-s.sumlockeddual[1]

                        #check acceptability of (sigma,kappa0,kappa1)
                        ok=0
                        for p in range(max(0,k0-sumdual[0]),min(gamma0[0]+1,k0+1)):
                            
                            if ok==1: break
                            for q in range(max(0,k1-sumdual[1]),min(gamma0[1]+1,k1+1)):
                                sigmastar=seqpair.get(k0-p,k1-q)
                                # print "CHECK", k0-p,k1-q,seq,sigmastar
                                if sigma<=sigmastar+s.speclocked:
                                    ok=1
                                    if debug & DDP:
                                        print ("acc ",p,q,k0-p,k1-q,'->',sigmastar,end='')
                                    break
                        if not ok:
                            if debug & DDP: print ("not acc")
                            continue


                        cur=0.0+s.l.ML[relkappa0]+s.r.ML[relkappa1]+logP(s.v,delta)
                        if debug & DDP:
                            print (" : %.2f %.2f %.2f" % ( s.l.ML[relkappa0],s.r.ML[relkappa1],cur),end='')

                        if (not curdef) or (curml<cur):
                            curdef=1
                            curml=cur
                            if debug & DDP:
                                print (' NEW! M=%f' % curml)
                            s.MLreconstr[kappa]=[ (sigma,kappa0,kappa1) ]
                        else:
                            if curml==cur:                                
                                if debug & DDP: print( " *")
                                s.MLreconstr[kappa].append( (sigma,kappa0,kappa1) )
                            else:
                                if debug & DDP: print()
                                

        if curdef:
            
            s.ML[kappa-s.minkappa]=curml
            if verbose>=VERBOSEGEN:
                drml.log('%d/%d. INT: kappa=%d/%d logml=%f\n' % (s.num,drml.st.root.num,kappa,max(kappas),s.ML[kappa-s.minkappa]))
            if debug & DDP: print()


def mlrDP(gt,st,drml,firstrun=1,lazy=0):
    """ Main DP algorithm """

    def closure(g):
        """ Compute closure of locked nodes in g """
        # some nodes of gt are locked
        # g.lock - true if locked
        # g.inspec - true/false - speciation/duplication node
        if g.isleaf(): return
        if g.lock==1:
            for c in [ g.l, g.r ]:
                if c.isleaf() or c.lock==1 or c.spec: continue
                if g.spec: # speciation closure
                    if c.lcamapping.parent==g.lcamapping:
                        c.lock = 1
                        closure(c)
                else: # dupl. closure
                    if c.lcamapping == g.lcamapping:
                        c.lock = 1
                        closure(c)


    
    
        

    computelogfact(len(gt.postall)+1)

    if firstrun:
        lcamapping(gt,st)
        st.setpostfixnum()
        # for checking if sth needs to be recomputed
        for s in st.postall: s.diff = 1
        for g in gt.postall: s.oldlock = 0


    setlcadupspec(gt)

    for s in st.postall:
        s.dup=0  # number of lca-duplications
        s.phi=0  # number of spec+dupl nodes below s in the lca-embedding
        s.spec=0 # is speciation
        s.speclocked=0
        s.duplocked=0
        s.locked=set()

    # unlocked speciations (-1)
    for g in gt.postint:
        if g.lock==-1:
            g.spec=0  # new duplication
        else:
            g.spec=g.lcaspec   # lca spec remains

    # propagate locked nodes (closure of L)
    for g in gt.postint:
        if g.lock==1: closure(g)

    # initialize locked nodes for s
    for g in gt.postint:
        g.lcamapping.phi+=1
        if g.lock==1:
            g.lcamapping.locked|=set(g._internalpostfix)

    # Propagate locked nodes
    for s in st.postint:
        s.locked|=s.l.locked|s.r.locked

    for g in gt.postint:
        if g.spec:
            g.lcamapping.spec+=1
            if g.lock==1: g.lcamapping.speclocked+=1
        else:
            g.lcamapping.dup+=1
            if g.lock==1: g.lcamapping.duplocked+=1

    if not firstrun:
        for s in st.postall: s.diff = 0
        for g in gt.postall:
            if g.lock!=g.oldlock or ((not g.isleaf()) and g.oldspec!=g.spec): g.lcamapping.diff=1
        # propagate diff values
        for s in st.postint:
            s.diff=s.diff or s.l.diff or s.r.diff

        #print g.labs(),g.lock, "SPEC" if g.spec else "DUP", g.lcaspec

    for g in gt.postall:
        g.ints=0   # number of internal nodes
        if g.isleaf(): g.ints=0
        else: g.ints=g.l.ints+g.r.ints+1

    for s in st.postint:
        s.phi=s.phi+s.l.phi+s.r.phi

    if debug & DDP:
        print ('  # dup spc phi sp-lck dp-lck s-node')
        for s in st.postall:
            print ('%3d %3d %3d %3d %4d   %4d   %s' % (s.num,s.dup,s.spec,s.phi,s.speclocked,s.duplocked,s.labs()),s.locked)

    # main dynamic programming loop
    for s in st.postall:
        if verbose>=VERBOSEGEN:
            drml.log('\nDYNPROG START: ' + s.labs()  + "\n")
        
        if lazy and s.diff:

            s.minkappa=len(s.locked)
            s.ML=[None]*(s.phi-s.minkappa+1)
            
            if s.isleaf():
                # DP for leaves
                for delta in range(s.minkappa,s.phi+1):
                    p=logP(float(s.v),delta)
                    if verbose>=VERBOSEGEN:
                        drml.log('%d/%d. LEAF: kappa=%d logml=%f\n' % (s.num,st.root.num,delta,p))
                    s.ML[delta-s.minkappa]=p
            else:
                # DP for internal nodes
                if s==st.root: rng=[s.phi]
                else: rng=list(range(s.minkappa,s.phi+1))
                mlrInternalDP(gt,s,rng,drml)

    for g in gt.postall:
        g.oldlock=g.lock
        if not g.isleaf(): g.oldspec=g.spec

    return st.root.ML[st.root.phi-s.minkappa]


# END OF ML ESTIMATION
################################################################

recpaircnt=0

def genreconcilable(gt,st,drml):
    """  Generate evolutionary scenarios from ds-settings """

    def grinit(s):
        """ Initialize node s for the enumeration """

        s.cur=list(range(s.up))
        s.curdup=s.grcand[:]

        if not s.isleaf():
            if not grinit(s.l): return 0
            if not grinit(s.r): return 0
            s.curdup.extend(s.l.grup)
            s.curdup.extend(s.r.grup)

        s.grup=[]

        while len(s.cur) > len(s.curdup):

            if debug & DDS:
                print ("WARNING unsufficient number of candidates!")

            if not grnext(s.l):
                grinit(s.l)
                if not grnext(s.r):
                    if debug & DDS: print ("No reconciliations")
                    return 0
            s.curdup=s.grcand[:]
            s.curdup.extend(s.l.grup)
            s.curdup.extend(s.r.grup)

        for i in s.cur: s.grup.append(s.curdup[i])

        return 1

        # the first is always correct (see TPOSTFIX in the generation)


    def grnext(s):

        """ Generate next variant of for the node s """

        s.grup=[]
        while 1:

# TODO: IMPROVE THIS PART BY FILTERING COMBINATIONS

            c=nextcombination(s.cur,s.up,len(s.curdup))

            if c:
                l=[]
                for i in s.cur: l.append(s.curdup[i])
                ok=1
                for g in l:
                    if not g.parent: continue # root
                    if g.parent in l: continue
                    if not g.parent in s.curdup:
                        if not (g.parent.spec and g.parent.speclock and (g.parent.lcamapping==s.parent)):
                            continue
    
    
                    ok=0
                    break

                if ok:
                    s.cur=c
                    s.grup=l
                    return 1

            else:
                break

        # no new combination
        # gen next for bottom nodes

        s.curdup=s.grcand[:]
        s.cur=list(range(s.up))
        res=1

        if s.isleaf(): res=0
        else:
            if not grnext(s.l):
                grinit(s.l)
                if not grnext(s.r):
                    res=0
            s.curdup=s.grcand[:]
            s.curdup.extend(s.l.grup)
            s.curdup.extend(s.r.grup)

        return res

    # lock nodes below locked speciation
    def lock(g,s):
        if g.lcamapping!=s: return
        g.speclock=1
        if g.isleaf(): return
        lock(g.l,s)
        lock(g.r,s)

    def grlab(s):
        return s.labs()

    def getdupnum(g,s):
        if g.isleaf() and s.isleaf(): return 0

        if g.grmap!=s and g.grmap.isdescendantof(s.r): return getdupnum(g,s.r)
        if g.grmap!=s and g.grmap.isdescendantof(s.l): return getdupnum(g,s.l)

        if not g.grspec and g.grmap==s: return getdupnum(g.l,s)+getdupnum(g.r,s)+1
        if g.grspec and g.grmap==s and g.l.grmap.isdescendantof(s.l) and g.r.grmap.isdescendantof(s.r):
            return getdupnum(g.l,s.l)+getdupnum(g.r,s.r)  # spec
        if g.grspec and g.grmap==s and g.r.grmap.isdescendantof(s.l) and g.l.grmap.isdescendantof(s.r):
            return getdupnum(g.r,s.l)+getdupnum(g.l,s.r) # spec

        return 0
        

    def getdls(g,s):
        """ print DLS tree - see the formula in the article """
        if g.isleaf() and s.isleaf(): return g.label

        if g.grmap!=s and g.grmap.isdescendantof(s.r):
            return '('+getdls(g,s.r)+','+grlab(s.l)+' - ) ~'  # loss/spec

        if g.grmap!=s and g.grmap.isdescendantof(s.l):
            return '('+getdls(g,s.l)+','+grlab(s.r)+' - ) ~'  # loss/spec


        if not g.grspec and g.grmap==s:
            return '('+getdls(g.l,s)+','+getdls(g.r,s)+' ) +'  # dupl
        if g.grspec and g.grmap==s and g.l.grmap.isdescendantof(s.l) and g.r.grmap.isdescendantof(s.r):
            return '('+getdls(g.l,s.l)+','+getdls(g.r,s.r)+' ) ~'  # spec
        if g.grspec and g.grmap==s and g.r.grmap.isdescendantof(s.l) and g.l.grmap.isdescendantof(s.r):
            return '('+getdls(g.r,s.l)+','+getdls(g.l,s.r)+' ) ~'  # spec

        return "error? UNKNOWN"

    ############################################################
    # genreconcilable starts here

    stnodes=st.postall
    recpaircnt=0
    for s in stnodes:
        s.grcand=[]
        s.up=s.phi-s.kappa

    for g in gt.postint:
        g.grspec=0 # all-duplications

    for g in gt.prefall:
        g.speclock=0
        g.grmap=g.lcamapping

    for s in st.postint:
        for g in s.specsetting:
            sproot=s.dual[g][0].parent
            sproot.grspec=1
            sproot.speclock=1
        for p in s.lockeddual:
            sproot=p[0].parent
            sproot.grspec=1
            sproot.speclock=1

    for s in st.postint:
        for g in s.specsetting:
            lock(s.dual[g][0],s.l)
            lock(s.dual[g][1],s.r)
        for p in s.lockeddual:
            lock(p[0],s.l)
            lock(p[1],s.r)

    for g in gt.nodes(TPREFIX,TINT):
        if not g.speclock:
            g.lcamapping.grcand.append(g)

    if not grinit(st.root): return 0

    if debug & DDS:
        for s in stnodes:
            print ("%s u%d k%d d%d s%d: " % ( s.labs(), s.up, s.kappa, s.delta, s.sigma ),end='')
            for g in s.grcand: print (g.printsmp(),end='')
            print()

    while 1:
        recpaircnt+=1
        if debug & DDS:
             print ("Reconcilable pair: %d.%d" % (dupspecsettingcnt, recpaircnt))
             for s in stnodes:
                 print ("%s:" % s.labs(),s.cur," --> ",end='')
                 for g in s.curdup: print  (g.printsmp(),end='')
                 print()

        if verbose>=VERBOSEDLS or drml:
            dls=getdls(gt.root,gt.root.grmap)
            dupnum=getdupnum(gt.root,gt.root.grmap)
            if verbose>=VERBOSEDLS:
                print ("DLSTree", drml.dpstats.logml, getdls(gt.root,gt.root.grmap))
                print()
            if drml: drml.newscenario(dls,dupnum)


        if recpaircnt>0 and drml.single: break # report single reconciliation

        if not grnext(st.root): break

    return recpaircnt

reconcilablecnt=0
specnodessettingcnt=0

def genfinalrec(gt,st, drml):
    """ Check speciation nodes and generated scenario """

    global dupspecsettingcnt
    global specnodessettingcnt

    specnodessettingcnt+=1
    if debug & DDS:
        print ("Spec nodes setting: %d.%d " % (dupspecsettingcnt,specnodessettingcnt))

        
    rej=0

    for s in st.postint:
        if debug & DDS:
            print  (s.labs(), ': dual=', s.dual, 'ss=', s.specsetting , ' si/spec=%d/%d' % (s.sigma, s.spec))

        alpha=s.sumlockeddual[:]
        for j in [0,1]:
            for ss in s.specsetting:
                alpha[j]+=s.dual[ss][j].ints

        if debug & DDS:
            print( " dual->spec: ", alpha)
            print( " kappas(s,lchild,rchild): ",s.kappa,s.l.kappa,s.r.kappa)
            print (" s(phi,delta):",s.phi,s.delta)

        if s.l.kappa < alpha[0]:
            rej=1
            if debug & DDS: print ("REJECTED (l-kappa) - too many locked in spec")
            break

        if s.r.kappa < alpha[1]:
            rej=1
            if debug & DDS: print ("REJECTED (r-kappa) - too many locked in spec")
            break

        if s.phi - alpha[0]-alpha[1]-s.sigma < s.delta:
            rej=1
            if debug & DDS: 
                print ("REJECTED (delta) - too few for duplications")
            break

        # alpha[0,1] - sum of number of internal nodes of in dual trees locked by speciations

    if not rej:
        if debug & DDS: print ("#%d accepted" % reconcilablecnt)
        return genreconcilable(gt,st,drml)


    return 0


################################################################################
# Speciation nodes setting generator

def gendsrec(gt,st,drml):
    """ Generate speciation nodes in gt from the speciation setting """

    def initspecsetting(s):
        """ Initialize setting of speciation nodes"""
        return list(range(s.sigma-s.speclocked))


    def nextspecsetting(s,c):
        """ Generate next speciation nodes setting """
        if s.sigma:
            c=nextcombination(c, s.sigma-s.speclocked, len(s.dual))
            if c: return (c, 0)
        return ([], 1)


    # Main procedure starts here

    cnt=0
    internal=st.internal()
    #internal.sort(cmp=lambda x,y: x.size-y.size)
    internal = sorted(internal, key=cmp_to_key(lambda x,y: x.size-y.size))
    for s in internal:
        s.specsetting = initspecsetting(s)

    # Main loop for enumeration of speciation settings
    while 1:

        cnt+=genfinalrec(gt,st,drml)
        if cnt>0 and drml.single: break # report single reconciliation

        found=0
        for s in internal:
            (specsetting,eoss) = nextspecsetting(s, s.specsetting)
            if not eoss:
                found=1
                break
            else:
                s.specsetting = initspecsetting(s)  # new

        if not found: break

    return cnt




def gends(gt,st,drml):
    """ DS settings enumeration based on DP output """

    def initdssettings(s,kappa,pos):
        """ Initialize DS settings generator for node s """

        s.kappa=kappa
        if s.isleaf():
            s.delta=kappa  # duplications for a leaf
            s.sigma=0      # leaf has no speciation nodes
            return True

        #print kappa,pos
        #print s.MLreconstr

        curset=s.MLreconstr[kappa][pos]

        s.delta=kappa-curset[0]-curset[1]-curset[2]
        s.sigma=curset[0]
        initdssettings(s.l,curset[1],0)
        initdssettings(s.r,curset[2],0)

    """ Main procedure starts here """

    global dupspecsettingcnt
    totalcnt=0
    ds0=0

    st.setsize()
    for s in st.prefall: s.pos=0
    internal=st.internal()
    #internal.sort(cmp=lambda x,y: x.size-y.size)
    internal = sorted(internal, key=cmp_to_key(lambda x,y: x.size-y.size))
    initdssettings(st.root,st.root.phi,0)

    if verbose>=VERBOSEGEN:
        print ("\nDuplication-speciation (DS) settings generator\n")

    while 1:

        dupspecsettingcnt+=1
        drml.newdssetting()
        if verbose>=VERBOSEDS:
            print ("DS setting #%d" % dupspecsettingcnt)
            print (st.printdr())
            print ()

        if debug & DDS:
            print ("d s k p node:len (only for non-zero dskp's)")
            for s in st.postall:
                if s.delta or s.sigma or s.kappa or s.phi:
                    print ("%d %d %d %d %s:%.1f" % (s.delta,s.sigma,s.kappa,s.phi,s.labs(),s.v))

        if drml.solve>=SOLVEDLS:

            recnt=gendsrec(gt,st,drml)
            if not recnt:
                drml.newds0setting()
                ds0+=1
            totalcnt+=recnt
            if verbose>=VERBOSEGEN:
                print ("%d reconciliation(s) found for DS setting #%d"  %  ( recnt,dupspecsettingcnt))
                print ()
            if recnt>0 and drml.single: break  # report one reconciliation
        found=0

        for s in internal:
            s.pos+=1
            if s.pos<len(s.MLreconstr[s.kappa]):
                initdssettings(s,s.kappa,s.pos)
                found=1
                break
            else:
                s.pos=0
                initdssettings(s,s.kappa,s.pos)

        if not found: break
    ds=dupspecsettingcnt
    dupspecsettingcnt=0
    return (totalcnt,ds,ds0)


def printprobtable():
    """ Print logP values """
    print ()
    print ("TABLE of PROB. VALUES logP(branchlen,duplications)")
    rdupl = list(range(0,10))
    rbl = list(range(1,10))
    print ('           dup=',end='')
    for dupl in rdupl: print ("%-4d " %dupl,end='')
    print ()

    for branchlen in rbl:
        print ("branchlen=%d  " % branchlen,end='')
        for dupl in rdupl:
            print ("%4.2f " %logP(branchlen, dupl),end='')
        print()


dpsteps=0
refactor=0

def mlr(gt,st,drml):
    """ Main procedure for computing estimations & reconciliations  """

    global verbose, verbosehi, debug, refactor, refactorfnd

    refactorfnd = 0

    def singlemle(gt,st,drml,hi=0):
        """ Single ml estimation with constraints
            drml - global drml
            hi - defined if hard instance
        """
        global verbose, verbosehi, debug
        oldverbose = verbose

        if hi:
            # hard instance processing
            verbose = verbosehi
            drml.bbsteps+=1
        ml=mlrDP(gt,st,drml,not hi,1)

        if hi and drml.stats.dls and ml < drml.stats.logml:
            # ml too small to have significant scenarios
            drml.log("(DLS generator skipped)")
        else:
            drml.initdpstats(ml)
            if drml.solve>=SOLVEDS:
                gends(gt,st,drml)
        verbose = oldverbose
        return (ml,drml.dpstats.dls)

    def searchforsetting(st,gt,specnodes,gtnodes,drml,depth):

        def testsinglenode(g,lval,drml,depth):
            global verbose
            drml.log(" (%d) lca-spec: " % depth, g.num, "raised" if lval==-1 else "locked", )
            for x in gtnodes: x.lock=0  # for all - important
            for x in specnodes: x.lock=x.curlock # set current locking positions
            g.lock=lval
            g.curlock=lval
            res=singlemle(gt,st,drml,1)+(g,lval)
            g.curlock=0
            return res

        global verbose, debug, refactorfnd

        bbstepno=drml.bbsteps
        drml.log("%d. Level %d start\n" % (bbstepno,depth))

        mspecnodes=[ g for g in specnodes if g.curlock==0 ]
        if not mspecnodes:
            drml.log("No speciation candidates (searchforsetting)\n")
            drml.log("Error?\n")
            sys.exit(-10)

        oldv=verbose
        ml=-10000e100
        res=[]

        scoreweight=1#/depth

        sn=mspecnodes

        if drml.bbalg & BBSCORES: 
            #sn.sort(cmp=lambda x,y: int(y.score-x.score))
            sn = sorted(sn, key=cmp_to_key(lambda x,y: int(y.score-x.score)))
        if drml.bbalg & BBRAND: sn= [ random.choice(sn) ]
        if drml.bbalg & BBRANDX:
            snx=[ i for i in sn if random.random()>0.5 ]
            if len(snx)>0: sn=[ snx[0] ]
            else: sn=[ sn[0] ]

        if drml.bbalg & BBPAIR: sn=[ sn[0] ]

        results=[]
        lastml=drml.stats.logml
        for g in sn:
            for i in [1,-1]:
                r=testsinglenode(g,i,drml, depth)
                results.append(r)
                if r[1]:
                    res.append(r)
                    drml.log( "Resolved %f\n" % r[0])
                else: drml.log("No DLS %f\n" % r[0])
            if len(res): break

        if len(res)==2:
            ml=max(res[0][0],res[1][0])
            drml.log( "%d. Level %d node %d: resolved %f\n" % (bbstepno,res[0][2].num,depth, ml))
            if drml.stats.dls and ml == drml.stats.logml:
                if lastml and lastml!=ml:
                    
                    
                    for x in specnodes: x.score=0.0

                res[0][2].score+=1.0*scoreweight


        elif len(res)==1:
            ml=res[0][0]
            g=res[0][2]
            str="%d. Level %d node %d:" % (bbstepno,depth,g.num)
            drml.log( "%s one case resolved %f\n" % ( str, ml))
            if drml.stats.dls and ml == drml.stats.logml:
                if lastml and lastml!=ml:
                    print ([ x.score for x in specnodes ])
                    print ("RESET")
                    for x in specnodes: x.score=0.0
                g.score+=0.5*scoreweight

            skip=0
            if res[0][3]==-1:  # locked was computed but with no scenario
                lockedml=results[len(results)-2][0]
                if drml.stats.dls and lockedml < drml.stats.logml:
                    # cut this case - no solutions
                    drml.log("%s raised case rejected (%f < %f)\n" % (str,lockedml,drml.stats.logml))
                    skip=1

            if not skip:
                drml.log( "%s continue with one constraint\n" % str)
                g=res[0][2]
                g.curlock=-res[0][3]   # reverse lock
                ml1=searchforsetting(st,gt,specnodes,gtnodes,drml,depth+1)
                if refactor and refactorfnd: return 0
                g.curlock=0
                ml=max(ml,ml1)

            drml.log( "%s resolved %f\n" % (str, ml))
        else:
            drml.log( "%d. Level %d: no solution\n" % (bbstepno,depth))
            drml.timelog()
            g=sn[0]
            str="%d. Level %d node %d:" % (bbstepno, depth,g.num)
            ml1=results[0][0]
            ml2=results[1][0]
            if drml.stats.dls and ml1 < drml.stats.logml:
                drml.log("%s locked case rejected (%f < %f)\n" % (str,ml1,drml.stats.logml))
            else:
                drml.log( "%s recursion with locked speciation\n" % str)
                g.curlock=1
                ml1=searchforsetting(st,gt,specnodes,gtnodes,drml,depth+1)

            if drml.stats.dls and ml2 < drml.stats.logml:
                drml.log("%s raised case rejected (%f < %f)\n" % (str,ml2,drml.stats.logml))
            else:
                drml.log( "%s recursion with raised speciation\n" % str)
                g.curlock=-1
                ml2=searchforsetting(st,gt,specnodes,gtnodes,drml,depth+1)


            g.curlock=0

            ml=max(ml1,ml2)

            drml.log( "%s resolved %f\n" % (str, ml))

        oldv=verbose
        return ml

    gt.setpostfixnum()
    drml.start(gt,st)

    gtnodes=gt.nodes(TPREFIX)  # faster termination

    for g in gtnodes: g.lock=0

    ml=singlemle(gt,st,drml,0)

    refactor=1

    if not drml.dpstats.dls and drml.solve==SOLVEHI:

        drml.log("\nHard instance detected\n")
        drml.log("Starting second phase - branch and bound algorithm\n")
        drml.log("BB algorithm: %d\n\n" % drml.bbalg)

        specnodes=[ g for g in gtnodes if g.lcaspec ]
        #specnodes.sort(cmp=lambda x,y: y.ints-x.ints)
        specnodes = sorted(specnode, key=cmp_to_key(lambda x,y: y.ints-x.ints))

        # initialize curlock
        for x in specnodes:
            x.curlock=0
            x.score=0.0

        # inform drml
        drml.hiprocessing()

        ml=searchforsetting(st,gt,specnodes,gtnodes,drml,1)

        if drml.himode==1:

            drml.stophiprocessing()

            drml.log("\nGenerating all optimal reconciliations\n")
            drml.log("\nFound %d configuration(s) of locked/raised lca-speciations\n" % len(drml.historage))

            drml.initdpstats(ml)
            drml.single=0
            drml.himode=0

            for (higt,hist) in drml.historage:
                gends(higt,hist,drml)

    else:

        if drml.dpstats.dls>0:
            drml.log("This is a valid MLE (reconciliations are found)\n")
        else:
            if drml.solve<=SOLVEDS:
                drml.log("The second phase of algorithm (-dh or -dha) is required to ensure correctnes of ML\n")

    drml.stop()

    return ml



# gen dls stats for a given set of dls trees and a species tree
def gendlsstats(tr,s,drml,proctype):
    global gsmatch
    oldgsmatch=gsmatch
    gsmatch=None

    s.root.setdls()
    for n in s.nodes(): 
        n.totaldlsstat=dict(zip(dlstypes,({},{},{},{})))
    gsmatch=oldgsmatch

    clusterdict=dict((frozenset(n.cluster),n) for n in s.nodes())

    mlfile="dls2ml.log"

    if 'm' in proctype:
        fm=open(mlfile,"w")

    for t in tr:
        
        computelogfact(len(t.postall)+1)
        
        for n in s.nodes(): 
            n.dlsstat=dict(zip(dlstypes,(0,0,0,0)))

        for n in t.nodes():
            sn=clusterdict[frozenset(n.cluster)]
            sn.dlsstat[n.dlstype]=sn.dlsstat[n.dlstype]+1

        if 'm' in proctype:
            curp=0.0
            for n in s.nodes(): 
                curp+=logP(float(n.v),n.dlsstat['+'])

            print (curp,t)
            fm.write("%.2f %f " %(curp,curp))
            for n in s.nodes(): 
                fm.write("%s=%d " % ("".join(n.cluster),n.dlsstat['+']))
            fm.write("%s\n" % t)
            
        for n in s.nodes():
            for tp in dlstypes:
                num=n.dlsstat[tp]
                n.totaldlsstat[tp][num]=n.totaldlsstat[tp].get(num,0)+1
                
    if 'm' in proctype:
        fm.close()
        drml.log("File %s created\n" % mlfile)   

    if 'f' in proctype:
        for n in s.nodes():
            print (" ".join(n.cluster))
            for tp in dlstypes:
                print (tp,n.totaldlsstat[tp])
        print()        

    if 'd' in proctype:        
        f=open("drml.txt","w")

        for t in tr:
            f.write("%s\n" %t)

        # write!
        drml.log("File drml.dls created\n")
        f.close()
    

# Print all mappings
def mapgen(gt,st):
    lcamapping(gt,st)
    for g in gt.prefall: g.m=g.lcamapping
    gint=gt.postint
    fnd=1
    cnt=1
    while fnd:
        print ("MAPPING",cnt)
        cnt+=1
        for g in gint:
            print (g.printsmp()," ->",g.m.labs())
        fnd=0
        for g in gint:
            if g.m!=st.root and ((g==gt.root) or (g!=gt.root and g.parent.m!=g.m)):
                g.m=g.m.parent
                fnd=1
                break
            g.m=g.l.m.lca(g.r.m)




def lcamlvalue(gt,st):   
    lcamapping(gt,st)
    setlcadupspec(gt)
    stnodes=st.postall
    for g in gt.nodes(): g.m=g.lcamapping
    for s in stnodes: s.mde=0
    tot=0
    for g in gt.postint:
        if not g.lcaspec: 
            g.m.mde+=1
            tot+=1

    curval=0
    for s in stnodes:
        curval+=logP(s.v,s.mde)

    return curval,tot

def fatmlvalue(gt,st):   
    curval=0
    stnodes=st.postall
    for s in stnodes:
        if s==st.root: 
            curval+=logP(s.v,len(gt.allleaves)-1)
        else:
            curval+=logP(s.v,0)
    return curval,len(gt.allleaves)-1


def exhaustiveml(gt,st,drml):
    """ Exhaustive ml - main function """

    def exhaustivemlst(stnodes):
        """  Compute log ml """
        curval=0
        stnodes=st.postall
        for s in stnodes:
            s.ml=[]
        
            for si in range(s.msi+1):
                s.ml.append(logP(s.v,s.mde+si))
            curval=curval+max(s.ml)
            s.mlsubtree=max(s.ml)
            if not s.isleaf():
                s.mlsubtree+=s.l.mlsubtree+s.r.mlsubtree

        return curval


    def exhaustivemlstopt(stnodes):
        curval=0
        for s in stnodes:
            if s.changed:
                ml=-100e1000
                for si in range(s.msi+1):
                    ml=max(ml,logP(s.v,s.mde+si))
                s.changed=0
                s.ml=ml

            curval+=s.ml

        return curval


    print()
    print ("Exhaustive ML")
    print()

    if debug & DEX:
        print( "-G \"",gt,"\"",end='')
        print("-S \"",st,"\"")

    lcamapping(gt,st)
    setlcadupspec(gt)
    stnodes=st.postall
    for g in gt.nodes(): g.m=g.lcamapping
    gint=gt.postint
    fnd=1
    cnt=0
    maxml=-100e100
    mlcnt=0

    if verbose>=VERBOSEGEN: print ("Mappings generator")
    for s in stnodes:
        s.msi=0
        s.mde=0
        s.changed=1
    for g in gint:
        if g.lcaspec:
            g.m.msi+=1
            g.mspec=1
        else:
            g.m.mde+=1
            g.mspec=0


    while fnd:
        cnt+=1
        if verbose>=VERBOSEGEN: print ("%3d." % cnt,end='')
        if verbose>=VERBOSEGEN:
            print ("Mapping")
            for g in gint: print (' ',g.printsmp()," ->",g.m.labs())

        if verbose==VERBOSEGEN: curml=exhaustivemlstopt(stnodes)
        else: curml=exhaustivemlst(stnodes)
        if curml==maxml:
            mlcnt+=1
        if curml>maxml:
            maxml=curml
            mlcnt=1
            if verbose<VERBOSEGEN:
                print ("%d. logML: " % cnt,maxml)

        if debug & DEX:
            print ("\n s de si ml[lessdup..alldup] P(edge) mlsubtree")
            for s in stnodes:
                print (' ',s.labs(),s.mde,s.msi,s.ml,max(s.ml),s.mlsubtree)

        if verbose==VERBOSEGEN: print ("logML=%.4f %.4f" %(curml,maxml))

        fnd=0
        for g in gint:
            curmap=g.m
            if curmap!=st.root and ((g==gt.root) or (g!=gt.root and g.parent.m!=curmap)):
                dmap=g.m.parent
                if g.mspec:   # new dupl
                    curmap.msi-=1
                    g.mspec=0
                else:
                    curmap.mde-=1
                curmap.changed=1
                g.m=dmap
                if g!=gt.root:
                    gpar=g.parent
                    if gpar.mspec and dmap==gpar.m:  # parent into dupl
                        gpar.mspec=0
                        gpar.m.msi-=1
                        gpar.m.mde+=1
                        gpar.m.changed=1
                dmap.mde+=1
                dmap.changed=1
                fnd=1
                break

            dmap=g.l.m.lca(g.r.m)

            if dmap==curmap: continue  # nop
            g.m=dmap
            if not g.mspec:
                curmap.mde-=1
                if g.lcaspec and dmap==g.lcamapping and g.l.m!=dmap and g.r.m!=dmap:
                    # dupl->spec
                    g.mspec=1
                    dmap.msi+=1
                    # check parent
                else:
                    # dupl-> lower dupl
                    dmap.mde+=1
                dmap.changed=1
                curmap.changed=1
                if g!=gt.root:
                    gpar=g.parent
                    if not gpar.mspec and gpar.lcaspec and gpar.m==gpar.lcamapping and gpar.l.m!=gpar.m and gpar.r.m!=gpar.m:
                        # par dup->spec
                        gpar.mspec=1
                        gpar.m.msi+=1
                        gpar.m.mde-=1
                        gpar.m.changed=1

    if verbose >= VERBOSEDLS:
        print ('Mappings processed: %d\n' % cnt)
    drml.log('# Exhaustive ML' )
    drml.log(gt.__repr__())
    drml.log('\n')
    drml.log(st.__repr__())
    drml.log('\nMappings: %d\n' % cnt)
    drml.log('\nLogML: %f\n' % maxml)

    return maxml



def main(argv=None):
    """ Main routine with program arguments """

    global verbose, debug, verbosehi

    # initialize the list of species names - for random tree generators
    def initLVS(n):
        lvs=[]
        for i in range(n):
            m=i/26
            c=chr(ord('a')+i%26)
            if m>0: lvs.append(c+"%d" % m)
            else: lvs.append(c)
        return lvs

    if argv is None:
        argv = sys.argv

    gtsize=5
    stsize=5
    randlen=6
    lvs=initLVS(stsize)
    gt=0
    st=0
    readdlsfiles=0

    tt=[]

    drml=Drml("drml", SOLVEDLS, 1, BBPAIR+BBSCORES)

    try:
        opts, args = getopt.getopt(sys.argv[1:], "D:b:d:fM:mo:x:ht:r:S:G:P:p:s:g:Rl:Tenv:V:Fc:D:L:X:y")
    except getopt.GetoptError:
        
        usage()
        sys.exit(2)

    if not opts:
        usage()
        sys.exit(2)

    #drml.basename=None
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()

        elif opt == '-M':
            global gsmatch
            gsmatch=arg

        elif opt =='-o':
            drml.initbasename(arg)        

        elif opt == '-d':
            if arg=='m':
                drml.solve=SOLVEML
            elif arg=='s':
                drml.solve=SOLVEDS
            elif arg=='r':
                drml.solve=SOLVEDLS
                drml.single=1
            elif arg=='h':
                drml.solve=SOLVEHI
                drml.single=1
            elif arg=='ra':
                drml.solve=SOLVEDLS
                drml.single=0
            elif arg=='ha':
                drml.solve=SOLVEHI
                drml.single=0
            else:
                print( "Unknown opt in -d. Expected: m s r h ra ha.")
                sys.exit(-1)

        elif opt == '-X':
            debug=int(arg)

        elif opt == '-b':
            drml.bbalg=int(arg)

        elif opt == '-D':
            readdlsfiles=arg

        elif opt == '-l':
            randlen=int(arg)

        elif opt == '-v':
            verbose=int(arg)
        elif opt == '-V':
            verbosehi=int(arg)
        elif opt == '-S':
            stsize=int(arg)
            lvs=initLVS(stsize)
        elif opt == '-G':
            
            gtsize=int(arg)

        elif opt == '-L':
            ##-dptype=2
            global lmbd
            lmbd=float(arg)

    for opt, arg in opts:
        if opt == '-g':
            gt=Tree(arg)
        elif opt =='-s':
            st=Tree(arg)
        elif opt =='-P' or opt =='-p':
            if opt=='-p' and drml.basename=="drml":
                # new basename
                if arg.rfind('.txt')!=-1:
                    drml.initbasename(arg[0:arg.rfind('.txt')])

            f=open(arg)
            lt=[]
            p=0
            for ln in f:
                
                if len(ln)==0: continue
                if ln[0]=='#':
                    
                    continue
                if ln[0] in '\r\n#%/':
                    continue  # skip comments, empty lines
                if p:
                    tt.append( (p,Tree(ln)) )
                    p=0
                else: p=Tree(ln)
            f.close()
        elif opt =='-r':
            for i in range(int(arg)):
                gtx=randomtree(lvs,0,gtsize)
                gtx.sortnodes()

                stx=randomtree(lvs[:],1)
                stx.sortnodes()
                stx.randombranchlenghts(1,randlen)
                tt.append((gtx,stx))


    if readdlsfiles:
        st=None # only the first will be processed
        stln=''
        tr=[]
        for fn in args:
            print ("Reading file: ",fn)
            f=open(fn)
            lt=[]
            p=0
            for ln in f:
                
                ln=ln.strip()
                if len(ln)==0: continue
                
                if ln[0]=='#':
                    
                    continue
                if ln[0] in '\r\n#%/':
                    continue  # skip comments, empty lines

                # each gene tree is ignored
                if p==1:
                    if not st: 
                        st=Tree(ln)  # first species tree 
                        stln=ln
                    else:
                        if ln!=stln:
                            print ("Error in file %s: found different species trees in input files" % fn)
                            print ("#%s#" % st)
                            print ("#%s#" % stln)
                            sys.exit(-1)

                # dls tree found
                if p>1: tr.append(DlsTree(ln))
                p+=1
            f.close()

        gendlsstats(tr,st,drml,readdlsfiles)

    if gt and st: tt.append( (gt,st) )
    gt=0
    st=0
    computelogfact(3000)

    for opt, arg in opts:

        # mappings only
        if opt == '-e':
            for (gt,st) in tt: mapgen(gt,st)

        # exhaustive ml
        elif opt == '-n':
            gt,st=tt[0]
            ml=exhaustiveml(gt,st,drml)
            print ("Exhaustive logml=",ml)

        # ml
        elif opt == '-m':
            cnt=0
            gt,st=tt[0]
            mlr(gt, st, drml)
            if len(tt)>1:
                print ("\nPlease run separately other pairs from the input file")

        elif opt == '-R':
            for gt, st in tt:
                print (gt)
                print (st)


        # compare exhaustive and dpml
        elif opt == '-c':
            for gt,st in tt:

                def check(res,ml):
                    diff=0
                    for i in res:
                        if abs(i-ml)>0.000001:
                            print (res,ml)
                            print ("Difference!",)
                            log.write("Warning! Difference")
                            if i<ml: log.write(" negative\n")
                            else: log.write(" positive\n")
                            diff=1
                            return 0
                    print ("Check: OK")
                    return res+[ml]

                alg=[]
                if arg=='a':
                    alg=[BBSCORES,BBSINGLETONS,BBPAIR,-1]
                elif arg=='s':
                    alg=[BBSCORES,BBSINGLETONS,BBPAIR]
                else:
                    alg=[bbalg,-1]

                res=[]
                for a in alg:
                    if len(res)>1: drml.log("Compare: ", res)
                    if a>0:
                        bbalg=a
                        res=check(res,mlr(gt,st,scendet))
                    else:
                        verbose=VERBOSEDLS
                        res=check(res,exhaustiveml(gt,st,drml))
                    if not res: sys.exit(-1)

                print ("Compare: ", res)
                log.write("Compare: ")
                log.write(res.__repr__())
                log.write(" - Equal\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())


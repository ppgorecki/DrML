# DrML - Maximum Likelihood Estimation in the Duplication Loss Model

DrML is a program for estimating maximum likelihood of evolutionary scenarios and computing optimal evolutionary scenarios in the duplication-loss model based on results presented in [1].

## Input/output

Drml generates at most two files:
drml.log - log file
drml.dls - output file that contains:
  gene tree (from the input)
  species tree (from the input)
  dlstree(s)

A DLS-tree is a tree, representing evolutionary scenario, written in the nested parenthesis notation with additional decoration: + denotes a duplication node, ~ denotes a speciation node, - denotes a loss node.

Both file names can be changed with -o or -p options. See below for more details.

## Options:

### Defining gene and species trees.

-s TREE - defines a rooted species tree with branch lengths
-g TREE - defines rooted gene tree

-P FILE - input file containing a pair of gene and species trees (# starts a comment):
    GENE TREE
    SPECIES TREE

-p FILE - like above, but the log file is FILE.log instead of drml.log

### Random trees

-r NUM - generate NUM random pairs of trees
-G NUM - number of leaves in random gene tree
-S NUM - number of leaves in random species tree
-l NUM - define max length of species tree branches (for random generator)

### General options

-v NUM - verbose level
    0 - summary only
    1 - summary and scenarios (DLS trees)
    2 - summary, DS-settings and DLS trees
    3 - summary, dynammic programming output, DS-settings and DLS trees

-V NUM - verbose level for hard instances processing

-o BASENAME
    Defines basename. Default is drml. If -p FILE[.txt] option is given then the basename is FILE.

-R print pairs of input trees (gene tree, species trees,...)

### Maximum likelihood estimation

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

### Exhaustive estimation and comparing results.

-n - ML by exhaustive enumeration, for testing on small trees only

### Gene-species tree labels matching.
   
-M type - how the labels of the leaves in the gene tree determine a species name 

This option should be used when the labels in gene trees are different from the species names.

default (no -M) - exact match 
pNUM[,LEN] - where NUM>0 - species names start at NUM position in each gene label, LEN (optional) is the length of species names
p-NUM[,LEN] - where NUM>0 - the position is determined from the end of gene labels
aDELIM[,LEN] - DELIM is a string, species names start after DELIM
bDELIM[,LEN] - species names are before DELIM

Example:

```
drml.py -s "(a,b)" -g "(a1,(b1,a2))" -m -Mp0,1
```

### Processing dls files.

-Df - print frequencies of events (duplications, losses and speciations)
-Dm - compute ml values for each DLS tree
-Dd - print dls trees only 

Note that the species tree should be the same for all dls files.

Example:

```
drml.py -Df *.dls
```

Or, if matching rule is needed:
```
python3  drml.py -Df  -Mp0,1 *.dls
```

### Exemplary datasets

The trees taken from version 7 of: http://treefam.org. Branch lengths in the species tree taken from http://www.timetree.org/. Rooted trees from unrooted TreeFam trees computed by UREC: https://bitbucket.org/pgor17/urec/.

In data dir:
tf7.txt - TreeFam 7 species tree with branch lengths.
tftrees.tgz - input files in the format "GeneTree EOLN SpeciesTree" (use with -p option) for each TreeFam 7 tree.

To execute on a single file, e.g., TF101051.txt:
```
drml.py -p TF101051.txt -m
```

### Bibliography

[1] Maximum likelihood models and algorithms for gene tree evolution with duplications and losses P. Górecki, G. Burleigh, O. Eulenstein, BMC bioinformatics 12 (Suppl 1), S15, 2011

[2] DrML: Probabilistic Modeling of Gene Duplications
Paweł Górecki and Oliver Eulenstein
Published Online:30 Dec 2013
https://doi.org/10.1089/cmb.2013.0078
 
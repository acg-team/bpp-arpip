---
layout: page
title: Examples
---

## 1: Inferring an ancestral sequences with extended substitution model from nucleotide input sequences

### Prepare the configuration file

We list all the parameters required by ARPIP in a text file named `conf.txt` 
(the order of the parameters is not relevant). For more information on the syntax of the parameter.

```
analysis_name = TEST 0
alphabet=DNA
input.sequence.file=../data/input/test_00/msa.fasta
input.sequence.sites_to_use=all
init.tree=user
input.tree.file=../data/input/test_00/tree.newick
model=PIP(model=JC69,lambda=10,mu=0.01)
opt.seed=1
opt.likelihood=0
opt.pip_param_estimate=0
rate_distribution=Constant
output.msa.file=../data/output/test_00/msa.fasta
output.tree.file=../data/output/test_00/tree.nwk
output.ancestral.file=../data/output/test_00/anc.fasta
output.node_rel.file=../data/output/test_00/node_rel.txt
output.mlindelpoints.file=../data/output/test_00/mlindelpoints.txt

```

### Execute the analysis
```
$ ARPIP params=../data/input/test_00/conf.txt
```
---
## 2: Inferring an ancestral sequences with extended substitution model from amino-acids input sequences

### Prepare the configuration file

We list all the parameters required by ARPIP in a text file named `conf.txt` 
(the order of the parameters is not relevant). For more information on the syntax of the parameter.

```
analysis_name = TEST 1
alphabet=Protein
input.sequence.file=../data/input/test_01/sim-0_msa_new.fasta
input.sequence.sites_to_use=all
init.tree=user
input.tree.file=../data/input/test_01/sim-0_new.newick
model=PIP(model=WAG01,lambda=10,mu=0.01)
opt.seed=1
opt.pip_param_estimate=0
output.msa.file=../data/output/test_01/msa.fasta
output.tree.file=../data/output/test_01/tree.nwk
output.ancestral.file=../data/output/test_01/anc.fasta
output.node_rel.file=../data/output/test_01/node_rel.txt
output.mlindelpoints.file=../data/output/test_01/mlindelpoints.txt

```

### Execute the analysis
```
$ ARPIP params=../data/input/test_01/conf.txt

```
---
## 3: Inferring an ancestral sequences with extended substitution model from amino-acids input sequences. The indel parameters are not provided and therefore should be inferred from the input sequences.  


The initial values of insertion and deletion rates of the PIP model are inferred from the aligned 
input sequences (MSA), using Brent method.  

### Prepare the configuration file

We list all the parameters required by ARPIP in a text file named `conf.txt` 
(the order of the parameters is not relevant). For more information on the syntax of the parameter.

```
analysis_name = TEST 2
alphabet=Protein
input.sequence.file=../data/input/test_02/sim-0_msa_new.fasta
input.sequence.sites_to_use=all
init.tree=user
input.tree.file=../data/input/test_02/sim-0_new.newick
model=PIP(model=WAG01)
opt.seed=1
opt.pip_param_estimate=1
output.msa.file=../data/output/test_02/msa.fasta
output.tree.file=../data/output/test_02/tree.nwk
output.ancestral.file=../data/output/test_02/anc.fasta
output.node_rel.file=../data/output/test_02/node_rel.txt
output.mlindelpoints.file=../data/output/test_02/mlindelpoints.txt
output.pipparams.file=../data/output/test_02/param.txt

```
Please note by using `opt.likelihood=1` after computing the loglikelihood  value the program would be terminated. 
This feature is for specific usage for likelihood computation.

### Execute the analysis
```
$ ARPIP params=../data/input/test_02/conf.txt
```
---
## 4: Inferring an ancestral sequences with extended substitution model from amino-acids input sequences. The tree is not provided and therefore inferred from the input sequences.  


In this case just the MSA and indel parameter is provided by user. 
These can be provided by the user when known. If the tree is not provided then the tool computes a 
distance matrix from the pairwise alignments and infers a rooted guide tree using the supported tree reconstructin method (wpgma, upgma, nj, bioinj) 
For the purpose of illustration in this example the Upgma algorithm was used.

### Prepare the configuration file

We list all the parameters required by ARPIP in a text file named `conf.txt` 
(the order of the parameters is not relevant). For more information on the syntax of the parameter.

```
analysis_name = TEST 3
alphabet=Protein
opt.seed=1
opt.likelihood=0
opt.pip_param_estimate=0
opt.combine_msa_asr=1
input.sequence.file=../data/input/test_03/sim0-msa_new.fasta
input.sequence.sites_to_use=all
init.tree=auto
init.tree=upgma
model=PIP(model=WAG01,lambda=10, mu=0.01)
rate_distribution=Constant
output.ancestral.file=../data/output/test_03/anc.fasta********
output.tree.file=../data/output/test_03/tree.nwk
output.node_rel.file=../data/output/test_03/node_rel.txt
output.mlindelpoints.file=../data/output/test_03/mlindelpoints.txt
```

By using `model=PIP(model=WAG01,lambda=10,mu=0.1)` we can use the known indel parameters.

### Execute the analysis
```
$ ARPIP params=../data/input/test_03/conf.txt
```

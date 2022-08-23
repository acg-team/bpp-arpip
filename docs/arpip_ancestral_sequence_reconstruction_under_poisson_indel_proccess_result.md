---
layout: page
title: Results
---
The output files could be interpreted as follows:

#### Output tree file:

The new rooted phylogenetic tree could be the same tree if the user provided it as input. 
Otherwise, the tree is the estimated from data. Please note that the output tree is rooted Bifurcated, no matter
what was the input tree. (since the current version of bio++ does not assign internal nodes name, one can find it in
 the relation file)

#### Output joint ASR file:

The ASR output file is in fasta format. Each row in the file represents the internal node named in postorder of the tree.
You can find the names in "node relation" file. The file contains the most likely 
ancestral value with indels (gap) according to the joint reconstruction in a manner similar to FastML.

#### Output Maximum Likelihood Indel Points:

The output is per site ordered based on the original MSA file. The format is as follow:

    node#:{insertion/Deletion)          node#:X meaning deletion and node#:I meaning insertion.
For example: V2:X;V1:X;root:I; means insertion at "root" and multiple deletions at V2 and V1.  
Please note that per site, we only have one insertion and zero or more deletion points.
   
#### Output tree ancestral relation file:

In the reconstructed phylogeny the relation between parent node and the children nodes are specified. 

#### Output PIP parameter file:

In the case that user used "opt.likelihood=1" meaning to infer the PIP parameter. This file contains "lambda", "mu" and
"Log likelihood" values per file.

    
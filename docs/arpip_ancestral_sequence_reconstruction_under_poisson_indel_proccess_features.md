---
layout: page
title: Features and project structure
---

#  Features and project structure 
---
# Input options
---
    
    
### Reading sequences

    input.sequence.file={path}                              The aligned sequence file to use. 

                                                            Will consider {n} random sites, with optional replacement.

*The following formats are currently supported:*

    Fasta(extended={bool}, strictNames={bool})              The fasta format. The argument extended, default to 'no' allows to enable the HUPO-PSI
                                                            extension of the format. The argument strict_names, default to 'no', specifies that
                                                            only the first word in the fasta header is used as a sequence names, the rest of the
                                                            header being considered as comments.
    Mase(siteSelection={chars})                             The Mase format (as read by Seaview and Phylo_win for instance), with an optional site
                                                            selection name.
    Phylip(order={interleaved|sequential}, type={classic|extended}, split={spaces|tab})
                                                            The Phylip format, with several variations. The argument order distinguishes between
                                                            sequential and interleaved format, while the option type distinguished between the
                                                            plain old Phylip format and the more recent extention allowing for sequence names
                                                            longer than 10 characters, as understood by PAML and PhyML. Finally, the split
                                                            argument specifies the type of character that separates the sequence name from the
                                                            sequence content. The conventional option is to use one (classic) or more (extended)
                                                            spaces, but tabs can also be used instead.
    Clustal(extraSpaces={int})                              The Clustal format.
                                                            In its basic set up, sequence names do not have space characters, and one space splits
                                                            the sequence content from its name. The parser can however be configured to allow
                                                            for spaces in the sequence names, providing a minimum number of space characters is
                                                            used to split the content from the name. Setting extraSpaces to 5 for instance, the
                                                            sequences are expected to be at least 6 spaces away for their names.
    Dcse()                                                  The DCSE alignment format. The secondary structure annotation will be ignored.
    Nexus()                                                 The Nexus alignment format. (Only very basic support is provided)
    GenBank()                                               The GenBank not aligned sequences format.
                                                            Very basic support: only retrieves the sequence content for now, all features are
                                                            ignored.

### Reading trees

    input.tree.file={path}                                  The phylogenetic tree file to use.
    input.tree.format={Newick|Nexus|NHX}                    The format of the input tree file.
---
#  Alphabet options
---


    alphabet={DNA|RNA|Protein)|Codon(letter={DNA|RNA},type={Standard|EchinodermMitochondrial|InvertebrateMitochondrial|VertebrateMitochondrial})}
                                                            The alphabet to use when reading sequences. DNA and RNA alphabet can in addition take
                                                            an argument: **bangAsgap={bool}**
                                                            Tell is exclamation mark should be considered as a gap character. The default
                                                            is to consider it as an unknown character such as 'N' or '?'.    
---
#  Options
---

    opt.seed={real}                                         Sets the seed value of the random number generator.
 
    opt.likelihood={0|1}                                    1: The user wants to know what is the value of joint likelihood of tree and MAS under PIP
                                                            0: Deactive this option. By default it is 0.
                                        
    opt.pip_param_estimate={0|1}                            1: The user does not know what are the evelutionary parameter (i.e. lambda and mu) and
                                                            wants program to compute them. 0: other way.
    opt.tree.scale={real}                                   Set the scale value to scale the branch lengths.
    
    opt.unknown_as_gap={0|1}                                1: This software does not support ambiguity in characters.
                                                            We are kindly ask users to remove the unknown chars, o.w. the 
                                                            software change them to gap and in the next step will remove 
                                                            all the only-gap columns from MSA file. By default this flag is 0.
    
    opt.combine_msa_asr={0|1}                               1: The user can see the result along with their corresponding MSA.
                                                            It is recommended that user to activate this flag when using 
                                                            'unknown_as_gap'. In the case of having column full of gap the length
                                                            of input MSA and ASR are not the same. By default this flag is 0.                                                        
    
    
---
#  Initial tree options
---

    init.tree={user|auto}                                   Set the method for the initial tree to use.
                                                            The user option allows you to use an existing file passed via input.tree.file
                                                            This file may have been built using another method like neighbor joining or
                                                            parsimony for instance. The random option picks a random tree, which is handy
                                                            to test convergence.  This may however slows down significantly the optimization
                                                            process.
    init.distance.method={wpgma|upgma|nj|bionj}             When distance method is required, the user can specify which algorithm to use.
                                                            With the option 'infere_distance_matrix' the algorithm calculates the distance 
                                                            matrix from the pairwise alignments and then with bionj algorithm infers the tree.
    

If the `init.tree metod=user`, then refer to the option you find in "Reading trees".


---
#  Evolutionary model options
---

### Substitution models


    model={string}                                          A description of the substitution model to use, using the keyval syntax.

*The following nucleotide models are currently available as a core model:*

    JC69
    K80([kappa={real>0}])
    F84([kappa={real>0}, theta={real]0,1[}, theta1={real]0,1[},theta2={real]0,1[} ,"equilibrium frequencies"])
    HKY85([kappa={real>0}, theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[} ,"equilibrium frequencies"])
    T92([kappa={real>0}, theta={real]0,1[} ,"equilibrium frequencies"])
    TN93([kappa1={real>0}, kappa2={real>0}, theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[} ,"equilibrium frequencies"])
    GTR([a={real>0}, b={real>0}, c={real>0}, d={real>0}, e={real>0}, theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[} ,"equilibrium frequencies"])
    L95([beta={real>0}, gamma={real>0}, delta={real>0}, theta={real]0,1[}, theta1={real]0,1[}, theta2={real]0,1[} ,"equilibrium frequencies"])
    SSR([beta={real>0}, gamma={real>0}, delta={real>0}, theta={real]0,1[}])
    RN95([thetaR={real]0,1[}, thetaC={real]0,1[}, thetaG={real]0,1[}, kappaP={real[0,1[}, gammaP={real[0,1[}, sigmaP={real>1}, alphaP={real>1}])
    RN95s([thetaA={real]0,0.5[}, gamma={real]0,0.5[}, alphaP={real>1}])

*The following protein models are currently available as a core model:*

    JC69
    DSO78
    JTT92
    WAG01
    LG08
    
*The following meta models are currently available:*

    PIP13(model={model description}, {lambda={real>0}, mu={real>0})

If you leave the 'lambda' and 'mu' empty then the program would estimate them using Brent's method. 
Please note that this algorithm is designed to work the 'PIP13' model. 


### Rate across site distribution

    rate_distribution={rate distribution description}       Specify the rate across sites distribution

*Only Constant rate is currently available:*

    Constant                                                Uses a constant rate across sites
    
---
#  Output options
---

### Output tree file

    output.tree.file={path}                                       The phylogenetic tree file to write to.
    output.tree.format={Newick|Nexus|NHX}                         The format of the output tree file.
    output.trees.file={path}                                      The file that will contain multiple trees.
    output.trees.format={Newick|Nexus|NHX}                        The format of the output tree file.

### Output alignment file

    output.msa.file={path}                                      Alignment used in the study.


### Output inferred file

    output.ancestral.file={path}                              Write ancestral seuqences inferred by algorithm.
    output.node_rel.file={path}                               Write the relation of nodes. It is important to idendifying the internal nodes.
    output.mlindelpoints.file={path}                          Write the inferred indel points.
    output.pipparams.file={path}                              Write the estimated PIP parameters if the user set opt.likelihood=1
# AlphaRING

AlphaRING is a customised implementation of [AlphaFold2](https://github.com/google-deepmind/alphafold) that uses the [RING4](https://ring.biocomputingup.it/) standalone package to capture non-covalent interactions at the atomic level in monomeric protein models. By piping together AlphaFold2 and RING4, AlphaRING extends upon both package's capabilities to predict the pathogenicity of any given missense variant based on the predicted changes in non-covalent bond formation between the wild-type monomeric protein model and its missense variant counterpart.


An associated AlphaRING paper is currently under review by the [RECOMB 2025](https://recomb.org/recomb2025/index.html) conference and will be made available upon completion. Additionally, the AlphaRING benchmark dataset and the code used to generate it will be made available soon. 

## Overview

![](https://github.com/loggy01/alpharing/blob/main/images/fig_1.png)
<p align="center">
  <b>Figure 1</b> Overview of the AlphaRING workflow
</p>


AlphaRING conducts the following workflow:

&nbsp;&nbsp;&nbsp;**1. Accept FASTAs:**
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;In this step, AlphaRING accepts two FASTAs: firstly, a FASTA representing a monomeric wild-type protein, and secondly, a FASTA representing a monomeric missense variant protein. We define a missense variant as differing from its paired wild-type by exactly one residue. AlphaRING will not accept anything else in place of this.
 

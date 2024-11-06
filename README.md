# AlphaRING

AlphaRING is a package for the prediction of pathogenicity of any given missense variant. The package is a customised implementation of [AlphaFold2](https://github.com/google-deepmind/alphafold) that models a monomeric wild-type protein and a variant counterpart, and captures their non-covalent bonds using [RING4](https://ring.biocomputingup.it/). AlphaRING uses the differences in non-covalent bond formation to predict the pathogenicity of the missense variant.

An AlphaRING manuscript is currently under review by [RECOMB 2025](https://recomb.org/recomb2025/index.html). AlphaRING benchmarking data will be made available soon.

## Overview

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="./images/fig_1.png">
  <source media="(prefers-color-scheme: light)" srcset="./images/fig_1.png">
  <img alt="Shows the AlphaRING pipeline, which is described below." src="./images/fig_1.png">
</picture>

<p align='center'> <strong>Figure 1</strong> Overview of the AlphaRING workflow </p>

For any given missense variant, AlphaRING conducts the following workflow:

1. **Accept FASTAs**: 

   In this step, AlphaRING accepts two FASTAs: a monomeric wild-type protein and a counterpart missense variant differing by a residue.

2. **Predict structures**: 

   In this step, AlphaFold2 is used to predict the structure of the wild-type and variant proteins. The best model of each is relaxed and extracted.

3. **Perform residue interaction network (RIN) analysis**

   In this step, RING4 is used to generate a RIN of both the wild-type and variant models, capturing their non-covalent interactions at the atomic level.

> [!NOTE]
> RING4 refers to bonds as "edges" and residues as "nodes". AlphaRING's source code extensively uses this terminology.

4. **Calculate residue weightings**

   In this step, AlphaRING assigns each non-covalent bond (hydrogen, ionic, π-cation, π-π stacking, and π-hydrogen) a weighting of importance to protein stability.

   Weightings are calculated using novel bond-type specific formulas that take into account bond-specific energy and geometry values provided by RING4. As the energy values provided by RING4 are 
   fixed for a given bond-type, our formulas consider the variable bond distance and angle to multiply energy by a value between 0 and 2. Both distance and angle can contribute a value between 0 
   and 1 towards the multiplier. More favourable distances and angles result in a larger multiplier. Therefore, a higher bond weighting indicates greater importance to protein stability. 

   Our bond-type specific formulas come in three flavours. The first flavour is used for instances where a shorter distance and smaller angle is favourable (π-cation and π-π stacking bonds):
   
   $$
   Bond_{weight} = energy \times \left( \left(1 - \left(\frac{distance}{distance_{max}}\right)\right) + \left(1 - \left(\frac{angle}{angle_{max}}\right)\right) \right)
   $$

   The second flavour of bond-type specific formulas is used for instances where a shorter distance and larger angle is favourable (hydrogen bonds):

   $$
   Bond_{weight} = energy \times \left( \left(1 - \left(\frac{distance}{distance_{max}}\right)\right) + \left(\frac{angle}{angle_{max}}\right)\right)
   $$

   The third and final flavour of bond-type specific formulas is used for instances where a shorter distance is favourable but the contribution of angle is negligible (ionic and π-hydrogen 
   bonds):

   $$
   Bond_{weight} = energy \times 2 \times \left(1 - \left(\frac{distance}{distance_{max}}\right)\right)
   $$

   Each residue's weight in the wildtype and variant proteins is then calculated by summing the weight of all its bonds.
   
7. **Calculate fold change (FC)**:

   In this step, AlphaRING calculates the fold change between the weight of the wild-type and variant residue at the position of residue substiution. Therefore, values futher from 1 indicate a 
   greater change in weighting.


9. **Calculate AlphaRING score**:

    In this step, the absolute log2 of the FC is taken to provide a final numeric value of pathogenicity, the AlphaRING score. This results in a minimum score of 0. Therefore, higher scores 
    indicate greater predicted pathogenicity.

## Installation

## Usage

## Downstream analysis





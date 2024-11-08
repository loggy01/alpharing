# AlphaRING

AlphaRING is a package designed to predict the pathogenicity of any given missense variant. It is a customised implementation of [AlphaFold2](https://github.com/google-deepmind/alphafold), which models a monomeric wild-type protein and a missense variant counterpart. Using [RING4](https://ring.biocomputingup.it/), AlphaRING captures their non-covalent bonds and analyses the differences in bond formation between the wild-type and variant proteins to predict the pathogenicity of the missense variant.

> [!NOTE]
> AlphaRING's source code refers to bonds as edges and residues as nodes.
>  
> An AlphaRING manuscript is currently under review by [RECOMB 2025](https://recomb.org/recomb2025/index.html). 
>
> AlphaRING benchmarking data will be made available soon...

## Overview

<picture>
  <source srcset="./images/fig_1.png">
  <img alt="Shows the AlphaRING workflow." src="./images/fig_1.png">
</picture>

<p align='center'> <strong>Figure 1</strong> Overview of the AlphaRING workflow </p>

For any given missense variant, AlphaRING conducts the following workflow:

1. **Accept FASTAs**: 

   In this step, AlphaRING accepts two FASTA files: one for the monomeric wild-type protein and one for the missense variant, which differs by a single residue.

2. **Predict structures**: 

   In this step, AlphaFold2 is used to predict the structure of the wild-type and variant proteins. The best model of each is relaxed and used going forward.

3. **Generate residue interaction networks (RINs)**

   In this step, RING4 is used to generate a RIN of both the wild-type and variant models, capturing their non-covalent bonds.

4. **Calculate residue weightings**

   In this step, AlphaRING assigns each non-covalent bond a weighting of importance to protein stability. Weightings are calculated using novel bond-type specific formulas that take into account bond-specific energy and geometry as calculated by RING4. As the energy values provided by RING4 are fixed for a given bond-type, our formulas consider the bond's variables, distance and angle, to multiply energy by a value between `0–2`. Both distance and angle 
   can contribute a value between `0–1` towards the multiplier. More favourable distances and angles result in a larger multiplier. Therefore, a higher bond weighting indicates greater importance to protein stability.

   Our bond-type specific formulas come in three flavours. The first flavour is used when a shorter distance and smaller angle is favourable (π-cation and π-π stacking):

>```math
>Bond_{weight} = energy \times \left( \left(1 - \left(\frac{distance}{distance_{max}}\right)\right) + \left(1 - \left(\frac{angle}{angle_{max}}\right)\right) \right)
>```

   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The second flavour is used when a shorter distance and larger angle is favourable (hydrogen):

>```math
> Bond_{weight} = energy \times \left( \left(1 - \left(\frac{distance}{distance_{max}}\right)\right) + \left(\frac{angle}{angle_{max}}\right)\right)
> ```

   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The third flavour is used when a shorter distance is favourable and angle is negligible (ionic and π-hydrogen):

>```math
> Bond_{weight} = energy \times 2 \times \left(1 - \left(\frac{distance}{distance_{max}}\right)\right)
> ```

   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Each residue's weight in the wild-type and variant proteins is calculated by summing the weight of its bonds.
   
5. **Calculate fold change (FC)**:

   In this step, AlphaRING calculates the FC between the weight of the wild-type and variant residue at the position of substitution. Therefore, values further from `1` indicate a 
   greater change in weighting.

6. **Calculate AlphaRING score**:

    In this step, the absolute log2 of the FC is taken to provide a final metric, the AlphaRING score. This results in a minimum AlphaRING score of `0`. Therefore, higher 
    AlphaRING scores indicate greater predicted pathogenicity.

## Dependencies

Before installation, you will need a machine running Linux and a modern NVIDIA GPU. In addition, ensure you have the following dependencies:  

- [AlphaFold2 databases](https://github.com/google-deepmind/alphafold/tree/f251de6613cb478207c732bf9627b1e853c99c2f#installation-and-running-your-first-prediction) (databases for AlphaFold-Multimer are not required) 

- [RING4](https://biocomputingup.it/services/download/) (version v4.0-2-ge939f57)  

- [Git](https://git-scm.com/downloads)

- [wget](https://www.tecmint.com/install-wget-in-linux/)
 
- [Miniconda](https://docs.anaconda.com/miniconda/)

## Installation

Firstly, `git clone` the `alpharing` repository into the same parent directory as RING4 and `cd` into it:

```bash
git clone --recurse-submodules https://github.com/loggy01/alpharing
cd ./alpharing
```

`wget` `stereo_chemical_props.txt` into the `alphafold` subdir:

```bash
wget -P alphafold/alphafold/common/ https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt
```

`cp` your RING4 directory into the `alpharing` directory:

```bash
cp -r ../<dir> ring
```
> [!NOTE]
> Replace `<dir>` with the name of your RING4 directory.

`conda create` the `alpharing` environment:

```bash
conda create -n alpharing -c bioconda -c conda-forge hhsuite hmmer kalign2 openmm=8.0.0 pdbfixer python=3.10
```

Finally, `conda activate` the `alpharing` environment and `pip install` necessary packages:

```bash
conda activate alpharing
pip install absl-py==1.0.0 biopython==1.79 chex==0.1.86 dm-haiku==0.0.12 dm-tree==0.1.8 immutabledict==2.0.0 jax==0.4.25 ml-collections==0.1.0 numpy==1.24.3 pandas==2.0.3 plotly==5.15.0 scipy==1.11.1 tensorflow-cpu==2.16.1 jaxlib==0.4.25+cuda11.cudnn86 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
```

## Usage

`conda activate` the `alpharing` environment and run the script `run_alpharing.py` in a shell script as follows:

```bash
conda activate alpharing
<path to run_alpharing.py> \
  --fasta_paths=<path to wild-type FASTA>,<path to variant FASTA> \
  --max_template_date=<yyyy-mm-dd e.g. 2020-01-15> \
  --data_dir=<path to AlphaFold2 databases dir> \
  --output_dir=<path to dir to save all results> \
  --uniref90_database_path=<path to uniref90.fasta> \
  --mgnify_database_path=<path to mgy_clusters_2022_05.fa> \
  --template_mmcif_dir=<path to mmcif_files dir> \
  --obsolete_pdbs_path=<path to obsolete.dat> \
  --bfd_database_path=<path to bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt> \
  --uniref30_database_path=<path to UniRef30_2021_03> \
  --pdb70_database_path=<path to pdb70> \
  --use_gpu_relax=<whether to relax with GPU i.e. True or False>
```
> [!NOTE]
> Replace argument values with your own.

> [!WARNING]
> These are the minimum requirements to run AlphaRING. Specific users (particularly HPC users) may need to include other modules and more in their shell script.  

## Downstream analysis

For any given missense variant, AlphaRING will generate two subdirs within `output_dir`: one for and named after the wild-type FASTA, and one for and named after the variant FASTA.

In addition to the regular AlphaFold2 and RING4 output, each subdir will contain a version of the `ringEdges` and `ringNodes` files with a `Weight` column attached. The latter file type is plotted and saved as an interactive residue `weightPlot` HTML:

<picture>
  <source srcset="./images/fig_2.png">
  <img alt="Shows a pair of wild-type and variant residue weightPlots." src="./images/fig_2.png">
</picture>

<p align='center'> <strong>Figure 2</strong> Example of a pair of wild-type and variant weightPlots</p>

The HTML can be used to visually compare weights of a substituted residue between a wild-type and variant, or be used for exploratory analysis of unsubstituted residues. Additionally, the HTML file can be saved as a PNG.

Finally, the variant subdir will contain an `alpharing_score.txt` file, which only contains the AlphaRING score of the missense variant.

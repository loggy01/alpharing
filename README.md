# AlphaRING

AlphaRING is a customised implementation of [AlphaFold](https://github.com/google-deepmind/alphafold) and [RING](https://ring.biocomputingup.it/) designed to predict the pathogenicity of any missense variant. 

To predict the pathogenicity of a missense variant, AlphaRING performs the following steps: 

1. Predicts the wild-type and variant protein structures.
2. Captures the wild-type and variant non-covalent bonds.
3. Uses changes in non-covalent bond formation by the variable residue to predict pathogenicity.

Therefore, AlphaRING is both scalable and transparent, enabling researchers to obtain explicit, biophysically-based predictions of the pathogenicity of missense variants.

> [!NOTE]
> Version 2 (AlphaRING-X) coming soon...here is a sneak peek from our ISMB/ECCB 2025 presentation: https://iscb.junolive.co/ISMB/live/exhibitor/ISMBECCB2025_poster_1714

## Overview

Here is an overview of the three-step AlphaRING workflow:

1. **Predict the wild-type and variant protein structures**:

   Firstly, AlphaRING accepts two protein FASTA files representing a wild-type:variant pair and predicts the structures of both using AlphaFold. The FASTA files must differ by exactly one residue.

   Here is an example of a pair of FASTA files representing proteins called "wildtype" and "variant", which differ only at residue three:

   ```
   >wildtype
   MKTAY
   ```
   
   ```
   >variant
   MKMAY
   ```

2. **Capture the wild-type and variant non-covalent bonds**:

   Secondly, AlphaRING passes the predicted structures of the wild-type and variant proteins to RING to capture the energies, distances, and geometries of the non-covalent bonds within both structures.

3. **Use changes in non-covalent bond formation by the variable residue to predict pathogenicity**:

   Thirdly, AlphaRING uses the changes in non-covalent bond formation by the variable residue (along with the predicted accuracy of the wild-type structure at the position of the variable residue) in novel formulas to calculate a singular numerical metric: AlphaRING score. A greater change in non-covalent bond formation corresponds to a higher AlphaRING score, i.e., greater missense variant pathogenicity.

For more information on the AlphaRING workflow, please refer to the associated [publication](https://www.biorxiv.org/content/10.1101/2024.11.12.623182v2), which will be updated soon.

## Installation

Before installation, ensure you have a Linux machine equipped with a modern NVIDIA GPU, CUDA 11.x, and cuDNN 8.6.x, and then do the following:

- [Download full AlphaFold genetic databases](https://github.com/google-deepmind/alphafold?tab=readme-ov-file#genetic-databases) 
- [Install RING (v4.0-2-ge939f57)](https://biocomputingup.it/services/download/)
- [Install Miniconda](https://docs.anaconda.com/miniconda/)

### Users

To use AlphaRING, please do the following:

1. Gather the necessary dependencies by creating an AlphaRING environment using Miniconda:

   ```bash
   conda create -n alpharing -c bioconda -c conda-forge python==3.10 hmmer kalign2 pdbfixer hhsuite==3.3.0 openmm==8.0.0
   ```

2. Activate the AlphaRING environment and install AlphaRING:

   ```bash
   conda activate alpharing
   pip install alpharing
   pip install jaxlib==0.4.25+cuda11.cudnn86 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
   ```

3. Copy your installation of RING into the AlphaRING environment (replace `/path/to/your/` with your actual paths): 

   ```bash
   cp -r /path/to/your/ring-4.0/* /path/to/your/miniconda3/envs/alpharing/lib/python3.10/site-packages/ring/.
   ```

You are now ready to use AlphaRING 🎉

### Developers

To modify and develop AlphaRING, please 
[add your SSH key to your GitHub account](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account) and do the following:

<details>
<summary><b>Instructions</b></summary>

<br>

1. Clone the AlphaRING GitHub repository:

   ```bash
   git clone --recurse-submodules git@github.com:loggy01/alpharing.git
   cd alpharing
   ```

2. Create the AlphaRING environment as described for users:

   ```bash
   conda create -n alpharing -c bioconda -c conda-forge python==3.10 hmmer kalign2 pdbfixer hhsuite==3.3.0 openmm==8.0.0
   ```

3. Activate the AlphaRING environment and install AlphaRING and the AlphaFold submodule (only needs to be done once):

   ```bash
   conda activate alpharing
   pip install -e .
   pip install jaxlib==0.4.25+cuda11.cudnn86 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
   pip install -e alphafold --no-deps
   ```

4. Begin development. Any changes made in the activated environment will be recognised and should be tested:

   ```bash
   pip install pytest
   pytest tests/
   ```

5. Before pushing to the remote or submitting a pull request, ensure that you install and test AlphaRING:

   ```bash
   pip install .
   pip install jaxlib==0.4.25+cuda11.cudnn86 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
   pytest tests/
   ```

</details>

## Usage

To predict the pathogenicity of a missense variant, activate the AlphaRING environment and execute the `run_alpharing` command as follows:

```bash
conda activate alpharing
run_alpharing \
  --fasta_paths=/path/to/wildtype/fasta,/path/to/variant/fasta \
  --output_dir=/path/to/output/ \
  --data_dir=/path/to/alphafold/database/ \
  --max_template_date=yyyy-mm-dd \
  --use_gpu_relax=True or False \
```

> [!WARNING]
> For correct functioning, ensure the wild-type FASTA file is listed before the variant in `--fasta_paths`. For more information on running AlphaRING, execute `run_alpharing --helpfull` in the activated environment.

## Downstream

AlphaRING stores the output for the wild-type and variant proteins in separate subdirectories within the directory specified by `--output_dir`. Each subdirectory is named after the basename of the FASTA file from which it was produced. Both subdirectories contain the default AlphaFold and RING outputs. Additionally, the variant subdirectory contains the file `alpharing_score.txt`, which stores the AlphaRING score of the variant.

Here is the structure of a wild-type or variant output subdirectory produced from a FASTA file named "wildtype.fa" or "variant.fa" respectively:

```
{wildtype|variant}/              
   alpharing_score.txt                                # only in variant
   confidence_model_{1,2,3,4,5}_pred_0.json           # AlphaFold
   features.pkl                                       # AlphaFold
   msas/                                              # AlphaFold
      bfd_uniref_hits.a3m                             # AlphaFold
      mgnify_hits.sto                                 # AlphaFold
      pdb_hits.hhr                                    # AlphaFold
      uniref90_hits.sto                               # AlphaFold
   ranked_{0,1,2,3,4}.cif                             # AlphaFold
   ranked_{0,1,2,3,4}.pdb                             # AlphaFold
   ranking_debug.json                                 # AlphaFold
   relaxed_model_{1|2|3|4|5}_pred_0.cif               # AlphaFold
   relaxed_model_{1|2|3|4|5}_pred_0.pdb               # AlphaFold
   relaxed_model_{1|2|3|4|5}_pred_0.pdb_ringEdges     # RING (edges = bonds)
   relaxed_model_{1|2|3|4|5}_pred_0.pdb_ringNodes     # RING (nodes = residues)
   relax_metrics.json                                 # AlphaFold
   result_model_{1,2,3,4,5}_pred_0.pkl                # AlphaFold
   timings.json                                       # AlphaFold
   unrelaxed_model_{1,2,3,4,5}_pred_0.cif             # AlphaFold
   unrelaxed_model_{1,2,3,4,5}_pred_0.pdb             # AlphaFold
```
> [!NOTE]
> To save time and avoid redundancy, when running a prediction, AlphaRING will check both FASTA files for an existing corresponding output subdirectory containing a relaxed model PDB file, and will skip the AlphaFold stage only for the FASTA files that do. This feature is especially useful for consecutive AlphaRING predictions when the wild-type FASTA file remains unchanged but the variant FASTA file varies.

## Citation

If you have used any aspect of the AlphaRING package, please cite the associated [publication](https://www.biorxiv.org/content/10.1101/2024.11.12.623182v2). In addition, please cite the AlphaFold [publication](https://www.nature.com/articles/s41586-021-03819-2) and RING [publication](https://academic.oup.com/nar/article/52/W1/W306/7660079).

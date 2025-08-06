# AlphaRING v2 (AlphaRING-X)

AlphaRING is a package designed for interpretable, protein structure-based prediction of missense variant deleteriousness.

To predict the deleteriousness of a missense variant, AlphaRING performs the following steps:

1. Predicts the structure of the wild-type protein using [AlphaFold](https://github.com/google-deepmind/alphafold) to extract the pLDDT of the substituted wild-type residue.
2. Converts wild-type structure into a residue interaction network using [RING](https://ring.biocomputingup.it/) to extract the degree of the substituted wild-type residue.
3. Uses wild-type structure to calculate the ΔΔG of the substitution using [FoldX](https://foldxsuite.crg.eu/) and extracts it.
4. Calculates the relative substitution position (RSP) along the protein.
5. Feeds pLDDT, degree, ΔΔG, and RSP into an in-house [XGBoost](https://github.com/dmlc/xgboost) classifier trained to classify missense variant deleteriousness.
6. Outputs the probability of deleteriousness and feature SHAP values to explain the prediction mechanistically.

## Installation

Before installation, ensure you have a Linux machine equipped with a modern NVIDIA GPU and the following:

1. [Full AlphaFold v2 genetic database](https://github.com/google-deepmind/alphafold?tab=readme-ov-file#genetic-databases)
2. [RING v4](https://biocomputingup.it/services/download/)
3. [FoldX v5.1](https://foldxsuite.crg.eu/academic-license-info)
4. [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main)

To install AlphaRING, please do the following:

1. Create an environment for AlphaRING using Miniconda:

   ```bash
   conda create -n alpharing -c bioconda -c conda-forge python==3.10 hmmer kalign2 pdbfixer hhsuite==3.3.0 openmm==8.0.0
   ```

2. Activate the environment and install AlphaRING:

   ```bash
   conda activate alpharing
   pip install alpharing
   ```

## Usage

To predict the deleteriousness of a missense variant, activate the AlphaRING environment and execute the `alpharing` command as follows:

```bash
alpharing \
  --fasta_path=... \
  --substitutions=... \
  --output_dir=... \
  --data_dir=... \
  --ring_exe_path=... \
  --foldx_exe_path=...
```

Argument breakdown:

- `--fasta_path`: path to a FASTA file representing the wild-type protein. See [here](https://github.com/loggy01/alpharing/tree/main/tests/test_data/input/protein.fa) for an example.
- `--substitutions`: list of one or more single-residue substitutions to be individually applied to the wild-type protein. Please represent the substitutions in FoldX format, e.g., WA70Y, where W is the wild-type residue, A is the chain, 70 is the substitution position, and Y is the variant residue. For AlphaRING, the chain should always be A. If multiple substitutions are provided, please separate them with commas only, e.g., WA70Y,WA80F.
- `--output_dir`: path to the directory that will store the output. 
- `--data_dir`: path to the directory of the full AlphaFold database.
- `--ring_exe_path`: path to the RING executable. This should remain in the original installation.
- `--foldx_exe_path`: path to the FoldX executable. This should remain in the original installation.

## Downstream

AlphaRING stores the output in a subdirectory within the directory specified by `--output_dir`. The subdirectory is named after the basename of the FASTA file specified by `--fasta_path`. In addition to the default AlphaFold, RING, and FoldX outputs, the subdirectory contains the file `alpharing_scores.txt`, which summarises the feature values, deleteriousness probability, and feature SHAP values of each substitution specified by `--substitutions`.

> [!NOTE]
> For efficiency, when running a prediction, AlphaRING will check if the FASTA file specified by `--fasta_path` already has its corresponding output subdirectory with an AlphaFold relaxed model file, and, if found, will skip running AlphaFold.

#!/usr/bin/env python

import glob
import os
import pickle
import re
import subprocess
from pathlib import Path

import pandas as pd
import shap
from absl import app, flags, logging
from Bio.PDB import PDBParser

flags.DEFINE_string(
    'fasta_path',
    None,
    'Path to a FASTA file representing the wild-type protein.'
)
flags.DEFINE_list(
    'substitutions',
    None,
    'List of one or more single-residue substitutions to be individually applied '
    'to the wild-type protein. Please represent the substitutions in FoldX format, '
    'e.g., WA70Y, where W is the wild-type residue, A is the chain, 70 is the '
    'substitution position, and Y is the variant residue. For AlphaRING, the '
    'chain should always be A. If multiple substitutions are provided, please'
    'separate them with commas only, e.g., WA70Y,WA80F.'
)
flags.DEFINE_string(
    'output_dir',
    None,
    'Path to the directory that will store the output.'
)
flags.DEFINE_string(
    'data_dir',
    None,
    'Path to the directory of the full AlphaFold database.'
)
flags.DEFINE_string(
    'ring_exe_path',
    None,
    'Path to the RING executable. This should remain in the original installation.'
)
flags.DEFINE_string(
    'foldx_exe_path',
    None,
    'Path to the FoldX executable. This should remain in the original installation.'
)

FLAGS = flags.FLAGS


def run_alphafold(fasta_path: str, output_dir: str, data_dir: str) -> str:
    """
    Generates/retrieves an AlphaFold relaxed model for the wild-type 
    protein and returns the file path.

    Args:
        fasta_path: Path to a FASTA file representing the wild-type protein.
        output_dir: Path to the directory that will store the output.
        data_dir: Path to the directory of the full AlphaFold database

    Returns:
        Path to the relaxed model file for the wild-type protein.
    """
    # Look for AlphaFold relaxed model file for wild-type protein. Run if not found
    logging.info('Checking for existing AlphaFold relaxed model file to use')
    subdir_name = Path(fasta_path).stem
    model_path_pattern = os.path.join(output_dir, subdir_name, 'relaxed*.pdb')
    if glob.glob(model_path_pattern):
        logging.info(
            f'Existing AlphaFold relaxed model file for "{fasta_path}" found at ' 
            f'"{glob.glob(model_path_pattern)[0]}". Skipping AlphaFold'
        )
    else:
        logging.info(
            f'No existing AlphaFold relaxed model file for "{fasta_path}". Running ' 
            'AlphaFold'
        )
        logging.info('Building AlphaFold command')
        alphafold_command = (
            f'python3 {os.path.join(Path(__file__).parent.parent, "alphafold", "run_alphafold.py")} '
            f'--fasta_paths={fasta_path} '
            f'--output_dir={output_dir} '
            f'--data_dir={data_dir} '
            f'--uniref90_database_path={os.path.join(data_dir, "uniref90", "uniref90.fasta")} '
            f'--mgnify_database_path={os.path.join(data_dir, "mgnify", "mgy_clusters_2022_05.fa")} '
            f'--bfd_database_path={os.path.join(data_dir, "bfd", "bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt")} '
            f'--uniref30_database_path={os.path.join(data_dir, "uniref30", "UniRef30_2021_03")} '
            f'--pdb70_database_path={os.path.join(data_dir, "pdb70", "pdb70")} '
            f'--template_mmcif_dir={os.path.join(data_dir, "pdb_mmcif", "mmcif_files")} '
            f'--obsolete_pdbs_path={os.path.join(data_dir, "pdb_mmcif", "obsolete.dat")} '
            '--max_template_date=2050-01-01 '
            '--use_gpu_relax=False'
        )
        logging.info(f'Running AlphaFold with command: "{alphafold_command}"')
        subprocess.run(alphafold_command, shell=True, check=True)
        logging.info('AlphaFold finished running')
    
    # Retrieve the AlphaFold relaxed model file
    logging.info('Retrieving AlphaFold relaxed model file')
    model_path = glob.glob(model_path_pattern)[0]
    logging.info(f'AlphaFold relaxed model file retrieved at "{model_path}"')
    
    return model_path


def run_ring(model_path: str, ring_exe_path: str) -> str:
    """
    Generates a RING residues (nodes) file from the AlphaFold relaxed 
    model of the wild-type protein and returns the file path. 

    Args:
        model_path: Path to the relaxed model file for the wild-type protein.
        ring_exe_path: Path to the RING executable.

    Returns:
        Path to the residues file for the wild-type protein.
    """
    # Run RING
    logging.info('Building RING command')
    ring_command = (
        f'{ring_exe_path} '
        f'-i {model_path} '
        f'--out_dir {Path(model_path).parent} '
        '--no_add_H '
        '--all_edges '
        '--relaxed'
    )
    logging.info(f'Running RING with command: "{ring_command}"')
    subprocess.run(ring_command, shell=True, check=True)
    logging.info('RING finished running')

    # Retrieve th RING residues file
    logging.info('Retrieving RING residues file')
    residues_pattern_path = f'{model_path}_ringNodes'
    residues_path = glob.glob(residues_pattern_path)[0]
    logging.info(f'RING residues file retrieved at "{residues_path}"')

    return residues_path


def run_foldx(model_path: str, substitutions: list, foldx_exe_path: str) -> str:
    """
    Generates a FoldX ΔΔG (scanning) file from the AlphaFold relaxed model 
    of the wild-type protein and given single-residue substitutions and 
    returns the file path.

    Args:
        model_path: Path to the relaxed model file for the wild-type protein.
        substitutions: List of one or more single-residue substitutions to be 
                       applied to the wild-type protein 
        foldx_exe_path: Path to the FoldX executable.
    
    Returns:
        Path to the ΔΔG file for the wild-type protein with the given substitutions.
    """
    # Run FoldX
    logging.info('Building FoldX command')
    model_path_parent = Path(model_path).parent
    model_path_name = Path(model_path).name
    model_path_stem = Path(model_path).stem
    substitutions = ','.join(substitutions)
    foldx_command = (
        f'{foldx_exe_path} '
        '--command=PositionScan '
        f'--pdb={model_path_name} '
        f'--pdb-dir={model_path_parent} '
        f'--positions={substitutions}'
    )
    logging.info(f'Running FoldX with command: "{foldx_command}"')
    subprocess.run(foldx_command, shell=True, check=True, cwd=model_path_parent)
    logging.info('FoldX finished running')

    # Retrieve the FoldX ΔΔGs file
    logging.info('Retrieving FoldX ΔΔGs file')
    ΔΔGs_pattern_path = os.path.join(
        model_path_parent, 
        f'PS_{model_path_stem}_scanning_output.txt'
    )
    ΔΔGs_path = glob.glob(ΔΔGs_pattern_path)[0]
    logging.info(f'FoldX ΔΔGs file retrieved at "{ΔΔGs_path}"')

    return ΔΔGs_path


def calculate_alpharing_score(
    model_path: str, 
    residues_path: str, 
    ΔΔGs_path: str,
    substitutions: list
) -> None:
    """
    Calculates the AlphaRING score and SHAP values for each substitution 
    and saves the results.

    Args:
        model_path: Path to the relaxed model file for the wild-type protein.
        residues_path: Path to the RING residues file for the wild-type protein.
        ΔΔGs_path: Path to the FoldX ΔΔG file for the wild-type protein 
                   with the given substitutions.
        substitutions: List of one or more single-residue substitutions to be 
                       applied to the wild-type protein.
    Returns:
        None (results are saved to alpharing_scores.txt in the output directory)
    """
    # For each substitution, extract the pLDDT feature of the wild-type residue
    logging.info('Extracting pLDDT of each substituted wild-type residue')
    substitution_positions = [int(re.search(r'\d+', s).group()) for s in substitutions]
    model_df = PDBParser().get_structure('', model_path)[0]['A']
    plddts = [
        ({r.get_id()[1]: r.child_list[0].get_bfactor() for r in model_df})[p] 
        for p in substitution_positions
    ]

    # For each substitution, extract the degree feature of the wild-type residue
    logging.info('Extracting degree of each substituted wild-type residue')
    residues_df = pd.read_csv(residues_path, sep='\t')
    degrees = pd.to_numeric(
        residues_df.loc[[p - 1 for p in substitution_positions], 'Degree']
    ).tolist()

    # For each substitution, extract the ΔΔG feature
    logging.info('Extracting ΔΔG of each substitution')
    ΔΔGs_df = pd.read_csv(ΔΔGs_path, sep='\t', header=None)
    ΔΔGs = pd.to_numeric(ΔΔGs_df.iloc[:, 1]).iloc[1::2].tolist()
    
    # For each substitution, calculate the RSP feature
    logging.info('Calculating RSP of each substitution')
    aa_sequence_length = len(residues_df)
    rsps = [p / aa_sequence_length for p in substitution_positions]

    # For each substitution, calculate the probability of being deleterious
    logging.info('Calculating AlphaRING score for each substitution')
    features = pd.DataFrame(
        {'pLDDT': plddts, 'Degree': degrees, 'ΔΔG': ΔΔGs, 'RSP': rsps},
        index=substitutions
    )
    classifier_pkg_path = os.path.join(
        Path(__file__).parent, 
        'alpharing_classifier_pkg.pkl'
    )
    with open(classifier_pkg_path, 'rb') as p:
        classifier_pkg = pickle.load(p)
    classifier = classifier_pkg['classifier']
    background = classifier_pkg['background']
    probabilities = classifier.predict_proba(features)[:, 1]

    # For each subsitution, calculate each feature's SHAP value
    logging.info('Calculating SHAP values for each substitution')
    shap_explainer = shap.Explainer(
        lambda features: classifier.predict_proba(features)[:, 1],
        background
    )
    shap_values = shap_explainer(features, silent=True)

    # Create a DataFrame with the results and save it
    logging.info('Saving results to file')
    results_df = pd.DataFrame({
        'Substitution': substitutions,
        'pLDDT': plddts,
        'Degree': degrees,
        'ΔΔG': ΔΔGs,
        'RSP': rsps,
        'label': [
            'Neutral' if p <= 0.2270
            else 'Deleterious'
            if p >= 0.2740 else 'Ambiguous'
            for p in probabilities
        ],
        'Probability': probabilities,
        'pLDDT SHAP': shap_values.values[:, 0],
        'Degree SHAP': shap_values.values[:, 1],
        'ΔΔG SHAP': shap_values.values[:, 2],
        'RSP SHAP': shap_values.values[:, 3],
    })
    results_df.to_csv(
        os.path.join(Path(model_path).parent, 'alpharing_scores.txt'),
        sep='\t',
        index=False
    )
    logging.info(f'AlphaRING has finished running')


def main(_):
    model_path = run_alphafold(FLAGS.fasta_path, FLAGS.output_dir, FLAGS.data_dir)
    residues_path = run_ring(model_path, FLAGS.ring_exe_path)
    ΔΔGs_path = run_foldx(model_path, FLAGS.substitutions, FLAGS.foldx_exe_path)
    calculate_alpharing_score(model_path, residues_path, ΔΔGs_path, FLAGS.substitutions)


def entry_point():
    flags.mark_flags_as_required([
        'fasta_path',
        'substitutions',
        'output_dir',
        'data_dir',
        'ring_exe_path',
        'foldx_exe_path'
    ])
    app.run(main)

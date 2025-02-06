#!/usr/bin/env python

import glob
import math
import os
import subprocess
from pathlib import Path

import pandas as pd
from absl import app, flags, logging
from Bio import SeqIO
from Bio.PDB import PDBParser

flags.DEFINE_list(
    'fasta_paths',
    None,
    'Paths to two FASTA files, the first being the wild-type and the second the ' 
    'variant. The paths should be separated by a comma. Both FASTA paths must ' 
    'have a unique basename as the basename is used for naming output.'
)
flags.DEFINE_string(
    'output_dir',
    None,
    'Path to the directory that will store all the results.'
)
flags.DEFINE_string(
    'data_dir',
    None,
    'Path to the directory of the AlphaFold database.'
)
flags.DEFINE_string(
    'max_template_date',
    None,
    'Maximum template release date for AlphaFold to consider. In the format of ' 
    'yyyy-mm-dd e.g. 2050-01-01. Important if using historical test sets.'
)
flags.DEFINE_boolean(
    'use_gpu_relax',
    None,
    'Whether to relax AlphaFold predictions on GPU. Relax on GPU can be much ' 
    'faster than CPU, so it is recommended to enable if possible. GPUs must be ' 
    'available if this setting is enabled.'
)

FLAGS = flags.FLAGS

BOND_WEIGHT_FORMULAS = {
    'HBOND': lambda bond: 
        bond['Energy'] * ((1 - (bond['Distance'] / 5.3)) + (bond['Angle'] / 180)),
    'IONIC': lambda bond: 
        bond['Energy'] * 2 * (1 - (bond['Distance'] / 4.5)),
    'PICATION': lambda bond: 
        bond['Energy'] * ((1 - (bond['Distance'] / 6.7)) + (1 - (bond['Angle'] / 45))),
    'PIPISTACK': lambda bond: 
        bond['Energy'] * ((1 - (bond['Distance'] / 7.3)) + (1 - (bond['Angle'] / 90))),
    'PIHBOND': lambda bond: 
        bond['Energy'] * 2 * (1 - (bond['Distance'] / 5.0))
}


def find_variable_residue(fasta_paths: list[str]) -> int:
    """
    Finds and returns the position of the variable residue between the 
    wild-type and variant FASTAs in the command-line flag 'fasta_paths'.

    Args: 
        fasta_paths: List of FASTA paths (wild-type and variant respectively).

    Returns:
        Variable residue position between wild-type and variant.

    Raises:
        ValueError: Exactly two FASTA paths have not been provided.
        ValueError: One or both FASTAs do not contain exactly one record.
        ValueError: One or both FASTAs contain an empty ID.
        ValueError: One or both FASTAs contain an empty sequence.
        ValueError: The two FASTAs have the same ID.
        ValueError: One or both FASTA sequences contain invalid residue(s).
        ValueError: The two FASTA sequences are not the same length.
        ValueError: The two FASTA sequences do not vary by exactly one residue.
    """
    # Check if the number of provided FASTA paths is not two
    if len(fasta_paths) != 2:
        raise ValueError('Exactly two FASTA paths have not been provided')
    
    # Check if both fastas do not contain exactly one non-empty record each 
    fastas = []
    for fasta_path in fasta_paths:
        fasta = list(SeqIO.parse(fasta_path, 'fasta'))
        if len(fasta) != 1:
            raise ValueError('One or both FASTAs do not contain exactly one ID and sequence')
        if not fasta[0].id.strip():
            raise ValueError('One or both FASTAs contain an empty ID')
        if not fasta[0].seq.strip():
            raise ValueError('One or both FASTAs contain an empty sequence')
        fastas.append(fasta[0])

    # Check if the two FASTAs contain the same ID
    wildtype_id, variant_id = fastas[0].id, fastas[1].id
    if wildtype_id == variant_id:
        raise ValueError('The two FASTAs have the same ID')
    
    # Check if not all residues in both FASTA sequences are valid
    wildtype_sequence, variant_sequence = fastas[0].seq, fastas[1].seq
    valid_residues = set("ACDEFGHIKLMNPQRSTVWY")
    for sequence in [wildtype_sequence, variant_sequence]:
        if not all(residue in valid_residues for residue in sequence):
            raise ValueError('One or both FASTA sequences contain invalid residue(s)')
    
    # Check if the two FASTAs do not contain sequences of the same length
    if len(wildtype_sequence) != len(variant_sequence):
        raise ValueError('The two FASTA sequences are not the same length')
    
    # Check if the two FASTA sequences do not vary by exactly one residue
    sequence_differences = []
    for position, (a, b) in enumerate(zip(wildtype_sequence, variant_sequence)):
        if a != b:
            sequence_differences.append((position, a, b))
    if len(sequence_differences) != 1:
        raise ValueError('The two FASTA sequences do not vary by exactly one residue')
    variable_residue = sequence_differences[0][0] + 1

    return variable_residue


def run_alphafold(alphafold_input: dict[str, any]) -> list[str]:
    """
    Runs AlphaFold and returns the model file path for both the wild-type
    and variant predictions.

    Args:
        alphafold_input: Dictionary of input parameters for AlphaFold.
    
    Returns:
        List of model paths (wild-type then variant).
    """
    # Determine if either FASTA file already has a model
    fasta_paths_to_run = []
    model_path_patterns = {}
    for fasta_path in alphafold_input['fasta_paths']:
        subdir_name = Path(fasta_path).stem
        model_path_pattern = os.path.join(alphafold_input['output_dir'], subdir_name, 'relaxed*.pdb')
        model_path_patterns[fasta_path] = model_path_pattern
        if glob.glob(model_path_pattern):
            logging.info(f'Found and will use pre-computed model for "{fasta_path}"')
        else:
            fasta_paths_to_run.append(fasta_path)

    # Run AlphaFold for FASTA files that do not have a model
    if fasta_paths_to_run:
        alphafold_path = os.path.join(Path(__file__).parent.parent, 'alphafold', 'run_alphafold.py')
        alphafold_command = f'python {alphafold_path}'
        for key, value in alphafold_input.items():
            if key == 'fasta_paths':
                value = ','.join(fasta_paths_to_run)
            alphafold_command += f' --{key}={str(value)}'
        subprocess.run(alphafold_command, shell=True, check=True)

    # Retrieve the model path for both the wild-type and variant predictions
    model_paths = []
    for fasta_path in alphafold_input['fasta_paths']:
        model_path = glob.glob(model_path_patterns[fasta_path])[0]
        model_paths.append(model_path)

    return model_paths


def run_ring(model_paths: list[str]) -> list[str]:
    """
    Runs RING for the wild-type and variant models and returns the 
    calculated bonds (edges) file path for both.

    Args:
        model_paths: List of model paths (wild-type then variant).
    
    Returns:
        List of bonds paths (wild-type then variant).
    """
    # Define the RING path
    ring_path = os.path.join(Path(__file__).parent.parent, 'ring', 'out', 'bin', 'ring')
    
    # Define the RING command to run for the wild-type and variant models
    bonds_paths = []
    for model_path in model_paths:
        output_subdir = Path(model_path).parent
        ring_command = (
            f'{ring_path} '
            f'-i {model_path} '
            f'--out_dir {output_subdir} '
            '--no_add_H '
            '--all_edges '
            '--relaxed'
        )
        
        # Run the RING command for the current model
        subprocess.run(ring_command, shell=True, check=True)
        
        # Retrieve the bonds path for the current model
        bonds_path_pattern = os.path.join(output_subdir, '*.pdb_ringEdges')
        bonds_path = glob.glob(bonds_path_pattern)[0]
        bonds_paths.append(bonds_path)
        
    return bonds_paths


def calculate_alpharing_score(
    variable_residue: int,
    model_paths: list[str], 
    bonds_paths: list[str],
    bond_weight_formulas: dict[str, callable]
) -> None:
    """
    Calculates the AlphaRING score for the variant and outputs it into 
    alpharing_score.txt.

    Args:
        variable_residue: Variable residue position between wild-type and variant.
        model_paths: List of model paths (wild-type then variant).
        bonds_paths: List of bonds paths (wild-type then variant).
        bond_weight_formulas: Dictionary of bond weighting formulas.

    Returns:
        None
    """
    # Filter for bonds that involve the variable residue in the wild-type and variant
    residue_weights = []
    for bonds_path in bonds_paths:
        bonds = pd.read_csv(bonds_path, sep='\t')
        bonds['BondType'] = bonds['Interaction'].str.split(':').str[0]
        bonds['Residue1'] = bonds['NodeId1'].str.split(':').str[1].astype(int)
        bonds['Residue2'] = bonds['NodeId2'].str.split(':').str[1].astype(int)
        relevant_bonds = bonds[
            (bonds['Residue1'] == variable_residue) | 
            (bonds['Residue2'] == variable_residue)
        ].copy()

        # Calculate the weight of each relevant bond in the current protein
        relevant_bonds['Weight'] = relevant_bonds.apply(
            lambda bond: bond_weight_formulas.get(bond['BondType'], lambda _: pd.NA)(bond), 
            axis=1
        )
        
        # Calculate the weight of variable residue in the current protein
        residue_weight = relevant_bonds['Weight'].sum(skipna=True) + 1
        residue_weights.append(residue_weight)
    wildtype_weight, variant_weight = residue_weights

    # Extract the wildtype pLDDT of the variable residue
    wildtype_model = PDBParser().get_structure('', model_paths[0])[0]['A']
    wildtype_plddt = None
    for residue in wildtype_model:
        residue_position = residue.get_id()[1]
        if residue_position == variable_residue:
            wildtype_plddt = residue.child_list[0].get_bfactor()
            break
    
    # Calculate the AlphaRING score for the variant
    alpharing_score = abs(math.log2(variant_weight / wildtype_weight)) - 3 * math.log2(100 - wildtype_plddt)

    # Output the AlphaRING score into alpharing_score.txt
    variant_output_subdir = Path(model_paths[1]).parent
    alpharing_score_path = os.path.join(variant_output_subdir, 'alpharing_score.txt')
    with open(alpharing_score_path, 'w') as alpharing_score_file:
        alpharing_score_file.write(f'{alpharing_score}\n')
    
    return None


def main(_):
    alphafold_input = {
        'fasta_paths': 
            FLAGS.fasta_paths,
        'output_dir': 
            FLAGS.output_dir,
        'data_dir': 
            FLAGS.data_dir,
        'max_template_date': 
            FLAGS.max_template_date,
        'use_gpu_relax':
            FLAGS.use_gpu_relax,
        'uniref90_database_path': 
            os.path.join(FLAGS.data_dir, 'uniref90', 'uniref90.fasta'),
        'mgnify_database_path':
            os.path.join(FLAGS.data_dir, 'mgnify', 'mgy_clusters_2022_05.fa'),
        'bfd_database_path':
            os.path.join(FLAGS.data_dir, 'bfd', 'bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt'),
        'uniref30_database_path':
            os.path.join(FLAGS.data_dir, 'uniref30', 'UniRef30_2021_03'),
        'pdb70_database_path':
            os.path.join(FLAGS.data_dir, 'pdb70', 'pdb70'),
        'template_mmcif_dir':
            os.path.join(FLAGS.data_dir, 'pdb_mmcif', 'mmcif_files'),
        'obsolete_pdbs_path':
            os.path.join(FLAGS.data_dir, 'pdb_mmcif', 'obsolete.dat')
    }

    logging.info('(AlphaRING) Finding variable residue position')
    variable_residue = find_variable_residue(FLAGS.fasta_paths)

    logging.info('(AlphaRING) Running AlphaFold:')
    model_paths = run_alphafold(alphafold_input)

    logging.info('(AlphaRING) Running RING:')
    bonds_paths = run_ring(model_paths)

    logging.info('(AlphaRING) Calculating AlphaRING score\n')
    calculate_alpharing_score(
        variable_residue,
        model_paths,
        bonds_paths,
        BOND_WEIGHT_FORMULAS
    )

def entry_point():
    flags.mark_flags_as_required([
        'fasta_paths',
        'output_dir',
        'data_dir',
        'max_template_date',
        'use_gpu_relax'
    ])
    app.run(main)

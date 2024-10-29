#!/usr/bin/env python3

import glob
import os
import subprocess
from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
import plotly.io as pio
from absl import app, flags, logging

flags.DEFINE_list(
    'fasta_paths',
    None,
    'Paths to FASTA files, each containing a prediction target that will be '
    'folded one after another. Paths should be separated by commas. All FASTA '
    'paths must have a unique basename as the basename is used to name the output '
    'directories for each prediction.'
)

flags.DEFINE_string(
    'data_dir',
    None,
    'Path to directory of supporting data.'
)

flags.DEFINE_string(
    'output_dir',
    None,
    'Path to a directory that will store the results.'
)

flags.DEFINE_string(
    'uniref90_database_path',
    None,
    'Path to the Uniref90 database for use by JackHMMER.'
)

flags.DEFINE_string(
    'mgnify_database_path',
    None,
    'Path to the MGnify database for use by JackHMMER.'
)

flags.DEFINE_string(
    'bfd_database_path',
    None,
    'Path to the BFD database for use by HHblits.'
)

flags.DEFINE_string(
    'uniref30_database_path',
    None,
    'Path to the UniRef30 database for use by HHblits.'
)

flags.DEFINE_string(
    'pdb70_database_path',
    None,
    'Path to the PDB70 database for use by HHsearch.'
)

flags.DEFINE_string(
    'template_mmcif_dir',
    None,
    'Path to a directory with template mmCIF structures, each named <pdb_id>.cif.'
)

flags.DEFINE_string(
    'max_template_date',
    None,
    'Maximum template release date to consider. Important if folding historical '
    'test sets.'
)

flags.DEFINE_string(
    'obsolete_pdbs_path',
    None,
    'Path to file containing a mapping from obsolete PDB IDs to the PDB IDs of '
    'their replacements.'
)

flags.DEFINE_boolean(
    'use_gpu_relax',
    None,
    'Whether to relax on GPU. Relax on GPU can be much faster than CPU, so it is '
    'recommended to enable if possible. GPUs must be available if this setting is '
    'enabled.'
)

FLAGS = flags.FLAGS

FLAG_NAMES = [
    'fasta_paths',
    'data_dir',
    'output_dir',
    'uniref90_database_path',
    'mgnify_database_path',
    'bfd_database_path',
    'uniref30_database_path',
    'pdb70_database_path',
    'template_mmcif_dir',
    'max_template_date',
    'obsolete_pdbs_path',
    'use_gpu_relax'
]

WEIGHT_FORMULAS = {
    'HBOND': lambda edge: 
        edge['Energy'] * ((1 - (edge['Distance'] / 5.3)) + (edge['Angle'] / 180)),
    'IONIC': lambda edge: 
        edge['Energy'] * 2 * (1 - (edge['Distance'] / 4.5)),
    'PICATION': lambda edge: 
        edge['Energy'] * ((1 - (edge['Distance'] / 6.7)) + (1 - (edge['Angle'] / 45))),
    'PIPISTACK': lambda edge: 
        edge['Energy'] * ((1 - (edge['Distance'] / 7.3)) + (1 - (edge['Angle'] / 90))),
    'PIHBOND': lambda edge: 
        edge['Energy'] * 2 * (1 - (edge['Distance'] / 5.0))
}


def run_alphafold(flags, flag_names):
    """
    Runs AlphaFold with input and returns the best model for each FASTA.
    
    Args:
        flags: Object containing the values of the command-line flags.
        flag_names: List of the names of the command-line flags.
    
    Returns:
        List of paths to the best AlphaFold models.
    """
    # Define the base AlphaFold command
    alpharing_dir = Path(__file__).parent
    alphafold_path = os.path.join(alpharing_dir, 'alphafold', 'run_alphafold.py')
    alphafold_command = f'python3 {alphafold_path}'

    # Add flags to the AlphaFold command
    for flag_name in flag_names:
        flag_value = getattr(flags, flag_name)
        if flag_value is not None:
            if flag_name == 'fasta_paths':
                flag_value = ','.join(flag_value)
            else:
                flag_value = str(flag_value)
            alphafold_command += f' --{flag_name}={flag_value}'

    # Run the AlphaFold command
    subprocess.run(alphafold_command, shell=True, check=True)
    
    # Collect the best model paths
    model_paths = []
    fasta_paths = flags.fasta_paths
    output_dir = flags.output_dir
    for fasta_path in fasta_paths:
        subdir_name = Path(fasta_path).stem
        model_path_pattern = os.path.join(output_dir, subdir_name, 'relaxed*.pdb')
        model_path = glob.glob(model_path_pattern)[0]
        model_paths.append(model_path)

    return model_paths


def run_ring(model_paths):
    """
    Runs RING with the AlphaFold models and returns their edges and nodes.
    
    Args:
        model_paths: List of paths to the best AlphaFold models.
    
    Returns:
        Tuple of lists of paths to the RING edges and nodes.
    """
    # Define the RING path
    alpharing_dir = Path(__file__).parent
    ring_path = os.path.join(alpharing_dir, 'ring', 'out', 'bin', 'ring')
    
    edge_paths, node_paths = [], []
    # Define a RING command to run for each model
    for model_path in model_paths:
        output_dir = Path(model_path).parent
        ring_command = (
            f'{ring_path} '
            f'-i {model_path} '
            f'--out_dir {output_dir} '
            '--no_add_H '
            '--all_edges '
            '--relaxed'
        )
        
        # Run the RING command for the current model
        subprocess.run(ring_command, shell=True, check=True)
        
        # Collect the edge and node paths for the current model
        edge_path_pattern = f'{output_dir}/*.pdb_ringEdges'
        node_path_pattern = f'{output_dir}/*.pdb_ringNodes'
        patterns = [edge_path_pattern, node_path_pattern]
        paths_all = [edge_paths, node_paths]
        for pattern, paths in zip(patterns, paths_all):
            path = glob.glob(pattern)[0]
            paths.append(path)
        
    return edge_paths, node_paths


def calculate_edge_weights(edge_paths, weight_formulas):
    """
    Calculates the weights of the RING edges and return them.
    
    Args:
        edge_paths: List of paths to the RING edge files.
        weight_formulas: Dictionary mapping edge types to their weight formulas.
    
    Returns:
        List of RING edge pandas DataFrames with weights.
    """
    edges_all = []
    # Calculate the edge weights for each edges file
    for edge_path in edge_paths:
        edges = pd.read_csv(edge_path, sep='\t')
        edge_types = [interaction.split(':')[0] for interaction in edges['Interaction']]
        edges['Weight'] = [
            weight_formulas.get(edge_type, lambda _: pd.NA)(edge)
            for edge, edge_type in zip(edges.to_dict('records'), edge_types)
        ]
        
        # Save the current updated edges file
        edges.to_csv(edge_path, sep='\t', index=False)
        
        # Collect the current updated edges file
        edges_all.append(edges)
    
    return edges_all


def calculate_node_weights(node_paths, edges_all):
    """
    Calculates the weights of the RING nodes and return them.
    
    Args:
        node_paths: List of paths to the RING node files.
        edges_all: List of RING edge pandas DataFrames with weights.
        
    Returns:
        List of RING node pandas DataFrames with weights
    """
    nodes_all = []
    # Calculate the node weights for each nodes file
    for node_path, edges in zip(node_paths, edges_all):
        nodes = pd.read_csv(node_path, sep='\t', index_col='NodeId')
        nodes['Weight'] = nodes.index.map(
            edges.groupby('NodeId1')['Weight'].sum()
            .add(edges.groupby('NodeId2')['Weight'].sum(), fill_value=0)
        )
        nodes.reset_index(inplace=True)
        
        # Save the current updated nodes file
        nodes.to_csv(node_path, sep='\t', index=False)
        
        # Collect the current updated nodes file
        nodes_all.append(nodes)
        
    return nodes_all


def plot_node_weights(nodes_all, model_paths):
    """
    Plots the node weights for each model.
    
    Args:
        nodes_all: List of RING node pandas DataFrames with weights.
        model_paths: List of paths to the best AlphaFold models.
        
    Returns:
        None (but saves each node weight plot to file).
    """
    # Prepare the node weight plots for each model
    for nodes, model_path in zip(nodes_all, model_paths):
        node_total = len(nodes)
        node_names = [
            f"{node_id.split(':')[3][0]}"
            f"{node_id.split(':')[3][1:].lower()}"
            f"{node_id.split(':')[1]}"
            for node_id in nodes['NodeId']
        ]
        node_numbers = nodes['NodeId'].str.split(':').str[1].astype(int)
        maximum_weight = nodes['Weight'].max()
        minimum_weight = nodes['Weight'].min()
        normalised_weights = (
            [0.5] * node_total if maximum_weight == minimum_weight 
            else (nodes['Weight'] - minimum_weight) / (maximum_weight - minimum_weight)
        )
        
        # Create a node weight plot
        figure = go.Figure()
        figure.add_trace(go.Scatter(
            x=node_numbers,
            y=nodes['Weight'],
            mode='markers',
            marker=dict(color=normalised_weights, colorscale='Reds', line=dict(width=0.25)),
            text = [
                f'{node_name}<br>{node_weight:.2f}' 
                for node_name, node_weight in zip(node_names, nodes['Weight'])
            ],
            hoverinfo='text'
        ))
        
        # Set the title and axis labels of the node weight plot
        output_dir_name = Path(model_path).parent.name
        figure_title = f'{output_dir_name}'
        figure.update_layout(
            title={'text': figure_title},
            xaxis={'title': {'text': 'Residue'}},
            yaxis={'title': {'text': 'Weight'}}
        )
        
        # Add a normalised color bar to the node weight plot
        for index, normalised_weight in enumerate(normalised_weights):
            color = px.colors.sample_colorscale('Reds', normalised_weight)[0]
            figure.add_shape(
                type="rect",
                x0=index + 1,
                x1=index + 2,
                yref='paper',
                y0=0,
                y1=0.025,
                line=dict(width=0),
                fillcolor=color,
                opacity=1
            )
        
        # Save the node weight plot to file
        model_name = Path(model_path).name
        figure_name = f'{model_name}_weightPlot.html'
        output_dir = Path(model_path).parent
        figure_path = os.path.join(output_dir, figure_name)
        pio.write_html(figure, figure_path)


def main(_):
    logging.info('(AlphaRING) Running AlphaFold:')
    model_paths = run_alphafold(FLAGS, FLAG_NAMES)

    logging.info('(AlphaRING) Running RING:')
    edge_paths, node_paths = run_ring(model_paths)

    logging.info('(AlphaRING) Calculating edge weights')
    edges_all = calculate_edge_weights(edge_paths, WEIGHT_FORMULAS)
    
    logging.info('(AlphaRING) Calculating node weights')
    nodes_all = calculate_node_weights(node_paths, edges_all)
    
    logging.info('(AlphaRING) Plotting node weights')
    plot_node_weights(nodes_all, model_paths)


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'fasta_paths',
        'data_dir',
        'output_dir',
        'uniref90_database_path',
        'mgnify_database_path',
        'bfd_database_path',
        'uniref30_database_path',
        'pdb70_database_path',
        'template_mmcif_dir',
        'max_template_date',
        'obsolete_pdbs_path',
        'use_gpu_relax'
    ])

    app.run(main)

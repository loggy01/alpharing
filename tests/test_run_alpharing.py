"""Integration test for run_alpharing.py"""

import os
import shutil
import subprocess
from pathlib import Path
from unittest import mock

import pytest
import alpharing.run_alpharing as alpharing
import pandas as pd
import pandas.testing as pdt
from absl import flags


def test_end_to_end(tmp_path):
    test_data_dir = os.path.join(Path(__file__).parent, 'test_data')
    temp_test_data_dir = os.path.join(tmp_path, "test_data")
    shutil.copytree(test_data_dir, temp_test_data_dir)
    fasta_path = os.path.join(temp_test_data_dir, 'input', 'protein.fa')
    output_dir = os.path.join(temp_test_data_dir, 'output')

    with mock.patch.object(subprocess, 'run') as mock_subprocess:
        mock_subprocess.return_value = None
        flags.FLAGS.unparse_flags()
        flags.FLAGS([
            'alpharing',
            f'--fasta_path={fasta_path}',
            '--substitutions=YA229S,VA194A,TA188Q',
            f'--output_dir={output_dir}',
            '--data_dir=dummy_data_dir',
            '--ring_exe_path=dummy_ring_exe_path',
            '--foldx_exe_path=dummy_foldx_exe_path'
        ])
        with pytest.raises(SystemExit):
            alpharing.entry_point()

    results_df_path = os.path.join(output_dir, 'protein', 'alpharing_scores.txt')
    results_df = pd.read_csv(results_df_path, sep='\t')
    expected_results_df = pd.DataFrame({
        "Substitution": ["YA229S", "VA194A", "TA188Q"],
        "pLDDT": [89.93, 85.72, 94.34],
        "Degree": [12, 6, 15],
        "ΔΔG": [2.25407, 1.09669, -0.238198],
        "RSP": [0.7762711864406779,
                0.6576271186440678,
                0.6372881355932203],
        "label": ["Deleterious", "Deleterious", "Deleterious"],
        "Probability": [0.7206843495368958,
                        0.3487226963043213,
                        0.7832099795341492],
        "pLDDT SHAP": [0.014773664648334156,
                       -0.04466838430613285,
                       0.11913539444406826],
        "Degree SHAP": [0.11654547924796732,
                        -0.05072184357792139,
                        0.22942934321860467],
        "ΔΔG SHAP": [0.07618079070001843,
                     -0.07727733583499988,
                     -0.12450742968668538],
        "RSP SHAP": [-0.07628489354004467,
                     -0.06807904845724505,
                     -0.030316636922458784],
    })

    pdt.assert_frame_equal(results_df, expected_results_df, atol=1e-3)

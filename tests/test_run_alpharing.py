"""Integration test for run_alpharing.py"""

import shutil
import subprocess
from pathlib import Path
from unittest import mock

import pytest
import alpharing.run_alpharing as run_alpharing
from absl import flags


def test_end_to_end(tmp_path):
    test_data_dir = Path(__file__).parent / 'test_data'
    temp_test_data_dir = tmp_path / "test_data"
    shutil.copytree(test_data_dir, temp_test_data_dir)

    wildtype_fasta_path = temp_test_data_dir / 'input' / 'NM_007294.4_BRCA1.fa'
    variant_fasta_path = temp_test_data_dir / 'input' / 'NM_007294.4_BRCA1_Arg170Trp_Benign.fa'
    output_dir = temp_test_data_dir / 'output'

    with mock.patch.object(subprocess, 'run') as mock_subprocess:
        mock_subprocess.return_value = None
        flags.FLAGS.unparse_flags()
        flags.FLAGS([
            'run_alpharing',
            f'--fasta_paths={wildtype_fasta_path},{variant_fasta_path}',
            f'--output_dir={output_dir}',
            '--data_dir=dummy_data_dir',
            '--max_template_date=2050-01-01',
            '--use_gpu_relax=False'
        ])
        with pytest.raises(SystemExit):
            run_alpharing.entry_point()

    alpharing_score_path = output_dir / 'NM_007294.4_BRCA1_Arg170Trp_Benign' / 'alpharing_score.txt'
    with open(alpharing_score_path) as alpharing_score_file:
        alpharing_score = float(alpharing_score_file.read())
    assert alpharing_score == pytest.approx(-18.098286221441384)

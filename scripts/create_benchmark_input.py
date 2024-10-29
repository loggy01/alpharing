#!/usr/bin/env python3


### Part 1: Setup ###

import os
import pandas as pd
import re
import requests

vt_sum_path = './variant_summary.txt' # Downloaded from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
wt_fa_dir = './wt_fastas/' # create this directory
vt_fa_dir = './vt_fastas/' # create this directory

resis = {
    'Ala': 'A',
    'Arg': 'R',
    'Asn': 'N',
    'Asp': 'D',
    'Cys': 'C',
    'Glu': 'E',
    'Gln': 'Q',
    'Gly': 'G',
    'His': 'H',
    'Ile': 'I',
    'Leu': 'L',
    'Lys': 'K',
    'Met': 'M',
    'Phe': 'F',
    'Pro': 'P',
    'Ser': 'S',
    'Thr': 'T',
    'Trp': 'W',
    'Tyr': 'Y',
    'Val': 'V',
}


### Part 2: Filter for missense variants ###

vt_sum = pd.read_csv(
    vt_sum_path,
    sep='\t',
    dtype={'Chromosome': str}
)

snv_vts = vt_sum[
    (vt_sum['Type'] == 'single nucleotide variant') &
    (vt_sum['ReviewStatus'] == 'reviewed by expert panel') &
    (vt_sum['Assembly'] == 'GRCh38') &
    ((vt_sum['ClinicalSignificance'] == 'Pathogenic') | 
     (vt_sum['ClinicalSignificance'] == 'Benign'))
]

mis_vt_pat = re.compile(
    r'(' + '|'.join(resis.keys()) + r')\d+(' + '|'.join(resis.keys()) + r')'
)

mis_vts = snv_vts[
    snv_vts['Name'].str.extract(mis_vt_pat).notnull().all(axis=1)
].copy()


### Part 3: Create dataframes for unique missense variants and unique IDs ###

mis_vts['NuccoreAccession'] = mis_vts['Name'].str.split('(').str[0]

mis_vts['ResidueSubstitution'] = [
    mis_vt_pat.search(name).group(0) for name in mis_vts['Name'] 
]

unq_mis_vts = mis_vts.drop_duplicates(
    subset=['NuccoreAccession', 'ResidueSubstitution'] 
)

unq_mis_vt_ids = unq_mis_vts.drop_duplicates(subset=['NuccoreAccession'])


### Part 4: Use unique IDs df to generate wildtype FASTAs ###

ids_to_remove = set()

for id, genes in zip(
    unq_mis_vt_ids['NuccoreAccession'], 
    unq_mis_vt_ids['GeneSymbol']
):
    gene = genes.split(';')[0]
    
    efetch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    efetch_params = {
        'db': 'nuccore',
        'id': id,
        'rettype': 'gb',
        'retmode': 'text'
    }
    efetch_output = requests.get(efetch_url, params=efetch_params)
    
    gb_data = efetch_output.text
    
    cds = gb_data.find('CDS')
    prot_start = gb_data.find('/translation=', cds)
    prot_start += len('/translation=') + 1
    prot_end = gb_data.find('"', prot_start)
    prot = gb_data[prot_start:prot_end]
    prot = prot.replace('\n', '').replace(' ', '')
    
    if len(prot) < 16 or len(prot) > 2700:
        ids_to_remove.add(id)
        continue
    
    wt_fa_path = os.path.join(wt_fa_dir, f'{id}_{gene}.fa')
    with open(wt_fa_path, 'w') as wt_fa_file:
        wt_fa_file.write(f'>{id}_{gene}\n')
        wt_fa_file.write(f'{prot}\n')

unq_mis_vts = unq_mis_vts[
    ~unq_mis_vts['NuccoreAccession'].isin(ids_to_remove)
]


### Part 5: Use unique missense variants df to generate variant FASTAs ###

for id, genes, resi_sub, status in zip(
    unq_mis_vts['NuccoreAccession'],
    unq_mis_vts['GeneSymbol'],
    unq_mis_vts['ResidueSubstitution'],
    unq_mis_vts['ClinicalSignificance']
):  
    genes = genes.split(';')
    gene_match = None
    
    for gene in genes:
        wt_fa_path = os.path.join(wt_fa_dir, f'{id}_{gene}.fa')
        if os.path.exists(wt_fa_path):
            gene_match = gene
            break
    
    vt_resi = resis[resi_sub[-3:]]
    resi_pos = int(re.search(r'\d+', resi_sub).group()) - 1
    
    with open(wt_fa_path, 'r') as wt_fa_file:
        wt_fa = wt_fa_file.read()
    
    wt_seq = wt_fa.split('\n', 1)[1]
    vt_seq = wt_seq[:resi_pos] + vt_resi + wt_seq[resi_pos + 1:]
        
    vt_fa_path = os.path.join(vt_fa_dir, f'{id}_{gene_match}_{resi_sub}_{status}.fa')
    with open(vt_fa_path, 'w') as vt_fa_file:
        vt_fa_file.write(f'>{id}_{gene_match}_{resi_sub}_{status}\n')
        vt_fa_file.write(vt_seq)

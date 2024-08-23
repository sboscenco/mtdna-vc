#!/usr/bin/env python3

import os
import argparse
import re
import pandas as pd
import numpy as np
from collections import defaultdict

def create_mut_sig_file(resultsdir, tumor_id):

    equivalence = {
    'GA': 'CT', 
    'GC': 'CG',  
    'GT': 'CA',  
    'AC': 'TG',  
    'AG': 'TC',  
    'AT': 'TA'  
}
    
    filloutfile = pd.read_csv(f"{resultsdir}/{tumor_id}.bam.maf",sep = '\t',header = 0, comment='#', low_memory=False)
    SNVs = filloutfile[filloutfile["Variant_Type"] == "SNP"]
    # nested dict to hold motif counts
    motif_counts = defaultdict(lambda: defaultdict(int))

    for _, row in SNVs.iterrows():
        motif = row["flanking_bps"]
        ref = row["Reference_Allele"]
        alt = row["Tumor_Seq_Allele2"]
        mutation = f'{ref}{alt}'

        mut_equiv = equivalence.get(mutation, mutation)
        motif_counts[motif][mut_equiv] += 1

    mutsigfile = pd.DataFrame(motif_counts).fillna(0).astype(int)

    mutsigfile.to_csv(resultsdir + "/" + tumor_id + '_mutsig.tsv',sep = '\t')



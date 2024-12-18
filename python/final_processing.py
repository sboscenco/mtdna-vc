import os
import numpy as np
import pandas as pd
import argparse
import pysam
import subprocess

def process_indelmaf(outdir, tumorbam, normalbam):
    maf = pd.read_csv(outdir + "/" + tumorbam + '.maf',header = 0,sep = '\t',comment = '#')
    if(normalbam != ""):
        maf['N_AltFwd'] = np.nan
        maf['N_AltRev'] = np.nan

    maf['t_F1R2'] = '0,0'
    maf['t_F2R1'] = '0,0'
    print(maf.head())
    maf[['T_RefFwd', 'T_AltFwd']] = maf['t_F1R2'].str.split(',', expand=True)
    maf[['T_RefRev', 'T_AltRev']] = maf['t_F2R1'].str.split(',', expand=True)
    
    maf = maf[~maf['Start_Position'].isin(list(range(3106, 3107)))]
    maf.to_csv(outdir + tumorbam + ".maf",index = None,sep = '\t')

def process_maf(outdir, workingdir, tumorbam, normalbam = "", indel = False):
   
    # Read in files
    maf = pd.read_csv(outdir + "/" + tumorbam + '.maf',header = 0,sep = '\t',comment = '#')
    genepos = pd.read_csv(workingdir + '/reference/GenePositions_imported.csv', header = 0, index_col = 0)
    mitotip = pd.read_csv(workingdir + '/reference/All_combinations.csv', header = 0, sep = ",")
    apogee2 = pd.read_csv(workingdir + '/reference/MitImpact_db_3.1.2.txt', header = 0, sep = "\t", low_memory = False)

    maf.rename(columns = {"t_depth": "T_TotalDepth", "n_depth": "N_TotalDepth", "t_ref_count": "T_RefCount", "t_alt_count": "T_AltCount",
                          "n_ref_count": "N_RefCount", "n_alt_count": "N_AltCount"}, inplace = True, errors = "ignore")
    
    # modify the gene names to include rRNA and tRNA and change the control region symbols
    maf['Hugo_Symbol'] = genepos.loc[maf['Start_Position'],'Gene'].reset_index(drop = True)
    maf.loc[(maf['Start_Position'].map(int) <= 576) | (maf['Start_Position'].map(int) >= 16024),'Hugo_Symbol'] = 'ControlRegion'
    maf.loc[(maf['Start_Position'].map(int) <= 5798) & (maf['Start_Position'].map(int) >= 5721),'Hugo_Symbol'] = 'MT-OLR'
    
    list_cols = ["Hugo_Symbol", "Chromosome", "NCBI_Build", "Start_Position", "End_Position", "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele",
               "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2", "Gene",
               "Consequence", "HGVSc", "HGVSp", "HGVSp_Short", 'flanking_bps',
               "T_TotalDepth", "T_RefCount", "T_AltCount", "T_AltFwd", "T_AltRev"]

    # subsetting columns of interest
    if(normalbam != ""): list_cols += ["Normal_Sample_Barcode", "N_TotalDepth", "N_RefCount", "N_AltCount", 'N_AltFwd', 'N_AltRev']
    maf = maf[list_cols]

    # create short variant id
    maf['ShortVariantID'] = maf['Reference_Allele'] + maf['Start_Position'].map(str) + maf['Tumor_Seq_Allele2']

    # Compute Heteroplasmy 
    maf['TumorVAF'] = maf['T_AltCount']/maf['T_TotalDepth']
    if("N_TotalDepth" in maf.columns):
        maf['NormalVAF'] = maf['N_AltCount']/maf['N_TotalDepth']
        maf["N_TotalDepth"] = maf["N_TotalDepth"].fillna(0)
        maf['NormalVAF'] = maf['NormalVAF'].fillna(0)
        maf["N_AltCount"] = maf["N_AltCount"].fillna(0)
        maf["N_RefCount"] = maf["N_RefCount"].fillna(0)

        # FILTER GERMLINE VAFS
        germline_vars = (maf["N_TotalDepth"] >= 5) & (maf["NormalVAF"] >= 0.5)
        maf = maf[~germline_vars]
        
    else:
        maf['NormalVAF'] = np.nan

    # add mitotip pathogenicity
    mitotip.columns = ["Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "MitoTIP_Score"] 
    apogee2 = apogee2[["Start", "Ref", "Alt", "APOGEE2", "Respiratory_Chain_complex", "HelixMTdb_AF_hom", "HelixMTdb_AF_het"]]
    apogee2.rename(columns = {"Start": "Start_Position", "Alt": "Tumor_Seq_Allele2", "Respiratory_Chain_complex": "Complex", "Ref": "Reference_Allele"}, inplace = True)

    # merge files
    maf = pd.merge(maf, apogee2, how = "left", on = ["Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"])
    maf = pd.merge(maf, mitotip, how = "left", on = ["Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"])
   
    # mitotip pathogeneicity annotations based on scores 
    maf["MitoTIP"] = "unknown"
    maf["MitoTIP"] = maf["MitoTIP"].astype(object)
    maf.loc[maf["MitoTIP_Score"] > 16.25, "MitoTIP"] = "likely_pathogenic"
    maf.loc[(maf["MitoTIP_Score"] > 12.66) & (maf["MitoTIP_Score"] <= 16.25), "MitoTIP"] = "possibly_pathogenic"
    maf.loc[(maf["MitoTIP_Score"] > 8.44) & (maf["MitoTIP_Score"] <= 12.66), "MitoTIP"] = "possibly_benign"
    maf.loc[maf["MitoTIP_Score"] <= 8.44] = "likely_benign"

    # write out final files
    if(indel):
         maf.to_csv(outdir + "/" + tumorbam[:-13] + 'indel.bam.maf', index = None, sep = '\t')
    else:
         maf.to_csv(outdir + "/" + tumorbam + '.maf', index = None, sep = '\t')
    
   
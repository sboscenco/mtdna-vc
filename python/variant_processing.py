#!/usr/bin/env python3

import os
import argparse
import re
import pandas as pd
import numpy as np
import subprocess
#from Bio import SeqIO
#from Bio.Seq import Seq

def variant_processing(tumor_id,resultsdir):
    """
    Run MTvariantpipeline and MuTect2 on the filtered cells
    MTvariantpipeline: A simple variant calling and annotation pipeline for mitochondrial DNA variants.
    """
    print("Starting variant processing...")

    # Read in MTvariantpipeline result
    MTvarfile = pd.read_csv(resultsdir + "/MTvariant_results/" + tumor_id + ".bam.maf", sep = "\t", comment='#', low_memory=False)
    # Read in MuTect result
    mutectfile = pd.read_csv(resultsdir + "/MuTect2_results/" + tumor_id + ".bam.maf", sep = "\t", comment='#', header=0, low_memory=False)

    # Filter out variants falling in the repeat regions of 513-525, and 3105-3109 (black listed regions)
    mutectfile.drop(mutectfile[((mutectfile['Start_Position'].between(513, 525)) | (mutectfile['Start_Position'].between(3105,3109)))].index, inplace = True)

    MTvarfile.drop(MTvarfile[((MTvarfile['Start_Position'].between(513, 525)) | (MTvarfile['Start_Position'].between(3105,3109)))].index, inplace = True)
    MTvarfile.drop(MTvarfile[((MTvarfile['End_Position'].between(513, 525)) | (MTvarfile['End_Position'].between(3105,3109)))].index, inplace = True)

    mutectfile.reset_index(drop=True, inplace=True)
    MTvarfile.reset_index(drop=True, inplace=True)

    # Output the overlap as final maf file
    combinedfile = pd.merge(mutectfile, MTvarfile, how='inner', on=['Chromosome','Start_Position',
        'Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Variant_Type'])
    
    print(combinedfile.columns)
    # Filter out columns that are non-useful
    keep_cols = ["Hugo_Symbol", "Entrez_Gene_Id", "NCBI_Build", "Chrom", "Start_Position", "End_Position", "VariantClass", "Variant_Type", "Ref", "Alt", 
                     "dbSNP_RS", "Sample", "NormalUsed", "T_TotalDepth", "T_RefCount", "T_AltCount", "N_TotalDepth", "N_RefCount", "N_AltCount"
                     "T_AltFwd", "T_AltRev", "Gene", "Feature", "Feature_type", "Consequence", "flanking_bps", "Heteroplasmy", "ShortVariantID", "TumorVAF", "NormalVAF"]
    
    #combinedfile = combinedfile[keep_cols]

    #print(combinedfile)
    """
    
    # Fix INDELs in the same position i.e. A:11866:AC and A:11866:ACC
    aux = combinedfile.loc[combinedfile['Variant_Type'] == 'INS'].groupby('Start_Position').count()['Hugo_Symbol'].reset_index()
    positions = list(aux['Start_Position'].loc[aux['Hugo_Symbol'] > 1])
    variants = list(combinedfile['ShortVariantID'].loc[(combinedfile['Start_Position'].isin(positions)) & (combinedfile['Variant_Type'] == 'INS')])
    if len(positions) != 0:
        dff = combinedfile.loc[combinedfile['ShortVariantID'].isin(variants)]
        # Create an auxuliary file only with the last rows to keep: keep unique positions with the highest TumorVAF
        dffaux = dff.sort_values(by='TumorVAF', ascending = False)
        dffaux = dffaux.drop_duplicates('Start_Position', keep = 'first')
        for i in positions:
            vals = dff[['t_alt_count_y', 't_alt_count_x']].loc[dff['Start_Position'] == i].sum(axis = 0).reset_index()
            dvals = dict(zip(list(vals['index']),list(vals[0])))
            dffaux.loc[dffaux['Start_Position'] == i,'t_alt_count_y'] = dvals['t_alt_count_y']
            dffaux.loc[dffaux['Start_Position'] == i,'t_alt_count_x'] = dvals['t_alt_count_x']
        #Remove all variants with duplicated indels
        combinedfile = combinedfile.loc[(~combinedfile['ShortVariantID'].isin(variants))]
        # Add unique indel variants with new values
        combinedfile = pd.concat([combinedfile, dffaux])
        combinedfile = combinedfile.sort_values(by='Start_Position', ascending = True)
        
    # Final annotation
    filloutfile = combinedfile.loc[:,['Hugo_Symbol','Entrez_Gene_Id','Center','NCBI_Build','Chromosome',
        'Start_Position','End_Position','Strand','Variant_Classification','Variant_Type','Reference_Allele',
        'Tumor_Seq_Allele1','Tumor_Seq_Allele2','dbSNP_RS','dbSNP_Val_Status','Tumor_Sample_Barcode',
        'Matched_Norm_Sample_Barcode','Match_Norm_Seq_Allele1','Match_Norm_Seq_Allele2','Tumor_Validation_Allele1',
        'Tumor_Validation_Allele2','Match_Norm_Validation_Allele1','Match_Norm_Validation_Allele2','Verification_Status',
        'Validation_Status','Mutation_Status','Sequencing_Phase','Sequence_Source','Validation_Method','Score','BAM_File',
        'Sequencer','Tumor_Sample_UUID','Matched_Norm_Sample_UUID','HGVSc','HGVSp','HGVSp_Short','Exon_Number','t_depth',
        't_ref_count','t_alt_count','n_depth','n_ref_count','n_alt_count','t_alt_fwd','t_alt_rev','all_effects','Gene',
        'Feature','Feature_type','Consequence','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons',
        'Existing_variation','ALLELE_NUM','DISTANCE','STRAND_VEP','SYMBOL','SYMBOL_SOURCE','HGNC_ID','BIOTYPE','CANONICAL',
        'CCDS','ENSP','SWISSPROT','TREMBL','UNIPARC','RefSeq','SIFT','PolyPhen','EXON','INTRON','DOMAINS','AF','AFR_AF',
        'AMR_AF','ASN_AF','EAS_AF','EUR_AF','SAS_AF','AA_AF','EA_AF','CLIN_SIG','SOMATIC','PUBMED','MOTIF_NAME',
        'MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE','IMPACT','PICK','TSL','HGVS_OFFSET','PHENO','MINIMISED',
        'GENE_PHENO','FILTER','flanking_bps','vcf_id','vcf_qual','gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF',
        'gnomAD_ASJ_AF','gnomAD_EAS_AF','gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_OTH_AF','gnomAD_SAS_AF']]
    filloutfile.columns = ['Hugo_Symbol','Entrez_Gene_Id','Center','NCBI_Build','Chrom','Start','End_Position','Strand',
        'VariantClass','Variant_Type','Ref','Tumor_Seq_Allele1','Alt','dbSNP_RS','dbSNP_Val_Status','Sample','NormalUsed',
        'Match_Norm_Seq_Allele1','Match_Norm_Seq_Allele2','Tumor_Validation_Allele1','Tumor_Validation_Allele2',
        'Match_Norm_Validation_Allele1','Match_Norm_Validation_Allele2','Verification_Status','Validation_Status',
        'Mutation_Status','Sequencing_Phase','Sequence_Source','Validation_Method','Score','BAM_File','Sequencer',
        'Tumor_Sample_UUID','Matched_Norm_Sample_UUID','HGVSc','HGVSp','HGVSp_Short','Exon_Number','T_TotalDepth',
        'T_RefCount','T_AltCount','N_TotalDepth','N_RefCount','N_AltCount','T_AltFwd','T_AltRev','all_effects','Gene',
        'Feature','Feature_type','Consequence','cDNA_position','CDS_position','Protein_position','Amino_acids','Codons',
        'Existing_variation','ALLELE_NUM','DISTANCE','STRAND_VEP','SYMBOL','SYMBOL_SOURCE','HGNC_ID','BIOTYPE','CANONICAL',
        'CCDS','ENSP','SWISSPROT','TREMBL','UNIPARC','RefSeq','SIFT','PolyPhen','Exon','INTRON','DOMAINS','AF','AFR_AF',
        'AMR_AF','ASN_AF','EAS_AF','EUR_AF','SAS_AF','AA_AF','EA_AF','CLIN_SIG','SOMATIC','PUBMED','MOTIF_NAME','MOTIF_POS',
        'HIGH_INF_POS','MOTIF_SCORE_CHANGE','IMPACT','PICK','TSL','HGVS_OFFSET','PHENO','MINIMISED','GENE_PHENO','FILTER',
        'flanking_bps','vcf_id','vcf_qual','gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_ASJ_AF','gnomAD_EAS_AF',
        'gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_OTH_AF','gnomAD_SAS_AF']
    filloutfile.index = [str(filloutfile['Ref'][i]) + ':' + str(int(filloutfile['Start'][i])) + ':' + 
                         str(filloutfile['Alt'][i]) for i in range(len(filloutfile))]
    
    # Obtain the mutation signature
    # Initialize the counts and mutation sigature matrix
    motifs_C = ["ACA","ACC","ACG","ACT","CCA","CCC","CCG","CCT","GCA","GCC","GCG","GCT","TCA","TCC","TCG","TCT"]
    motifs_T = ["ATA","ATC","ATG","ATT","CTA","CTC","CTG","CTT","GTA","GTC","GTG","GTT","TTA","TTC","TTG","TTT"]
    mutsigfile = pd.DataFrame(index=['counts_CA','counts_CG','counts_CT','counts_TA','counts_TC','counts_TG'], columns=range(16))
    mutsigfile = mutsigfile.fillna(0)
    counts_CA = np.zeros(16)
    counts_CG = np.zeros(16)
    counts_CT = np.zeros(16)
    counts_TA = np.zeros(16)
    counts_TC = np.zeros(16)
    counts_TG = np.zeros(16)

    # Import the reference fasta file
    fasta_sequences = SeqIO.parse(open(reffile),'fasta')
    for fasta in fasta_sequences:
        currheader, currsequence = fasta.id, fasta.seq
        if 'MT' in currheader:
            sequence = [base for base in currsequence]
        if 'chrM' in currheader:
            sequence = [base for base in currsequence]
    varref = [variants[0] for variants in pd.Series(filloutfile.index.values).str.split(':')]
    varpos = [variants[1] for variants in pd.Series(filloutfile.index.values).str.split(':')]
    varalt = [variants[2] for variants in pd.Series(filloutfile.index.values).str.split(':')]
    mutsigmotifs = []
    for eachone in range(len(varpos)):
        prevpos = int(varpos[eachone])-2
        currpos = int(varpos[eachone])-1
        nextpos = int(varpos[eachone])
        motif = ''.join([sequence[prevpos],sequence[currpos],sequence[nextpos]])
        mutsigmotifs.append(motif)
        if varref[eachone] == 'C':
            if varalt[eachone] == 'A':
                counts_CA[motifs_C.index(motif)] += 1
            elif varalt[eachone] == 'G':
                counts_CG[motifs_C.index(motif)] += 1
            elif varalt[eachone] == 'T':
                counts_CT[motifs_C.index(motif)] += 1
        elif varref[eachone] == 'T':
            if varalt[eachone] == 'A':
                counts_TA[motifs_T.index(motif)] += 1
            elif varalt[eachone] == 'C':
                counts_TC[motifs_T.index(motif)] += 1
            elif varalt[eachone] == 'G':
                counts_TG[motifs_T.index(motif)] += 1
        elif varref[eachone] == 'G':
            motif = str(Seq(motif).complement())
            if varalt[eachone] == 'A':
                counts_CT[motifs_C.index(motif)] += 1
            elif varalt[eachone] == 'C':
                counts_CG[motifs_C.index(motif)] += 1
            elif varalt[eachone] == 'T':
                counts_CA[motifs_C.index(motif)] += 1
        elif varref[eachone] == 'A':
            motif = str(Seq(motif).complement())
            if varalt[eachone] == 'C':
                counts_TG[motifs_T.index(motif)] += 1
            elif varalt[eachone] == 'G':
                counts_TC[motifs_T.index(motif)] += 1
            elif varalt[eachone] == 'T':
                counts_TA[motifs_T.index(motif)] += 1
    mutsigfile.loc['counts_CA'] = counts_CA
    mutsigfile.loc['counts_CG'] = counts_CG
    mutsigfile.loc['counts_CT'] = counts_CT
    mutsigfile.loc['counts_TA'] = counts_TA
    mutsigfile.loc['counts_TC'] = counts_TC
    mutsigfile.loc['counts_TG'] = counts_TG
    
    # store the mutation signature info in variants file
    filloutfile['mutsig'] = mutsigmotifs
    # Saving the mutation signature
    mutsigfile.to_csv(resultsdir + "/" + tumor_id + '_mutsig.tsv',sep = '\t')
    # Calculate heteroplasmy
    filloutfile["Heteroplasmy"] = filloutfile['T_AltCount'].astype(int) / filloutfile['T_TotalDepth'].astype(int)
    # Combined matrices together
    filloutfile.to_csv(f"{resultsdir}/{tumor_id}.bam.maf",sep = '\t',na_rep='',index=False)
"""
variant_processing("chrMT.TC2", "~/Desktop/reznik/hcc/vc/TC2/")
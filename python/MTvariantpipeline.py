# imports
import os
import numpy as np
import pandas as pd
import argparse
import pysam
import subprocess
import final_processing
from Bio import SeqIO
from Bio.Seq import Seq
import re

## For calling snps   
def call_vars_normal(tumorbam, normalbam, mtchrom, minmapq, minbq, fasta, datadir, normaldir, vcfdir, ncbibuild):
    subprocess.run(f"bcftools mpileup --region {mtchrom} --count-orphans --no-BAQ --min-MQ {minmapq} --min-BQ {minbq} " \
                + "--ignore-RG --skip-any-set UNMAP,SECONDARY,QCFAIL,DUP --annotate DP,AD,ADF,ADR --gap-frac 0.005 " \
                + f"--tandem-qual 80 -L 1000000 -d 1000000 --open-prob 30 --skip-indels --fasta-ref {fasta} {datadir}/{tumorbam} {normaldir}/{normalbam} " \
                + f"| bcftools call --multiallelic-caller --ploidy {ncbibuild} --keep-alts -v " \
                + f"| bcftools query --format '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%DP\t%ADF\t%ADR]\n' > {vcfdir}/{tumorbam}_temp.maf", shell = True, check = True)

    print(f'Done with COUNT CALL for {tumorbam} + {normalbam}')
 
def call_vars(tumorbam, mtchrom, minmapq, minbq, fasta, datadir, vcfdir, ncbibuild):
    subprocess.run(f"bcftools mpileup --region {mtchrom} --count-orphans --no-BAQ --min-MQ {minmapq} --min-BQ {minbq} " \
                + "--ignore-RG --skip-any-set UNMAP,SECONDARY,QCFAIL,DUP --annotate DP,AD,ADF,ADR --gap-frac 0.005 " \
                + f"--tandem-qual 80 -L 1000000 -d 1000000 --open-prob 30 --skip-indels --fasta-ref {fasta} {datadir}/{tumorbam}  " \
                + f"| bcftools call --multiallelic-caller --ploidy {ncbibuild} --keep-alts -v " \
                + f"| bcftools query --format '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%DP\t%ADF\t%ADR]\n' > {vcfdir}/{tumorbam}_temp.maf", shell = True, check = True)
    
    print(f'Done with COUNT CALL for {tumorbam}')


def maf_call(tumorbam, fasta, vcfdir, workingdir, ncbibuild, vepcache, outdir):
    retaincols = ','.join(['N_AltRev', 'N_AltFwd', 'T_AltRev', 'T_AltFwd', 'Normal_Sample_Barcode'])
    
    subprocess.run(f"perl {workingdir}/vcf2maf/maf2maf.pl --vep-data {vepcache}/ --input-maf {vcfdir}/{tumorbam}_temp2.maf " \
                + f"--output-maf {outdir}/{tumorbam}.maf --retain-cols {retaincols} --ncbi-build {ncbibuild} --ref-fasta {fasta}", shell = True, check = True)
    print(f'Done with MAF call for {tumorbam}')
    

def process_tempmaf(vcfdir, tumorbam, normalbam, molecule, minstrand):
    #name cols
    tempmaf = pd.read_csv(vcfdir + tumorbam + "_temp.maf",header = None,sep = '\t')

     #remove rows with no alt reads
    tempmaf = tempmaf[tempmaf[3] != "."]

    if(normalbam != ""):
        tempmaf.columns = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "T_Allelic_Depth", "t_depth", "T_ForwardDepth", "T_ReverseDepth",
                          "N_Allelic_Depth", "n_depth", "N_ForwardDepth", "N_ReverseDepth"]
        tempmaf["Normal_Sample_Barcode"] = normalbam[:-4]

        # split up ref and alt reads 
        tempmaf[['n_ref_count', 'n_alt_count']] = tempmaf['N_Allelic_Depth'].str.split(',', expand=True)
        tempmaf[['N_RefFwd', 'N_AltFwd']] = tempmaf['N_ForwardDepth'].str.split(',', expand=True)
        tempmaf[['N_RefRev', 'N_AltRev']] = tempmaf['N_ReverseDepth'].str.split(',', expand=True)
        tempmaf.drop(columns = ["N_Allelic_Depth", "N_ForwardDepth", "N_ReverseDepth"], inplace = True)

    else:
        tempmaf.columns = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "T_Allelic_Depth", "t_depth", "T_ForwardDepth", "T_ReverseDepth"]
        tempmaf["Normal_Sample_Barcode"] = ""

    tempmaf["Tumor_Sample_Barcode"] = tumorbam[:-4]

    tempmaf[['t_ref_count', 't_alt_count']] = tempmaf['T_Allelic_Depth'].str.split(',', expand=True)
    tempmaf[['T_RefFwd', 'T_AltFwd']] = tempmaf['T_ForwardDepth'].str.split(',', expand=True)
    tempmaf[['T_RefRev', 'T_AltRev']] = tempmaf['T_ReverseDepth'].str.split(',', expand=True)
    tempmaf.drop(columns = ["T_Allelic_Depth", "T_ForwardDepth", "T_ReverseDepth"], inplace = True)
   
   #make sure enough reads supporting both forward and reverse 
   #check the rna condition!
    if molecule == 'dna':
        tempmaf = tempmaf[ (tempmaf['T_AltFwd'].map(int) >= minstrand) & (tempmaf['T_AltRev'].map(int) >= minstrand) ]
    elif molecule == 'rna':
        tempmaf = tempmaf[ (tempmaf['T_AltFwd'].map(int) >= minstrand) | (tempmaf['T_AltRev'].map(int) >= minstrand) ]

    tempmaf = tempmaf[~tempmaf['Start_Position'].isin(list(range(3106, 3107)))]

    # secondary tempoary maf file
    tempmaf.to_csv(vcfdir + tumorbam + "_temp2.maf",index = None,sep = '\t')
   

def variant_calling(datadir, tumorbam, normalbam, normaldir, vcfdir, outdir, workingdir, vepcache, fasta, minmapq, minbq, minstrand, genome, mtchrom, molecule):
    '''
    Variant calling
    '''
    print("Starting variant calling...")

    # need to test on mm10 still not sure if it will work 
    if(normalbam != ""): 
        call_vars_normal(tumorbam, normalbam, mtchrom, minmapq, minbq, fasta, datadir, normaldir, vcfdir, genome)
    else: 
        call_vars(tumorbam, mtchrom, minmapq, minbq, fasta, datadir, vcfdir, genome) 

    # process MAF
    process_tempmaf(vcfdir, tumorbam, normalbam, molecule, minstrand)
    maf_call(tumorbam, fasta, vcfdir, workingdir, genome, vepcache, outdir)
    final_processing.process_maf(outdir, workingdir, tumorbam, normalbam)

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("-d", "--datadir",type=str, help="directory for input BAM files", required=True)
    parser.add_argument("-v", "--vcfdir", type=str, help="directory for intermediate VCF files", required=True)
    parser.add_argument("-o","--outdir",type=str,help="directory for output MAF files", required=True)
    parser.add_argument("-w", "--workingdir", type=str, help="Working directory", required=True)
    parser.add_argument("-vc", "--vepcache", type=str, help="Directory for vep cache", required=True)
    parser.add_argument("-b", "--tumorbam", type=str, help="name of tumor bam, should not include full path", required=True)
    parser.add_argument("-n", "--normalbam", type=str, help="name of normal bam, should not include full path", default="")
    parser.add_argument("-nd", "--normaldir", type=str, help="directory that contains matched normal file",default="")
    parser.add_argument("-g","--genome",type=str,help="Genome build to use, default = GRCh37", default = "GRCh37")
    parser.add_argument("-f", "--fasta", type=str, help="path to fasta file", default="")
    parser.add_argument("-q","--mapq",type=int,help="minimum mapping quality, default = 10",default = 10)
    parser.add_argument("-Q","--baseq",type=int,help="minimum base quality, default = 10",default = 10)
    parser.add_argument("-s","--strand",type=int,help="Minimum number of reads mapping to forward and reverse strand to call mutation, default=2",default = 2)
    parser.add_argument("-m", "--mtchrom",type=str, help="Chromosome type", default="MT")
    parser.add_argument("-mo", "--molecule",type=str, help="Type of molecule (dna or rna), default=dna", default="dna")
    
    # read arguments
    args = parser.parse_args()
    datadir = args.datadir
    vcfdir = args.vcfdir
    outdir = args.outdir
    workingdir = args.workingdir
    vepcache = args.vepcache
    tumorbam = args.tumorbam
    normalbam = args.normalbam
    normaldir = args.normaldir
    genome = args.genome
    fasta = args.fasta
    minmapq = args.mapq
    minbq = args.baseq
    minstrand = args.strand
    mtchrom = args.mtchrom
    molecule = args.molecule
    
    # make output directories if they don't exist
    os.makedirs(vcfdir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True) 

    variant_calling(datadir, tumorbam, normalbam, normaldir, vcfdir, outdir, workingdir, vepcache, 
                    fasta, minmapq, minbq, minstrand, genome, mtchrom, molecule)
    
    print(f"DONE WITH MT VARIANT PIPELINE: {tumorbam}")

#!/usr/bin/env python3

import os
import argparse
import re
import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
import make_mutsig

# set parameters for genome build
def set_genome_params(genome, fasta, workingdir):
    valid_genomes = {"GRCm38": "/mm10/mm10_MT.fa", "mm10": "/mm10/mm10_MT.fa", "GRCh38": "/reference/GRCh38/genome_MT.fa", "GRCh37": "/reference/GRCh37/Homo_sapiens.GRCh37.dna.chromosome.MT.fa"}
    if genome not in valid_genomes.keys():
        raise ValueError(f"Genome {genome} is not supported")
    
    if(genome == "GRCm38"):
        genome = "mm10"
    
    ncbibuild = genome

    if fasta == "":
        fasta = workingdir + valid_genomes[genome]
    
    if(genome == "mm10"):
        species = "mus_musculus"
    else:
        species = "homo_sapiens"
    return genome, fasta, ncbibuild, species

def reference_detect(reffile):
    print("Determining the mtDNA chromosome name...")
    for sequence in SeqIO.parse(open(reffile), "fasta"):
        if re.search('MT', sequence.description.split(" ")[0]):
            return("MT")
        elif re.search('chrM', sequence.description.split(" ")[0]):
            return("chrM")
    raise Exception("Chromosome is neither MT nor chrM")

def merge_normal_tumor(resultsdir, tumor_id, normal_id):
    tumorfile = pd.read_csv(resultsdir + "/TEMPMAFfiles/tempMuTect2/" + tumor_id + ".bam.maf", sep = "\t", header=1, low_memory=False)
    normalfile = pd.read_csv(resultsdir + "/TEMPMAFfiles/tempMuTect2/" + normal_id + ".bam.maf", sep = "\t", header=1, low_memory=False)
    
    # keep only vars that are not in the normal 
    df_merged_outer = pd.merge(tumorfile[['Chromosome','Start_Position',
        'Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Variant_Type']], normalfile[['Chromosome','Start_Position',
        'Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Variant_Type']], how='outer', indicator=True)
    combinedfile = df_merged_outer[df_merged_outer['_merge'] == 'left_only'].drop(columns='_merge')
   
    # Combined matrices together
    combinedfile.to_csv(f"{resultsdir}/MuTect2_results/{tumor_id}.bam.maf",sep = '\t',na_rep='NA',index=False)


def variant_calling_normal(resultsdir,tumordir,tumor_id,reffile,genome,minmapq,minbq,minstrand,workingdir,vepcache,mtchrom,ncbibuild,species,normal_id,normaldir):
    os.makedirs(f"{resultsdir}/TEMPMAFfiles/tempMuTect2", exist_ok = True)
    os.makedirs(f"{resultsdir}/MuTect2_results", exist_ok = True)
    os.makedirs(f"{resultsdir}/MTvariant_results", exist_ok = True)

    # Running MTvariantpipeline with matched normal
    print("Running MTvariantpipeline with matched normal..")
    subprocess.run(f"python3 {workingdir}/MTvariantpipeline.py -d {tumordir}/ -v {resultsdir}/TEMPMAFfiles/ -w {workingdir}/ " +
        f"-o {resultsdir}/MTvariant_results/ -b {tumor_id}.bam -n {normal_id}.bam -nd {normaldir}/ -g {genome} -q {minmapq} " +
        f"-Q {minbq} -s {minstrand} -vc {vepcache} -f {reffile} -m {mtchrom} -mo {molecule} ", shell=True, check=True)

    # MuTect2 mitochondrial mode on tumor
    print("Running MuTect2 on tumor..")
    subprocess.run(f"gatk --java-options -Xmx4g Mutect2 -R {reffile} --mitochondria-mode true -L {mtchrom} " +
        f"-mbq {minbq} --minimum-mapping-quality {minmapq} --QUIET -I {tumordir}/{tumor_id}.bam " +
        f"-tumor {tumor_id.replace('-','_')} -O {resultsdir}/TEMPMAFfiles/tempMuTect2/{tumor_id}.bam.vcf.gz", shell=True, check=True)

    # MuTect2 mitochondrial mode on normal
    print("Running MuTect2 on normal..")
    subprocess.run(f"gatk --java-options -Xmx4g Mutect2 -R {reffile} --mitochondria-mode true -L {mtchrom} " +
        f"-mbq {minbq} --minimum-mapping-quality {minmapq} --QUIET -I {normaldir}/{normal_id}.bam " +
        f"-tumor {normal_id.replace('-','_')} -O {resultsdir}/TEMPMAFfiles/tempMuTect2/{normal_id}.bam.vcf.gz", shell=True, check=True)

    # Left align MuTect2 results (-m - is there for a reason)
    subprocess.run(f"bcftools norm -m - -f {reffile} {resultsdir}/TEMPMAFfiles/tempMuTect2/{tumor_id}.bam.vcf.gz " +
        f"-o {resultsdir}/TEMPMAFfiles/tempMuTect2/{tumor_id}.bam.vcf", shell=True, check=True)
    subprocess.run(f"bcftools norm -m - -f {reffile} {resultsdir}/TEMPMAFfiles/tempMuTect2/{normal_id}.bam.vcf.gz " +
        f"-o {resultsdir}/TEMPMAFfiles/tempMuTect2/{normal_id}.bam.vcf", shell=True, check=True)
    
    # Convert the MuTect2 result from vcf to maf file
    subprocess.run(f"perl {workingdir}/vcf2maf/vcf2maf.pl --species {species} --vep-data {vepcache} " +
        f"--ncbi-build {ncbibuild} --vep-overwrite --input-vcf {resultsdir}/TEMPMAFfiles/tempMuTect2/{tumor_id}.bam.vcf " +
        f"--output-maf {resultsdir}/TEMPMAFfiles/tempMuTect2/{tumor_id}.bam.maf --ref-fasta {reffile}", shell=True, check=True)
    subprocess.run(f"perl {workingdir}/vcf2maf/vcf2maf.pl --species {species} --vep-data {vepcache} " +
        f"--ncbi-build {ncbibuild} --vep-overwrite --input-vcf {resultsdir}/TEMPMAFfiles/tempMuTect2/{normal_id}.bam.vcf " + 
        f"--output-maf {resultsdir}/TEMPMAFfiles/tempMuTect2/{normal_id}.bam.maf --ref-fasta {reffile}", shell=True, check=True)
    
    merge_normal_tumor(resultsdir, tumor_id, normal_id)

def variant_calling(resultsdir,tumordir,tumor_id,reffile,genome,minmapq,minbq,minstrand,workingdir,vepcache,mtchrom,ncbibuild,species):

    os.makedirs(f"{resultsdir}/MuTect2_results", exist_ok = True)
    os.makedirs(f"{resultsdir}/MTvariant_results", exist_ok = True)
<<<<<<< HEAD
=======
    os.makedirs(f"{resultsdir}/TEMPMAFfiles/tempMuTect2", exist_ok = True)
>>>>>>> origin/master

    # Running MTvariantpipeline without matched normal
    print("Running MTvariantpipeline..")
    subprocess.run(f"python3 {workingdir}/MTvariantpipeline.py -d {tumordir}/ -v {resultsdir}/TEMPMAFfiles/ " +
        f"-o {resultsdir}/MTvariant_results/ -b {tumor_id}.bam -g {genome} -q {minmapq} -Q {minbq} -s {minstrand} " +
        f"-w {workingdir}/ -vc {vepcache} -f {reffile} -m {mtchrom} -mo {molecule} ", shell=True, check=True)

    # MuTect2 mitochondrial mode
    print("Running MuTect2..")
    subprocess.run(f"gatk --java-options -Xmx4g Mutect2 -R {reffile} --mitochondria-mode true -L {mtchrom} " +
        f"-mbq {minbq} --minimum-mapping-quality {minmapq} --QUIET -I {tumordir}/{tumor_id}.bam " +
        f"-tumor {tumor_id.replace('-','_')} -O {resultsdir}/MuTect2_results/{tumor_id}.bam.vcf.gz", shell=True, check=True)

    # Left align MuTect2 results (-m - is there for a reason)
    subprocess.run(f"bcftools norm -m - -f {reffile} {resultsdir}/MuTect2_results/{tumor_id}.bam.vcf.gz " +
        f"-o {resultsdir}/MuTect2_results/{tumor_id}.bam.vcf", shell=True, check=True)
    
    # Convert the MuTect2 result from vcf to maf file
    subprocess.run(f"perl {workingdir}/vcf2maf/vcf2maf.pl --species {species} --vep-data {vepcache} " +
        f"--ncbi-build {ncbibuild} --vep-overwrite --input-vcf {resultsdir}/MuTect2_results/{tumor_id}.bam.vcf " + 
<<<<<<< HEAD
        f"--output-maf {resultsdir}/MuTect2_results/{tumor_id}.bam.maf --ref-fasta {reffile}", shell=True, check=True)
=======
        f"--output-maf {resultsdir}/MuTect2_results/{tumor_id}_temp.bam.maf --ref-fasta {reffile}", shell=True, check=True)

    tumorfile = pd.read_csv(resultsdir + "/MuTect2_results/" + tumor_id + "_temp.bam.maf", sep = "\t", header = 1, low_memory = False)
    tumorfile = tumorfile[['Chromosome','Start_Position',
        'Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Variant_Type']]
    tumorfile.to_csv(f"{resultsdir}/MuTect2_results/{tumor_id}.bam.maf", sep = "\t", na_rep = "NA", index = False)
>>>>>>> origin/master

def variant_processing(tumor_id,resultsdir):
    """
    Run MTvariantpipeline and MuTect2 on the filtered cells
    MTvariantpipeline: A simple variant calling and annotation pipeline for mitochondrial DNA variants.
    """
    print("Starting variant processing...")

    # Overlap between MuTect and MTvariantpipeline
    MTvarfile = pd.read_csv(resultsdir + "/MTvariant_results/" + tumor_id + ".bam.maf", sep = "\t", comment='#', low_memory=False)
    mutectfile = pd.read_csv(resultsdir + "/MuTect2_results/" + tumor_id + ".bam.maf", sep = "\t", comment='#', header=0, low_memory=False)

    MTvarfile['Start_Position'] = MTvarfile['Start_Position'].astype(str)  
    mutectfile['Start_Position'] = mutectfile['Start_Position'].astype(str)


    # Filter out variants falling in the repeat regions of 513-525, and 3105-3109 (black listed regions)
    rmregions = list(range(513,524)) + list(range(3105,3109))
    if len(mutectfile['Start_Position'][mutectfile['Start_Position'].isin(rmregions)]) > 0:
        mutectfile = mutectfile[~mutectfile['Start_Position'].isin(rmregions)]
    if len(MTvarfile['Start_Position'][MTvarfile['Start_Position'].isin(rmregions)]) > 0:
        rmthese = MTvarfile['Start_Position'].isin(rmregions)
        MTvarfile = MTvarfile[~rmthese]
    mutectfile.index = range(len(mutectfile.index))
    MTvarfile.index = range(len(MTvarfile.index))
    
    # Output the overlap as final maf file
    filloutfile = pd.merge(mutectfile, MTvarfile, how='inner', on=['Chromosome','Start_Position',
        'Reference_Allele','Tumor_Seq_Allele2','Variant_Classification','Variant_Type'])
    
    filloutfile.index = [str(filloutfile['Reference_Allele'][i]) + ':' + str(int(filloutfile['Start_Position'][i])) + ':' + 
                         str(filloutfile['Tumor_Seq_Allele2'][i]) for i in range(len(filloutfile))]
    
    filloutfile.to_csv(f"{resultsdir}/{tumor_id}.bam.maf",sep = '\t',na_rep='',index=False)


if __name__ == "__main__":
    # Parse necessary arguments
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-t", "--tumor_id",type=str, help="(REQUIRED) Path of the tumor sample", required=True)
    parser.add_argument("-w", "--workingdir", type=str, help="(REQUIRED) Working directory with scripts", required=True)
    parser.add_argument("-re", "--resultsdir", type=str, help="(REQUIRED) Directory for results", required=True)
    parser.add_argument("-r", "--reffile",type=str, help="Reference fasta file",default="")
    parser.add_argument("-g", "--genome",type=str, help="Genome version, default=GRCh37",default = "GRCh37")    
    parser.add_argument("-q","--mapq",type=int,help="Minimum mapping quality, default = 20",default = 20)
    parser.add_argument("-Q","--baseq",type=int,help="Minimum base quality, default = 20",default = 20)
    parser.add_argument("-s","--strand",type=int,help="Minimum number of reads mapping to forward and reverse strand to call mutation, default=2",default = 2)
    parser.add_argument("-th","--threshold",type=int,help="The critical threshold for calling a cell wild-type, default=0.1",default = 0.1)
    parser.add_argument("-vc", "--vepcache", type=str, help="Directory for vep cache", default="$HOME/.vep")
    parser.add_argument("-n", "--normal_id", type=str, help="Path of the normal sample",default="")
    parser.add_argument("-mo", "--molecule",type=str, help="Type of molecule (dna or rna), default=dna", default="dna")

    # read in arguments
    args = parser.parse_args()
    tumor_id = args.tumor_id
    reffile = args.reffile
    genome = args.genome
    minmapq = args.mapq
    minbq = args.baseq
    minstrand = args.strand
    threshold = args.threshold
    workingdir = args.workingdir
    vepcache = args.vepcache
    resultsdir = args.resultsdir
    normal_id = args.normal_id
    molecule = args.molecule

    #print("Miminum mapping quality of " + str(minmapq))
    #print("Miminum base quality of " + str(minbq))
    #print("Miminum number of reads mapping to forward and reverse strand to call mutation of " + str(minstrand))
    #print("Reference file: " + reffile)

    genome, reffile, ncbibuild, species = set_genome_params(genome, reffile, workingdir)
    mtchrom = reference_detect(reffile)
    
    tumordir = os.path.dirname(tumor_id) if '/' in tumor_id else os.getcwd()
    tumor_id = os.path.splitext(os.path.basename(tumor_id))[0]

    if normal_id != "":
       normaldir = os.path.dirname(normal_id) if '/' in normal_id else os.getcwd()
       normal_id = os.path.splitext(os.path.basename(normal_id))[0]   
       variant_calling_normal(resultsdir,tumordir,tumor_id,reffile,genome,minmapq,minbq,minstrand,workingdir,vepcache,mtchrom,ncbibuild,species,normal_id,normaldir)
    else:
        variant_calling(resultsdir,tumordir,tumor_id,reffile,genome,minmapq,minbq,minstrand,workingdir,vepcache,mtchrom,ncbibuild,species)
    variant_processing(tumor_id,resultsdir)
    make_mutsig.create_mut_sig_file(resultsdir, tumor_id)
<<<<<<< HEAD
    print(f"DONE WITH BULKPIPELINE: {tumor_id}")
=======
    print(f"DONE WITH BULKPIPELINE: {tumor_id}")
>>>>>>> origin/master
# mtDNA-VC: Mitochondrial DNA Variant Calling Pipeline

This pipeline performs variant calling for mitochondrial DNA

## Installation (in the command line)

### Clone this repository:
```bash
git clone https://github.com/ziadbakouny18/mtdna-vc
``` 

## Install miniconda (at least Python v3.6):
https://docs.conda.io/en/latest/miniconda.html


Create a new conda environment and activate it (make sure it is activated before installing the dependencies):

```bash
conda env create -f mito.yml
conda activate mito
``` 

Verify installation of the mito environment
```bash
conda env list
``` 

## Setup VEP cache
Additionally, a VEP offline cache needs to be installed (NOTE: CACHE MUST BE SAME VERSION AS VEP VERSION and THE REFERENCE GENOME USED -- 112.0 in mito.yaml default). Please use https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/ to install vep cache. Due to the size of the caches, it will likely take several hours to install.

Then unzip
```bash
tar -xzvf homo_sapiens_vep_112_GRCh38.tar.gz
``` 
OR
```bash
tar -xzvf homo_sapiens_vep_112_GRCh37.tar.gz
``` 

## Choosing and downloading a reference genome
You will also need to ensure that the reference genome fasta has a dictionary. After downloading the correct reference genome (make sure the mitochondrial chromosome name matches the one in your bam file, i.e. chrM or MT).
You can use the reference genome that your reads were mapped to but make sure the mitochondrial genome you use uses rCRS as described here:
https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use

For example:
```bash
curl -O ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
tar -xzvf human_g1k_v37.fasta.gz
``` 

To make "genome.dict" file run the following in the same directory as your reference genome.fa file:
```bash
picard CreateSequenceDictionary \
    R=human_g1k_v37.fasta \
    O=human_g1k_v37.dict
``` 

## Prepare BAM by subsetting to mitochondrial reads
Subsetting for mitochondrial regions (where chrM can either be chrM or MT, depending on how it is annotated in the input BAM file):
```bash
conda activate mito
samtools view -b [library_id].bam chrM > [library_id]_chrM.bam
``` 

Index the input .bam file:
```bash
samtools index [library_id]_chrM.bam
``` 

# Running the pipeline
python3 "$script_dir/bulkpipeline.py" \
  -t "$bam_path" \                     # Path to input BAM or file list
  -w "$script_dir" \                   # Path to pipeline scripts
  -re "$results_dir/$sample_name" \   # Output directory for results
  -vc "$vep_cache_dir" \              # VEP cache directory
  -r "$ref_file" \                    # Reference genome FASTA
  -g "$genome"                        # Genome build (e.g., GRCh37 or GRCh38)

# Usage example
python3 "$script_dir/bulkpipeline.py" \
  -t /data1/reznike/elbakoz/data/tang_cellmetabolism_2023/example_mt.bam \
  -w /home/elbakoz/mtdna-vc/python \
  -re /data1/reznike/elbakoz/results/tang_cellmetabolism_2023/mitochondrial_variant_calling \
  -vc /data1/reznike/elbakoz/reference/grch37/human_g1k_v37.fasta \
  -r /data1/reznike/elbakoz/reference/vep_cache \
  -g "GRCh37"

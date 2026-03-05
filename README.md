# Project Name

HapBridge: A Methylation-Guided Approach for Correcting Switch Errors and Bridging Phased Blocks in Long-Read Phasing
## Introduction

> We propose HapBridge, which corrects and bridges SNP-based phased blocks by selecting haplotype-specific methylation and utilizing the consistency of methylation signals.

## Installation and Usage

### Prerequisites
Ensure you have the following tools installed:
```sh
- Python >= 3.8
- pysam >= 0.14
- samtools >= 1.18
- HapCUT2 = v1.3.4
- Whatshap = 1.7
- clair3 docker version >= v0.1-r12
- BGZIP (bgzip tabix)
```

### Installation

```sh
# Clone the repository
$ git clone https://github.com/Humonex/HapBridge.git

# Navigate to the project directory
$ cd HapBridge
$ git clone https://github.com/vibansal/HapCUT2.git

### Usage

```sh 
# Run the main python.py
```


## Usage Examples


```python
# Import the library
conda create -n HapBridge python=3.9
conda activate HapBridge
conda install pysam
conda install whatshap
conda config --add channels conda-forge
conda config --add channels bioconda
conda install whatshap
whatshap --version (1.7)

# Use a specific feature
python bridge.py phased.vcf.gz haplotagged.bam HapBridge_phased.vcf -t 30

# Required arguments:
phased.vcf.gz --The SNP-based phasing vcf file (Hapcut2 or Whatshap)
haplotagged.bam --The read tagged by whatshap haplotag
HapBridge_phased.vcf --The phasing result by HapBridge
-t  --thread --HapBridge supports multi-threaded operation
```

## A complete experimental process


```sh
# Call SNP

The SNP called by Clair3 https://github.com/HKU-BAL/Clair3

# SNP Phasing (whatshap)

whatshap phase -o whatshap_phased.vcf --reference=reference.fasta input.vcf input.bam --ignore-read-groups
bgzip phased.vcf
tabix -p vcf phased.vcf.gz

# or SNP Phasing (HapCUT2)
cd HapBridge/HapCUT2
./build/extractHAIRS  --ont 1 --bam input.bam --VCF input.vcf --out fragment_file --ref reference.fasta
./build/HAPCUT2 --fragments fragment_file --VCF input.vcf --output haplotype_output_file_hapcut2

# Whatshap haplotag: Tagging reads by haplotype

whatshap haplotag -o haplotagged.bam --reference reference.fasta haplotype_output_file_hapcut2.vcf.gz(or whatshap_phased.vcf.gz) input.bam
samtools index haplotagged.bam

# HapBridge
python bridge.py haplotype_output_file_hapcut2.vcf.gz haplotagged.bam HapBridge_phased.vcf -t 30

```
## NOTE

Note:
1. The Original Bam file must contain MM, ML (methylation) tags
2. Please check whether the bam file after Whatshap haplotag has PS and HP tags
3. if pysam has MM tag reading error, please try this command: samtools view -bF 2304 -o output.bam input.bam

## License

This project is distributed under the [GNU General Public License v3.0](LICENSE).
Preprints can be obtained from the following link：https://doi.org/10.1101/2025.11.07.687303


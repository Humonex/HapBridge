# Project Name

HapBridge: Improving Phasing Performance and Correcting Switch Errors Using Methylation in Long Reads

## Introduction

> We propose HapBridge, which corrects and bridges SNP-based phased blocks by selecting haplotype-specific methylation and utilizing the consistency of methylation signals.

## Installation and Usage

### Prerequisites

```sh
- Python >= 3.8
- pysam >= 0.14
- samtools >= 1.18
- Whatshap >= 1.7
- clair3 docker version >= v0.1-r12
```

### Installation

```sh
# Clone the repository
$ git clone https://github.com/Humonex/HapBridge.git

# Navigate to the project directory
$ cd HapBridge

### Usage

```sh 
# Run the main python.py
$ python bridge.py
```


## Usage Examples


```python
# Import the library
conda create -n HapBridge python=3.9
conda install pysam

# Use a specific feature
python bridge.py phased.vcf.gz haplotagged.bam HapBridge_phased.vcf -t 30

# Required arguments:
phased.vcf.gz --The SNP-based phasing vcf file (Hapcut2 or Whatshap)
haplotagged.bam --The read tagged by whatshap harlot
HapBridge_phased.vcf --The phasing result bt HapBridge
-t  --thread --HapBridge supports multi-threaded operation
```

## A complete experimental process


```sh
# Call SNP

The SNP called by Clair3 https://github.com/HKU-BAL/Clair3

# SNP Phasing

whatshap phase -o phased.vcf --reference=reference.fasta input.vcf input.bam --ignore-read-groups
bgzip phased.vcf
tabix -p vcf phased.vcf.gz

# Whatshap haplotag: Tagging reads by haplotype

whatshap haplotag -o haplotagged.bam --reference reference.fasta phased.vcf.gz alignments.bam
samtools index haplotagged.bam

# HapBridge
python bridge.py phased.vcf.gz haplotagged.bam HapBridge_phased.vcf -t 30

```
## NOTE

Note:
1. The Original Bam file must contain MM, ML (methylation) tags
2. Please check whether the bam file after Whatshap haplotag has PS and HP tags
3. if pysam has MM tag reading error, please try this command: samtools view -bF 2304 -o output.bam input.bam

## License

This project is distributed under the [GNU General Public License v3.0](LICENSE).


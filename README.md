# Project Name

HapBridge: impoving phase performance using methylation in long reads

## Introduction

> We propose HapBridge, which corrects and bridges SNP-based phased blocks by selecting haplotype-specific methylation and utilizing the consistency of methylation signals.

## Installation and Usage

### Prerequisites

```sh
- Python >= 3.8
- pysam >= 0.14
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

phased.vcf.gz --The SNP-based phasing vcf file (Hapcut2 or Whatshap)
haplotagged.bam --The read tagged by whatshap harlot
HapBridge_phased.vcf --The phasing result bt HapBridge
-t  --thread --HapBridge supports multi-threaded operation


## License

This project is distributed under the [GNU General Public License v3.0](LICENSE).

---

Thank you for your support! If you have any questions or suggestions, please reach out via [issues](https://github.com/username/project_name/issues).

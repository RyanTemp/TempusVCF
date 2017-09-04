# TempusVCF

## Technical Challenge Description
For this challenge, you are asked to prototype a variant annotation tool. We will provide you with a VCF file, and you will create a small software program to output a table annotating each variant in the file. Each variant must be annotated with the following pieces of information:
1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.) If there are multiple possibilities, annotate with the most deleterious possibility.
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from Broad Institute ExAC Project API
(API documentation is available here: http://exac.hms.harvard.edu/)
6. Additional optional information from ExAC that you feel might be relevant.
For this project please upload all relevant code (written in whatever language you like) along with the annotated VCF file to a Github account and provide the link to the below email address. Please note that work will be assessed based on quality of code and documentation more-so than the annotation.

## Tool Description
Here, I developed a small tool in response to the above described challenge. This tool will take in the input vcf file, then using SnpEff, it will annotate the deleteriousness and impact of each variants. 

## Dependencies
1. This tool is written in python2.7. Please first make sure python is installed.

2. This tool requires the latest version of SnpEff.
Download latest version of SnpEff and install it under `HOME` directory
```bash
cd ~
wget https://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
```
Download the `GRCh37.75` database.
```bash
java -jar $HOME/snpEff/snpEff.jar download GRCh37.75
```

3. In addition to a few general packages `json` `urllib2` `os` `sys` `getopt` `time`, this tool also depends on PyVCF to parse and write vcf files.
```bash
pip install pyvcf
```
## Installation
1. Download `annotation.py` script.
2. In the `main()` function, modify `SNPEFF_DIR` variable to indicate the path to the directory containing the `snpEff.jar` file. Usually, it is the unzipped directory from the downloaded snpEff_latest_core.zip file. As recommended above, this directory will be `$HOME/snpEff/`.
3. In the directory containing `annotation.py` script, run
```bash
python annotation.py -i Challenge_data.vcf -o Challenge_data.ann.vcf
```

# Ecoli_wgs_typing

# Streptococcus suis whole genome assembly and typing pipeline


Pipeline for whole genome assembly and analysis of E.coli. Works only for Oxford Nanopore reads



### Usage
Requires input directory containg sub-directories with the fastq files and output directory. Outputs several intermediate files with a html report with AMR,MLST,Serotyping and virulence factors found in the sample.
```
nextflow run main.nf --input dirwith_subdir_ --outdir Results -profile docker --trim_barcodes
```
```
Parameters:

--input		Input directory containg sub-sirectories with fastq files
--out_dir	Output directory
optional
--trim_barcodes barcode and adapter trimming using porechop
```
### Dependencies
* nextflow
* docker
* wsl2
### Software and references used
* dragonflye (https://github.com/rpetit3/dragonflye)
* porechop (https://github.com/rrwick/Porechop)
* minimap2 (Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191)
* abricate (https://github.com/tseemann/abricate)
* mlst (https://github.com/tseemann/mlst,This publication made use of the PubMLST website)
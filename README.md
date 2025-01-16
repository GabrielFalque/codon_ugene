# codon_ugene

## Description
An R script for analysing codon usage. It calculates, for each gene, the Relative Synonymous Codon Usage (RSCU), the Codon Adaptation Index (CAI), the Effective Number of Codons (ENC), the GC content, the values of GC1, GC2, GC3, GC12, the Relative Codon Deoptimization Index (RCDI), the frequency of dinucleotides and their Relative Dinucleotide Abundance, the proportion of aromatic amino acids, and the GRAVY scrore.

It also performs a Correspondence Analysis (CA) on the RSCU values.

## Installation et Preparation
1. **Clone this repo** :
```bash
git clone https://github.com/GabrielFalque/codon_ugene.git
cd codon_ugene
```
   
2. **Install dependencies** :
```bash
Rscript setup.R
```

3. **Execute script** :
```bash
Rscript codon_ugene.R
```

## Usage
```bash
Rscript ./codon_ugene.R --dir genes/ --outdir outdir/


Options:
	-d DIR, --dir=DIR
		Directory with gene sequences (mandatory).
Files must be named like gene_{gene_name}.fasta with gene_name being Orf1 for example. (Ex : gene_Orf1.fasta

	-o OUTDIR, --outdir=OUTDIR
		Directory where to store output files (mandatory).

	-r FILE, --ref=FILE
		Species codon usage CSV (by default :'human_codon_usage.csv'). 
You can download RSCU table from the desired species directly from http://codonstatsdb.unr.edu/index.html.

	--cov
		If the sequences are from SARS-CoV-2.

	--verbose
		Verbose mode.

	-h, --help
		Show this help message and exit

```

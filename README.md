# codon_ugene

## Description
An R script for analysing codon usage. It calculates, for each gene, the Relative Synonymous Codon Usage (**RSCU**), the Codon Adaptation Index (**CAI**), the Effective Number of Codons (**ENC**), the **GC content**, the values of **GC1**, **GC2**, **GC3**, **GC12**, the Relative Codon Deoptimization Index (**RCDI**), the **dinucleotides frequency** and their **Relative Dinucleotide Abundance**, the **proportion of aromatic amino acids**, and the **GRAVY score**.

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
Sans précision, ce script comparera les différentes données obtenues à celles du génome humain. Si vous voulez le comparer au génome d'une autre espèce, il faut fournir en argument `--ref my_species.csv`. Le fichier *csv* en question doit contenir les informations sur l'usage des codons et doit être dans le même format que le fichier *codon_ugene/data/human_codon_usage.csv*. Les tables d'usage de codons de nombreuses espèces, dans le bon format, peuvent être trouvées sur __[Codon Statistics Database](http://codonstatsdb.unr.edu/index.html)__ (Subramanian et al. 2022 [[1]](#1)).

## References
<a id="1">[1]</a> 
[Subramanian et al. 2022](https://academic.oup.com/mbe/article/39/8/msac157/6647594?login=false). 


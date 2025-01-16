# codon_ugene

## Description
An R script for analysing codon usage. It calculates, for each gene, the Relative Synonymous Codon Usage (**RSCU**), the Codon Adaptation Index (**CAI**), the Effective Number of Codons (**ENC**), the **GC content**, the values of **GC1**, **GC2**, **GC3**, **GC12**, the Relative Codon Deoptimization Index (**RCDI**), the **dinucleotides frequency** and their **Relative Dinucleotide Abundance**, the **proportion of aromatic amino acids**, and the **GRAVY score**.

It also performs a Correspondence Analysis (CA) on the RSCU values.

## Installation and Preparation
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
The directory containing the gene sequences must therefore contain one *.fasta* file per gene, the name of which must be *gene_{gene_name}.fasta* with gene_name being Orf1 for example. (Ex: *gene_Orf1.fasta*). The fasta file must contain several sequences of the same gene. The more sequences there are, the more accurate the results will be. Sequences with too many undefined nucleotides (‘n’) will not be taken into account (more than 15% ‘n’). In addition, sequences for which anomalies are detected (no START, STOP, internal STOP codon, number of nucleotides not a multiple of 3 in the coding part) will be removed from the analysis directly by the script. If the sequences come from an alignment, the gaps will be removed.

Please keep the same sequence name (in the header) between the different gene FASTA files.

Without precision, this script will compare the different data obtained to that of the human genome. If you want to compare it to the genome of another species, you need to supply `--ref my_species.csv` as an argument. The *csv* file in question must contain information on codon usage and must be in the same format as the *codon_ugene/data/human_codon_usage.csv* file. Codon usage tables for many species, in the correct format, can be found at __[Codon Statistics Database](http://codonstatsdb.unr.edu/index.html)__ (Subramanian et al. 2022 [[1]](#1)).

### Gene sequences preparation

Personally, I use complete virus genome sequences. I choose a reference genome (one with the fewest missing parts and ambiguous nucleotides), and retrieve its annotated gene file from __[NCBI](https://www.ncbi.nlm.nih.gov/)__ [[2]](#2)). I align my sequences (via __[MAFFT](https://mafft.cbrc.jp/alignment/server/index.html)__ [[3]](#3)) for example). Then I separate my alignment by gene. To do this I use a *python* script (which I'll put on my github soon) such as :
```bash
python ./extract_genes_aln.py --ref_file virus_annotated_genes.txt --aligned_file complete_genome_sequences.mafft.fasta --seq_id seq1 --output_directory my/results/
```
with `--seq_id seq1` specifying the name of the sequence in my fasta file from which I took the annotation file *virus_annotated_genes.txt* from NCBI.

## Outputs



## References
<a id="1">[1]</a> 
Krishnamurthy Subramanian, Bryan Payne, Felix Feyertag, David Alvarez-Ponce, The Codon Statistics Database: A Database of Codon Usage Bias, Molecular Biology and Evolution, Volume 39, Issue 8, August 2022, msac157, [https://doi.org/10.1093/molbev/msac157](https://doi.org/10.1093/molbev/msac157). 

<a id="2">[2]</a> 
Benson DA, Cavanaugh M, Clark K, Karsch-Mizrachi I, Lipman DJ, Ostell J, Sayers EW. GenBank. Nucleic Acids Res. 2013 Jan;41(Database issue):D36-42. doi: 10.1093/nar/gks1195. Epub 2012 Nov 27. PMID: 23193287; PMCID: PMC3531190, [https://pubmed.ncbi.nlm.nih.gov/23193287/](https://pubmed.ncbi.nlm.nih.gov/23193287/). 

<a id="3">[3]</a> 
Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol Biol Evol. 2013 Apr;30(4):772-80. doi: 10.1093/molbev/mst010. Epub 2013 Jan 16. PMID: 23329690; PMCID: PMC3603318, [https://pubmed.ncbi.nlm.nih.gov/23329690/](https://pubmed.ncbi.nlm.nih.gov/23329690/). 

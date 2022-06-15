# DNA Methylase Finder
DNA Methylase Finder is a tool to detect/predict DNA methylase genes of prokaryotes. More specifically, it will detect genes with DNA methylase DOMAINS by querying their amino acid sequences. This tool was able to correctly call all "gold standard" DNA methlyases in the REBASE [database](http://rebase.neb.com/cgi-bin/rebgoldlist) (proteins with `M.` prefix). De novo discovery of methylases was conducted in hundreds of bacterial genomes and metagenomic assemblies. I have not noticed any systematic false positive patterns, but false negative rate would be difficult to assess.

## Input
Nucleotide contigs/genomes (`.fna`) or amino acids (`.faa`).

## Output (important files)
**All input types:**
`.DNA_methylases.combined.summary.tsv` : summary table of all detected DNA methylases, their coordinates, subtype, and specificity. Open in text editor or Excel, etc.
`.DNA_methylases.combined.faa` : multi-fasta of amino acid sequences of predicted DNA methylase genes (including merged fragments)

**Nucleotide inputs only:**
`neighborhood_annotations/*neighbors.gb` : DNA methylase gene neighborhood map of DNA methylase and surrounding genes. Good for determining if each DNA methylase gene is part of a Restriction Modification System. Open in GenBank file/plasmid viewer, e.g. SnapGene Viewer or UGENE.


## Installation instructions

NOTE: I have only tested this on a Linux system. I imagine it will only work on Linux.

1. Make sure you have conda installed

`conda -V`

2. Clone this repo

`git clone !`

3. Change to the `DNA_methylase_finder` directory:

`cd DNA_methylase_finder`


4. Create the conda environment called `dna_methylase_finder`:

`conda env create --file dna_methylase_finder.yml`

5. Download and unpack the databases (3.5 Gb compressed, 11 Gb decompressed):

```
wget ! 
tar -xvf DNA_methylase_finder_DBS_v1.0.tar.gz
```

## Usage


### Tips
* Setting `--neighborhoods False` will reduce the runtime by quite a bit and is suggested for finding DNA methylase genes in large datasets
* Increasing number of CPUs available with `-t` will make the tool run faster.
* See something weird? Open an Issue on this repo describing your problem in detail.


## Help Menu
```
usage: run_DNA_methylase_finder.py [-h] -it INPUT_TYPE 
                                        -f INPUT_FILE 
                                        -r RUN_TITLE 
                                        [-t CPU]
                                        [--meth_hmms METHYLASE_HMMS] [--cdd_plus_hmms CDD_PLUS_HMMS]
                                        [--legit_domains LEGIT_DOMAIN_LIST] [--motif_blastp MOTIF_ANNOTATE_BLASTP]
                                        [--subtype_hmms SUBTYPE_ANNOTATE_HMM] [--prod_args PROD_ARGS] [--pid PID]
                                        [--cov COV] [--s_subunit_hmms S_SUBUNIT_HMM] [--re_hmms RE_HMM]
                                        [--neighborhoods NEIGHBORHOODS] [--merge MERGE]

DNA Methylase Finder, v1.0

optional arguments:
  -h, --help            show this help message and exit

 REQUIRED ARGUMENTS for DNA Methylase Finder, v1.0 :
  -it INPUT_TYPE, --input_type INPUT_TYPE
                        OPTIONS: nucl, AA -- nucl PREFERRED! nucl is a nucleotide fasta file .fna extension. Each
                        header must be unique before the first space character. AA is an amino acid fasta file with
                        fasta file .fna extension. Each header must be unique before the first space character.
  -f INPUT_FILE, --input_file INPUT_FILE
                        nucl file with .fna extenstion, prodigal directory, or AA seq file with .faa extension.
  -r RUN_TITLE, --run_title RUN_TITLE
                        Name of this run. A directory of this name will be created. Must be unique from older runs or
                        older run will be renamed. Must consist of ONLY letters, numbers and underscores (_)

 OPTIONAL ARGUMENTS for DNA Methylase Finder, v1.0:
  -t CPU, --cpu CPU     Default: 4 -- Number of CPUs available for run.
  --meth_hmms METHYLASE_HMMS
                        Default: standard database -- Hmmer-formatted file of HMMs of putative DNA methylases
  --cdd_plus_hmms CDD_PLUS_HMMS
                        Default: standard database -- Hmmer-formatted file of HMMs of all CDD + putative DNA methylases
  --legit_domains LEGIT_DOMAIN_LIST
                        Default: standard list -- text file (1 entry per line) with names of DNA methyalse Hmmer models
  --motif_blastp MOTIF_ANNOTATE_BLASTP
                        Default: standard database -- BLASTP-formatted file of all REBASE DNA methylase proteins with
                        motif tag
  --subtype_hmms SUBTYPE_ANNOTATE_HMM
                        Default: standard database -- Hmmer-formatted file of HMMs of Olveira subtype HMMs with subtype
                        tag@@
  --prod_args PROD_ARGS
                        Default: -c -p meta -- arguments for prodigal in quotation marks. Only relvant for --input_type
                        nucl (-it nucl). Make sure to keep settings to produce AA, nucleotide, and gtf files from
                        prodigal step. Do not use memory or CPU arguments.
  --pid PID             Default: 80 -- minimum threshold for AA Percent Identity of predicted methylase gene to REBASE
                        homolog to predict motif specificity
  --cov COV             Default: 80 -- minimum threshold for alignment coverage of predicted methylase gene to REBASE
                        homolog to predict motif specificity
  --s_subunit_hmms S_SUBUNIT_HMM
                        Default: standard database -- Hmmer-formatted file of HMMs specificity subunit proteins
  --re_hmms RE_HMM      Default: standard database -- Hmmer-formatted file of HMMs restriction enzyme (endonuclease)
                        proteins
  --neighborhoods NEIGHBORHOODS
                        Default: True -- Make DNA methylase gene neighborhood maps? True - OR - False
  --merge MERGE         Default: True -- Merge adjacent DNA methylases with the assumption that they are a broken ORF?
                        True - OR - False
```
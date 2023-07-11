#!/usr/bin/env python

import argparse
import sys, os
import subprocess
import importlib.util

__version__='1.0.1'

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
        
pathname = os.path.dirname(sys.argv[0])  
find_methylase_script_path = os.path.abspath(pathname)      
print(find_methylase_script_path) 

parser = argparse.ArgumentParser(description='DNA Methylase Finder v' + str(__version__))

required_args = parser.add_argument_group(' REQUIRED ARGUMENTS for DNA Methylase Finder, v' + str(__version__))

required_args.add_argument("-it", "--input_type", dest="INPUT_TYPE", type=str, required=True, help='OPTIONS: nucl, AA -- nucl PREFERRED! nucl is a nucleotide fasta file .fna extension. Each header must be unique before the first space character. AA is an amino acid fasta file with fasta file .fna extension. Each header must be unique before the first space character.')
required_args.add_argument("-f", "--input_file", dest="INPUT_FILE", type=str, required=True, help='nucl file with .fna extenstion, prodigal directory, or AA seq file with .faa extension.')
required_args.add_argument("-r", "--run_title", dest="run_title", type=str, required=True, help='Name of this run. A directory of this name will be created. Must be unique from older runs or older run will be renamed. Must consist of ONLY letters, numbers and underscores (_)')

#-#-#
optional_args = parser.add_argument_group(' OPTIONAL ARGUMENTS for DNA Methylase Finder, v' + str(__version__))

optional_args.add_argument('--version', action='version', version=str(__version__))
optional_args.add_argument("-t", "--cpu", dest="CPU", type=int, default=4, help='Default: 4 -- Number of CPUs available for run. ')
optional_args.add_argument("--meth_hmms", dest="METHYLASE_HMMs", type=str, default=str(find_methylase_script_path) + '/methylase_hmms/meth_hmms_v1.0', help='Default: standard database -- Hmmer-formatted file of HMMs of putative DNA methylases')
optional_args.add_argument("--cdd_plus_hmms", dest="CDD_PLUS_HMMs", type=str, default=str(find_methylase_script_path) + '/cdd_plus_hmms/cdd_plus_hmms_v1.0', help='Default: standard database -- Hmmer-formatted file of HMMs of all CDD + putative DNA methylases')
optional_args.add_argument("--legit_domains", dest="LEGIT_DOMAIN_LIST", type=str, default=str(find_methylase_script_path) + '/legit_DNA_methylase_domain_model_list_v1.0.txt', help='Default: standard list -- text file (1 entry per line) with names of DNA methyalse Hmmer models ')
optional_args.add_argument("--motif_blastp", dest="MOTIF_ANNOTATE_BLASTP", type=str, default=str(find_methylase_script_path) + '/motif_protein_blastp/motif_blastp_v1.0', help='Default: standard database -- BLASTP-formatted file of all REBASE DNA methylase proteins with motif tag')
optional_args.add_argument("--subtype_hmms", dest="SUBTYPE_ANNOTATE_HMM", type=str, default=str(find_methylase_script_path) + '/subtype_hmms/subtype_hmms_olveira_v1.0', help='Default: standard database -- Hmmer-formatted file of HMMs of Olveira subtype HMMs with subtype tag@@')
optional_args.add_argument("--prod_args", dest="PROD_ARGS", type=str, default="-c -p meta", help='Default: -c -p meta -- arguments for prodigal in quotation marks. Only relvant for --input_type nucl (-it nucl). Make sure to keep settings to produce AA, nucleotide, and gtf files from prodigal step. Do not use memory or CPU arguments.')
optional_args.add_argument("--pid", dest="PID", type=str, default="80", help='Default: 80 -- minimum threshold for AA Percent Identity of predicted methylase gene to REBASE homolog to predict motif specificity')
optional_args.add_argument("--cov", dest="COV", type=str, default="80", help='Default: 80 -- minimum threshold for alignment coverage of predicted methylase gene to REBASE homolog to predict motif specificity')
optional_args.add_argument("--s_subunit_hmms", dest="S_SUBUNIT_HMM", type=str, default=str(find_methylase_script_path) + '/specificity_subunit_hmms/specificity_subunit_hmms_v1.0', help='Default: standard database -- Hmmer-formatted file of HMMs specificity subunit proteins')
optional_args.add_argument("--re_hmms", dest="RE_HMM", type=str, default=str(find_methylase_script_path) + '/restriction_enzyme_hmms/RE_hmms_v1.0', help='Default: standard database -- Hmmer-formatted file of HMMs restriction enzyme (endonuclease) proteins')
optional_args.add_argument("--neighborhoods", dest="NEIGHBORHOODS", type=str2bool, default="True", help='Default: True -- Make DNA methylase gene neighborhood maps? True - OR - False')
optional_args.add_argument("--merge", dest="MERGE", type=str2bool, default="True", help='Default: True -- Merge adjacent DNA methylases with the assumption that they are a broken ORF? True - OR - False')


args = parser.parse_args()

def is_tool(name):
	"""Check whether `name` is on PATH."""
	from distutils.spawn import find_executable
	return find_executable(name) is not None

if is_tool("bioawk") :
	print ("bioawk found")
else:
	print ("bioawk is not found. Exiting.")
	quit()

if is_tool("hmmscan") :
	print ("hmmscan found")
else:
	print ("hmmscan is not found. Exiting.")
	quit()
if is_tool("blastp") :
	print ("blastp found")
else:
	print ("blastp is not found. Exiting.")
	quit()
if is_tool("prodigal") :
	print ("prodigal found")
else:
	print ("prodigal is not found. Exiting.")
	quit()

package_name = 'Bio'
spec = importlib.util.find_spec(package_name)
if spec is None:
    print(package_name +" is not installed")
    quit()
else:
	print ("Biopython found")

package_name = 'BCBio'
spec = importlib.util.find_spec(package_name)
if spec is None:
    print(package_name +" is not installed")
    quit()
else:
	print ("BCBio found")


package_name = 'numpy'
spec = importlib.util.find_spec(package_name)
if spec is None:
    print(package_name +" is not installed")
    quit()
else:
	print ("numpy found")

subprocess.call(['bash', str(find_methylase_script_path) + '/DNA_methylase_finder_v1.0.sh', str(args.INPUT_TYPE), str(args.INPUT_FILE), str(args.run_title), str(args.CPU), str(args.METHYLASE_HMMs), str(args.CDD_PLUS_HMMs), str(args.LEGIT_DOMAIN_LIST), str(args.MOTIF_ANNOTATE_BLASTP), str(args.SUBTYPE_ANNOTATE_HMM), str(args.PROD_ARGS), str(args.PID), str(args.COV), str(args.S_SUBUNIT_HMM), str(args.RE_HMM), str(find_methylase_script_path), str(args.NEIGHBORHOODS), str(args.MERGE)])


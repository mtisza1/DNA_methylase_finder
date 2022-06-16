#!/bin/bash

# script to thoroughly and accurately identify all DNA methylase genes from a protein fasta file ".faa". Headers in ".faa" file must be unique from each other before the first space character.

echo "DNA Methylase Finder!"
echo "Version 1.0"

INPUT_TYPE=$1
INPUT_FILE=$2
OUTPUT_DIRECTORY=$3
CPU=$4
METHYLASE_HMMs=$5
CDD_PLUS_HMMs=$6
LEGIT_DOMAIN_LIST=$7
MOTIF_ANNOTATE_BLASTP=$8
SUBTYPE_ANNOTATE_HMM=$9
PROD_ARGS=${10}
PID=${11}
COV=${12}
S_SUBUNIT_HMM=${13}
RE_HMM=${14}
FIND_METHYLASE_DIR=${15}
NEIGHBORHOODS=${16}
MERGE=${17}
BASE_DIRECTORY=$PWD


# Making output folder
if [ ! -d "$OUTPUT_DIRECTORY" ]; then
	mkdir "$OUTPUT_DIRECTORY"
else
	rand_dir=$( head /dev/urandom | tr -dc A-Za-z0-9 | head -c 3 ; echo '' )
	DAY1=$( date +"%m-%d-%y" )
	mv ${OUTPUT_DIRECTORY}/ ${OUTPUT_DIRECTORY}_old_${DAY1}_${rand_dir} 
	mkdir "$OUTPUT_DIRECTORY"
fi 

if [ ${INPUT_TYPE} == "nucl" ] ; then
	INPUT_NUCL=${INPUT_FILE}
	if [ ${INPUT_NUCL: -4} == ".fna" ] && [ -s $INPUT_NUCL ]; then
		echo "A seemingly good .fna file found. Continuing"
	else
		echo "A .fna file not provided. Exiting."
		exit
	fi
	if [ -s ${BASE_DIRECTORY}/${INPUT_NUCL} ] ; then 
		echo ${BASE_DIRECTORY}/${INPUT_NUCL} ; 
	else  
		cp ${INPUT_NUCL} ${BASE_DIRECTORY}/ ; 
		INPUT_NUCL=$( basename $INPUT_NUCL ) 
	fi
	cd $OUTPUT_DIRECTORY
	## run prodigal
	MDYT=$( date +"%m-%d-%y---%T" )

	echo "running prodigal on nucleotide file -- $MDYT"
	#-# edit the input file argument to specify absolute path, so users can invoke a .fna file from anywhere.
	prodigal -i ../${INPUT_NUCL} -a ${INPUT_NUCL%.fna}.faa -d ${INPUT_NUCL%.fna}.genes.fna -f gff -o ${INPUT_NUCL%.fna}.out -q ${PROD_ARGS} 
	INPUT_AA=${INPUT_NUCL%.fna}.faa
elif [ ${INPUT_TYPE} == "AA" ] ; then
	INPUT_AA=${INPUT_FILE}
	if [ ${INPUT_AA: -4} == ".faa" ] && [ -s $INPUT_AA ]; then
		echo "A seemingly good .faa file found. Continuing"
	else
		echo "A .faa file not provided. Exiting."
		exit
	fi

	if [ -s ${BASE_DIRECTORY}/${INPUT_AA} ] ; then 
		cp ${BASE_DIRECTORY}/${INPUT_AA} ${BASE_DIRECTORY}/${OUTPUT_DIRECTORY}/  
	else  
		cp ${INPUT_AA} ${BASE_DIRECTORY}/${OUTPUT_DIRECTORY}/ ; 
		INPUT_AA=$( basename $INPUT_AA ) 
	fi
	cd $OUTPUT_DIRECTORY
else
	echo "Correct argument for --input_type not given. Use nucl, or AA. Exiting."
	exit
fi

MDYT=$( date +"%m-%d-%y---%T" )
echo "Step 1. Splitting AA file for parallelization -- $MDYT"

bioawk -c fastx '{print ">"$name"_" ; print $seq}' ${INPUT_AA} > ${INPUT_AA%.faa}.sort.faa

TOTAL_AA_SEQS=$( grep -F ">" ${INPUT_AA%.faa}.sort.faa | wc -l | bc )
if [ $TOTAL_AA_SEQS -ge 1 ] ; then 
	AA_SEQS_PER_FILE=$( echo "scale=0 ; $TOTAL_AA_SEQS / $CPU" | bc )
	if [ $AA_SEQS_PER_FILE = 0 ] ; then
		AA_SEQS_PER_FILE=1
	fi
	awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_INPUT_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < ${INPUT_AA%.faa}.sort.faa
	SPLIT_AA_LARGE=$( find * -maxdepth 0 -type f -name "SPLIT_INPUT_AA_*.fasta" )
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "Step 2. Querying input AA seqs against DNA methylase HMMs -- $MDYT"
	echo "$SPLIT_AA_LARGE" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --domtblout {}.hmmscan.out --cpu 1 -E 1e-4 --noali $METHYLASE_HMMs {}.fasta >/dev/null 2>&1
	cat SPLIT_INPUT_AA_*.hmmscan.out | grep -v "^#" | sed 's/ \+/	/g' | sort -u -k4,4 > INPUT_AA_COMBINED.hmmscan.sort.out
	if [ ! -z INPUT_AA_COMBINED.hmmscan.sort.out ] ; then
		echo "found some preimlinary hits"
		cut -f1,4,7,20,21 INPUT_AA_COMBINED.hmmscan.sort.out | while read LINE ; do
			MODEL=$( echo "$LINE" | cut -f1 )
			AA_NAME=$( echo "$LINE" | cut -f2 ) ; 
			LOW_START=$( grep "${AA_NAME}" SPLIT_INPUT_AA_*hmmscan.out | grep "$MODEL" | sed 's/ \+/	/g' | cut -f20 | sort -g | head -n1 )
			HIGH_END=$( grep "${AA_NAME}" SPLIT_INPUT_AA_*hmmscan.out | grep "$MODEL" | sed 's/ \+/	/g' | cut -f21 | sort -rg | head -n1 )
			LENGTH=$(( $HIGH_END-$LOW_START ))
			### needs a bioawk formatted AA seq file
			grep -A1 "$AA_NAME" ${INPUT_AA%.faa}.sort.faa | bioawk -v STQ="$LOW_START" -v LNQ="$LENGTH" -v ENDQ="$HIGH_END" -c fastx '{print ">"$name"_@#@_"STQ"-"ENDQ ; print substr($seq, STQ, LNQ)}' ; 
		done > Preliminary_DNA_methylase_calls1.trimmed_to_hit.faa

	else
		echo "No preliminary DNA methylase hits found. Exiting."
		exit
	fi
else
	echo "is $INPUT_AA misformatted or empty? Exiting."
	exit
fi

if [ -s Preliminary_DNA_methylase_calls1.trimmed_to_hit.faa ] ; then
	PRELIM_AA_SEQS=$( grep -F ">" Preliminary_DNA_methylase_calls1.trimmed_to_hit.faa | wc -l | bc )
	if [ $PRELIM_AA_SEQS -ge 1 ] ; then 
		AA_SEQS_PER_FILE=$( echo "scale=0 ; $PRELIM_AA_SEQS / $CPU" | bc )
		if [ $AA_SEQS_PER_FILE = 0 ] ; then
			AA_SEQS_PER_FILE=1
		fi
		awk -v seq_per_file="$AA_SEQS_PER_FILE" 'BEGIN {n_seq=0;} /^>/ {if(n_seq%seq_per_file==0){file=sprintf("SPLIT_PRELIM_AA_%d.fasta",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < Preliminary_DNA_methylase_calls1.trimmed_to_hit.faa
		SPLIT_AA_PRELIM=$( find * -maxdepth 0 -type f -name "SPLIT_PRELIM_AA_*.fasta" )
		MDYT=$( date +"%m-%d-%y---%T" )
		echo "Step 3. Cross-checking preliminary DNA methylase domain hits against all CDD -- $MDYT"
		echo "$SPLIT_AA_PRELIM" | sed 's/.fasta//g' | xargs -n 1 -I {} -P $CPU -t hmmscan --domtblout {}.hmmscan.out --cpu 1 -E 1e-2 --noali $CDD_PLUS_HMMs {}.fasta >/dev/null 2>&1
		cat SPLIT_PRELIM_AA_*.hmmscan.out | grep -v "^#" | sed 's/ \+/	/g' | sort -u -k4,4 > PRELIM_AA_VS_CDD_PLUS.hmmscan.sort.out
		if [ -s PRELIM_AA_VS_CDD_PLUS.hmmscan.sort.out ] ; then
			grep -f $LEGIT_DOMAIN_LIST PRELIM_AA_VS_CDD_PLUS.hmmscan.sort.out > AA_seqs_with_legit_DNA_methylase_domain.hmmscan.sort.out
			if [ -s AA_seqs_with_legit_DNA_methylase_domain.hmmscan.sort.out ] ; then
				cat AA_seqs_with_legit_DNA_methylase_domain.hmmscan.sort.out | cut -f4 | sed 's/_@#@_.*//g' | while read LINE ; do
					grep -A1 "${LINE}" ${INPUT_AA%.faa}.sort.faa
				done > ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa
				echo "One or more DNA methylase genes is predicted. See ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa"
			else
				echo "AA_seqs_with_legit_DNA_methylase_domain.hmmscan.sort.out not found"
			fi
		else
			echo "PRELIM_AA_VS_CDD_PLUS.hmmscan.sort.out not found."
		fi
	else
		echo "No seqs found in prelinary hits .faa file"
	fi
else
	echo "Preliminary_DNA_methylase_calls1.trimmed_to_hit.faa not found"
fi

rm -f SPLIT_INPUT_AA* SPLIT_PRELIM_*

MDYT=$( date +"%m-%d-%y---%T" )

echo "Done finding DNA methylases."
echo "Step 4. Looking for adjacent DNA methylase gene fragments to combine into 1 gene -- $MDYT"

if [ ${INPUT_TYPE} == "nucl" ] ; then
	if [ -s ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa ] ; then
		echo " " >> IDd_as_frags1.txt
		grep ">" ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa | sed 's/>//g' | while read LINE ; do
			##
			if grep -q ">${LINE}" IDd_as_frags1.txt ; then
				echo "already merged $LINE"
			else
				QSTRAND=$( grep "^>${LINE%_} " ${INPUT_AA} | sed 's/ //g' | cut -d "#" -f4 )
				NUMB=$( echo "$LINE" | sed 's/.*_\([0-9]\{1,9\}\)_/\1/' ) 
				NONUMB=$( echo "$LINE" | sed 's/\(.*\)_[0-9]\{1,9\}_/\1/' )
				PLUS1=$((${NUMB}+1))
				#echo "$LINE, $NUMB, $PLUS1"
				if grep -q ">${NONUMB}_${PLUS1}_" ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa ; then
					PLUS1STRAND=$( grep "^>${NONUMB}_${PLUS1} " ${INPUT_AA} | sed 's/ //g' | cut -d "#" -f4 )
					if [[ $PLUS1STRAND == $QSTRAND ]] ; then
						PLUS1_frag="yes"
						echo ">${NONUMB}_${PLUS1}_",">${LINE}" >> IDd_as_frags1.txt
					else 
						PLUS1_frag="no"
					fi
				else
					PLUS1_frag="no"
				fi
				PLUS2=$((${NUMB}+2))
				if grep -q ">${NONUMB}_${PLUS2}_" ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa ; then
					PLUS2STRAND=$( grep "^>${NONUMB}_${PLUS2} " ${INPUT_AA} | sed 's/ //g' | cut -d "#" -f4 )
					if [[ $PLUS2STRAND == $QSTRAND ]] ; then
						PLUS2_frag="yes"
						echo ">${NONUMB}_${PLUS2}_",">${LINE}" >> IDd_as_frags1.txt
					else 
						PLUS2_frag="no"
					fi
				else
					PLUS2_frag="no"
				fi
				#echo "$LINE, $PLUS1_frag, $PLUS2_frag"
				if [ "$PLUS1_frag" == "yes" ] && [ "$PLUS2_frag" == "yes" ] ; then 
					###
					LINE_AA=$( grep -A1 ">${LINE}" ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa | tail -n1 | sed 's/\*//g' )
					PLUS1_AA=$( grep -A1 ">${NONUMB}_${PLUS1}_" ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa | tail -n1 | sed 's/\*//g' )
					PLUS2_AA=$( grep -A1 ">${NONUMB}_${PLUS2}_" ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa | tail -n1 | sed 's/\*//g' )
					if [ "$QSTRAND" == 1 ] ; then
						echo "$LINE" | awk -v LAA="$LINE_AA" -v P1AA="$PLUS1_AA" -v P2AA="$PLUS2_AA" -v lq="$LINE" -v PL1="${NONUMB}_${PLUS1}_" -v PL2="${NONUMB}_${PLUS2}_" '{ print ">"lq"@"PL1"@"PL2"#merged" ; print LAA"XXXX"P1AA"XXXX"P2AA }' >> ${INPUT_AA%.faa}.DNA_methylases.merged_frags.faa
					elif [ "$QSTRAND" == -1 ] ; then
						echo "$LINE" | awk -v LAA="$LINE_AA" -v P1AA="$PLUS1_AA" -v P2AA="$PLUS2_AA" -v lq="$LINE" -v PL1="${NONUMB}_${PLUS1}_" -v PL2="${NONUMB}_${PLUS2}_" '{ print ">"lq"@"PL1"@"PL2"#merged" ; print P2AA"XXXX"P1AA"XXXX"LAA }' >> ${INPUT_AA%.faa}.DNA_methylases.merged_frags.faa
					else
						echo "orientation issue $LINE"
					fi
				elif [ "$PLUS1_frag" == "yes" ] ; then
					LINE_AA=$( grep -A1 ">${LINE}" ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa | tail -n1 | sed 's/\*//g' )
					PLUS1_AA=$( grep -A1 ">${NONUMB}_${PLUS1}_" ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa | tail -n1 | sed 's/\*//g' )
					if [ "$QSTRAND" == 1 ] ; then
						echo "$LINE" | awk -v LAA="$LINE_AA" -v P1AA="$PLUS1_AA" -v lq="$LINE" -v PL1="${NONUMB}_${PLUS1}_" '{ print ">"lq"@"PL1"#merged" ; print LAA"XXXX"P1AA }' >> ${INPUT_AA%.faa}.DNA_methylases.merged_frags.faa
					elif [ "$QSTRAND" == -1 ] ; then
						echo "$LINE" | awk -v LAA="$LINE_AA" -v P1AA="$PLUS1_AA" -v lq="$LINE" -v PL1="${NONUMB}_${PLUS1}_" '{ print ">"lq"@"PL1"#merged" ; print P1AA"XXXX"LAA }' >> ${INPUT_AA%.faa}.DNA_methylases.merged_frags.faa
					else
						echo "orientation issue $LINE"
					fi
				elif [ "$PLUS2_frag" == "yes" ] ; then 
					LINE_AA=$( grep -A1 ">${LINE}" ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa | tail -n1 | sed 's/\*//g' )
					PLUS2_AA=$( grep -A1 ">${NONUMB}_${PLUS2}_" ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa | tail -n1 | sed 's/\*//g' )
					if [ "$QSTRAND" == 1 ] ; then
						echo "$LINE" | awk -v LAA="$LINE_AA" -v P2AA="$PLUS2_AA" -v lq="$LINE" -v PL2="${NONUMB}_${PLUS2}_" '{ print ">"lq"@"PL2"#merged" ; print LAA"XXXX"P2AA }' >> ${INPUT_AA%.faa}.DNA_methylases.merged_frags.faa
					elif [ "$QSTRAND" == -1 ] ; then
						echo "$LINE" | awk -v LAA="$LINE_AA" -v P2AA="$PLUS2_AA" -v lq="$LINE" -v PL2="${NONUMB}_${PLUS2}_" '{ print ">"lq"@"PL2"#merged" ; print P2AA"XXXX"LAA }' >> ${INPUT_AA%.faa}.DNA_methylases.merged_frags.faa
					else
						echo "orientation issue $LINE"
					fi
				fi
			fi
		done
		grep ">" ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa | sed 's/>//g' | while read LINE ; do
			if grep -q ">${LINE}" IDd_as_frags1.txt ; then
				echo "$LINE is being merged"
			else
				grep -A1 ">${LINE}" ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa | sed '/--/d' >> ${INPUT_AA%.faa}.DNA_methylases.combined.faa
			fi
		done
		if [ -s ${INPUT_AA%.faa}.DNA_methylases.merged_frags.faa ] ; then
			cat ${INPUT_AA%.faa}.DNA_methylases.merged_frags.faa >> ${INPUT_AA%.faa}.DNA_methylases.combined.faa
		fi
	else
		echo "${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa not found"
	fi
elif [ ${INPUT_TYPE} == "AA" ] ; then
	echo "You chose AA input, so fragments will not be joined"
	cp ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa ${INPUT_AA%.faa}.DNA_methylases.combined.faa
fi


MDYT=$( date +"%m-%d-%y---%T" )
echo "Step 5.  Guessing (low confidence) motif specificity and subtype of putative methylases -- $MDYT"

if [ -s ${INPUT_AA%.faa}.DNA_methylases.combined.faa ] ; then 
	echo "Motif guessing part"
	blastp -query ${INPUT_AA%.faa}.DNA_methylases.combined.faa -db $MOTIF_ANNOTATE_BLASTP -outfmt '6 std qlen slen' -max_target_seqs 100 -num_threads $CPU -out ${INPUT_AA%.faa}.DNA_methylases.combined.w_motifs.blastp.out
	python ${FIND_METHYLASE_DIR}/anicalc.py -i ${INPUT_AA%.faa}.DNA_methylases.combined.w_motifs.blastp.out -o ${INPUT_AA%.faa}.DNA_methylases.combined.w_motifs.aaicalc.out
	awk -v pidq="$PID" -v covq="$COV" '{OFS="\t"}{FS="\t"}{ if ($4>=pidq && $5>=covq) {print}}' ${INPUT_AA%.faa}.DNA_methylases.combined.w_motifs.aaicalc.out | while read LINE ; do 
		MOTIF=$( echo "$LINE" | sed 's/.*___\([A-Z]\{2,20\}\)__.*/\1/' ) ; 
		echo "$LINE" | awk -v motq="$MOTIF" '{OFS="\t"}{FS="\t"}{ if (/qname/) {print $1, "Motif", $2, $3, $4, $5, $6} else {print $1, motq, $2, $3, $4, $5, $6}}'  ; 
	done > ${INPUT_AA%.faa}.DNA_methylases.combined.w_motifs.${PID}_${COV}.out
	sort -u -k1,1 ${INPUT_AA%.faa}.DNA_methylases.combined.w_motifs.${PID}_${COV}.out > ${INPUT_AA%.faa}.DNA_methylases.combined.w_motifs.${PID}_${COV}.top.out
	echo "Subtype guessing part"
	hmmscan --tblout ${INPUT_AA%.faa}.DNA_methylases.combined.type_hmmscan.out --cpu $CPU -E 1e-3 --noali $SUBTYPE_ANNOTATE_HMM ${INPUT_AA%.faa}.DNA_methylases.combined.faa >/dev/null 2>&1
	grep -v "^#" ${INPUT_AA%.faa}.DNA_methylases.combined.type_hmmscan.out | sed 's/ \+/	/g' | sort -u -k3,3 > ${INPUT_AA%.faa}.DNA_methylases.combined.type_hmmscan.sort.out
else
	echo "Combined Predicted Methylase .faa file not found."
fi

MDYT=$( date +"%m-%d-%y---%T" )
echo "Step 6.  Making summary table -- $MDYT"

if [ -s ${INPUT_AA%.faa}.DNA_methylases.combined.faa ] && [ -s ${INPUT_AA%.faa}.DNA_methylases.combined.type_hmmscan.sort.out ] && [ -s ${INPUT_AA%.faa}.DNA_methylases.combined.w_motifs.${PID}_${COV}.top.out ] ; then
	echo "gene name	sub-type	contig	start position	end position	orientation	motif guess	AAI to top hit	AF to top hit" > ${INPUT_AA%.faa}.DNA_methylases.combined.summary.tsv
	grep ">" ${INPUT_AA%.faa}.DNA_methylases.combined.faa | sed 's/>//g' | while read LINE ; do 
		if grep -q "$LINE" ${INPUT_AA%.faa}.DNA_methylases.combined.type_hmmscan.sort.out ; then 
			SUBTYPE=$( grep "$LINE" ${INPUT_AA%.faa}.DNA_methylases.combined.type_hmmscan.sort.out | cut -f1 | sed 's/@@.*//g' ) ; 
		else 
			SUBTYPE="Unknown" ; 
		fi ; 
		
		if [ -s ${INPUT_NUCL%.fna}.genes.fna ] ; then
			if echo "$LINE" | grep -q "#merged" ; then
				NUMBER_ATS=$( echo "$LINE" | grep -o "@" | wc -l )
				#echo "$LINE $NUMBER_ATS"
				CONTIG=$( echo "$LINE" | cut -d "@" -f1 | sed 's/\(.*\)_[0-9]\{1,9\}_/\1/' )
				if [ "$NUMBER_ATS" == 1 ] ; then
					echo "$LINE" | sed 's/#merged//g ; s/@/ /g' | while read ONE TWO ; do
						STARTQ=$( grep -e "${ONE::-1} # " -e "${TWO::-1} # " ${INPUT_NUCL%.fna}.genes.fna | sed 's/ //g' | cut -d "#" -f2 | sort -g | head -n1 )
						ENDQ=$( grep -e "${ONE::-1} # " -e "${TWO::-1} # " ${INPUT_NUCL%.fna}.genes.fna | sed 's/ //g' | cut -d "#" -f3 | sort -g | tail -n1 )
						ORIENT=$( grep "${ONE::-1} # " ${INPUT_NUCL%.fna}.genes.fna | sed 's/ //g' | cut -d "#" -f4 | sed 's/-1/-/g ; s/1/+/g' )
					done
				else
					echo "$LINE" | sed 's/#merged//g ; s/@/ /g' | while read ONE TWO THREE ; do
						STARTQ=$( grep -e "${ONE::-1} # " -e "${TWO::-1} # " -e "${THREE::-1} # " ${INPUT_NUCL%.fna}.genes.fna | sed 's/ //g' | cut -d "#" -f2 | sort -g | head -n1 )
						ENDQ=$( grep -e "${ONE::-1} # " -e "${TWO::-1} # " -e "${THREE::-1} # " ${INPUT_NUCL%.fna}.genes.fna | sed 's/ //g' | cut -d "#" -f3 | sort -g | tail -n1 )
						ORIENT=$( grep "${ONE::-1} # " ${INPUT_NUCL%.fna}.genes.fna | sed 's/ //g' | cut -d "#" -f4 | sed 's/-1/-/g ; s/1/+/g' )
					done
				fi
			else
				CONTIG=$( echo "$LINE" | sed 's/\(.*\)_[0-9]\{1,9\}_/\1/' )
				PROD_FMT=${LINE::-1}
				STARTQ=$( grep "${PROD_FMT} # " ${INPUT_NUCL%.fna}.genes.fna | sed 's/ //g' | cut -d "#" -f2 )
				ENDQ=$( grep "${PROD_FMT} # " ${INPUT_NUCL%.fna}.genes.fna | sed 's/ //g' | cut -d "#" -f3 )
				ORIENT=$( grep "${PROD_FMT} # " ${INPUT_NUCL%.fna}.genes.fna | sed 's/ //g' | cut -d "#" -f4 | sed 's/-1/-/g ; s/1/+/g' )

			fi
		else
			STARTQ="NA"
			ENDQ="NA"
			ORIENT="NA"
		fi

		if grep -q "$LINE" ${INPUT_AA%.faa}.DNA_methylases.combined.w_motifs.${PID}_${COV}.top.out ; then 
			MOTIF=$( grep "$LINE" ${INPUT_AA%.faa}.DNA_methylases.combined.w_motifs.${PID}_${COV}.top.out | awk '{OFS="\t"}{FS="\t"}{print $2, $5" AAI", $6" AF"}' ) ; 
		else 
			MOTIF="not_found	not_found	not_found" ; 
		fi ; 
		echo -e "$LINE\t$SUBTYPE\t${CONTIG}\t${STARTQ}\t${ENDQ}\t${ORIENT}\t$MOTIF" ; 
	done >> ${INPUT_AA%.faa}.DNA_methylases.combined.summary.tsv
else 
	echo "not all files required for summary table found."
fi


if [ "$NEIGHBORHOODS" == "True" ] ; then
	MDYT=$( date +"%m-%d-%y---%T" )
	echo "Step 7. Annotating DNA methylase gene neighborhoods -- $MDYT"

	if [ ${INPUT_TYPE} == "nucl" ] ; then
		if [ -s ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa ] ; then
			###
			mkdir neighborhood_annotations
			grep ">" ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa | sed 's/>//g' | while read LINE ; do
				NONUMB=$( echo "$LINE" | sed 's/\(.*\)_[0-9]\{1,9\}_/\1/' )
				bioawk -c fastx '{print ">"$name, $seq}' ${INPUT_AA%.faa}.sort.faa | grep -C5 "${LINE}" | grep "^>${NONUMB}_" | awk '{ print $1 ; print $2}' > neighborhood_annotations/${LINE}neighbors.faa
				hmmscan --tblout neighborhood_annotations/${LINE}neighbors.S_subunit.out --cpu $CPU -E 1e-5 --noali $S_SUBUNIT_HMM neighborhood_annotations/${LINE}neighbors.faa >/dev/null 2>&1
				hmmscan --tblout neighborhood_annotations/${LINE}neighbors.RE.out --cpu $CPU -E 1e-5 --noali $RE_HMM neighborhood_annotations/${LINE}neighbors.faa >/dev/null 2>&1
				grep ">" neighborhood_annotations/${LINE}neighbors.faa | sed 's/>//g' | while read NEIGHBOR ; do
					if [ "$NEIGHBOR" != "$LINE" ] ; then
						###
						if grep -q "$NEIGHBOR" neighborhood_annotations/${LINE}neighbors.S_subunit.out neighborhood_annotations/${LINE}neighbors.RE.out ; then
							continue
						else
							grep -A1 "$NEIGHBOR" neighborhood_annotations/${LINE}neighbors.faa >> neighborhood_annotations/${LINE}neighbors.for_CDD.faa
						fi
					fi
				done
				if [ -s neighborhood_annotations/${LINE}neighbors.for_CDD.faa ] ; then
					hmmscan --tblout neighborhood_annotations/${LINE}neighbors.CDD.out --cpu $CPU -E 1e-5 --noali $CDD_PLUS_HMMs neighborhood_annotations/${LINE}neighbors.for_CDD.faa >/dev/null 2>&1
				fi

				HMM_OUTs=$( find .  -type f -name "${LINE}neighbors*.out" )
				if [ -n "$HMM_OUTs" ] ; then
					cat neighborhood_annotations/${LINE}neighbors*.out | grep -v "^#" | sed 's/ \+/	/g' | sort -u -k3,3 > neighborhood_annotations/${LINE}neighbors.combined_HMM.sort.out
				fi
				# make neighborhood fasta and gff part
				echo "##gff-version 3" > neighborhood_annotations/${LINE}neighbors.gff
				START_FNA=$( grep ">" neighborhood_annotations/${LINE}neighbors.faa | while read HEADER ; do grep "${HEADER%_} " ${INPUT_AA} ; done | sed 's/ //g' | cut -d "#" -f2 | head -n1 | bc ) ; 
				END_FNA=$( grep ">" neighborhood_annotations/${LINE}neighbors.faa | while read HEADER ; do grep "${HEADER%_} " ${INPUT_AA} ; done | sed 's/ //g' | cut -d "#" -f3 | tail -n1 | bc )
				END1_FNA=$((${END_FNA}+1))
				N_CONTIG=$( echo "$LINE" | sed 's/\(.*\)_[0-9]\{1,9\}_$/\1/' )

				bioawk -vstq="$START_FNA" -venq="$END1_FNA" -vneigh="$N_CONTIG" -c fastx '{lenq=(enq-stq)}{ if ($name == neigh) {print ">"$name"@"stq"-"enq ; print substr($seq, stq, lenq)}}' ../${INPUT_NUCL} > neighborhood_annotations/${LINE}neighbors.fna
				grep ">" neighborhood_annotations/${LINE}neighbors.faa | sed 's/>//g' | while read HEADER ; do
					HSTART=$( grep "${HEADER%_} " ${INPUT_AA} | sed 's/ //g' | cut -d "#" -f2 )
					HEND=$( grep "${HEADER%_} " ${INPUT_AA} | sed 's/ //g' | cut -d "#" -f3 )
					NSTART=$(( $HSTART-$START_FNA+1 ))
					NEND=$(( $HEND-$START_FNA+1 ))
					HSTRAND=$( grep "${HEADER%_} " ${INPUT_AA} | sed 's/ //g' | cut -d "#" -f4 | sed 's/-1/-/g ; s/1/+/g' )
					NUMB=$( echo "$HEADER" | sed 's/.*_\([0-9]\{1,9\}\)_/\1/' ) 
					NONUMB=$( echo "$HEADER" | sed 's/\(.*\)_[0-9]\{1,9\}_/\1/' )
					MINUS1=$((${NUMB}-1))
					MINUS2=$((${NUMB}-2))
					if grep -q "^${HEADER}" ${INPUT_AA%.faa}.DNA_methylases.combined.summary.tsv ; then
						NATT=$( grep "^${HEADER}" ${INPUT_AA%.faa}.DNA_methylases.combined.summary.tsv | awk '{FS="\t"}{OFS="\t"}{ print "DNA Methylase of subtype: "$2" motif guess: "$7 }' )
						echo -e "${N_CONTIG}@${START_FNA}-${END1_FNA}\tfind_methylases\tCDS\t${NSTART}\t${NEND}\t.\t${HSTRAND}\t.\tcodon_start=1;transl_table=11;product=${NATT};ID=${HEADER}" >> neighborhood_annotations/${LINE}neighbors.gff
					elif grep -q "^${NONUMB}_${MINUS1}_" ${INPUT_AA%.faa}.DNA_methylases.combined.summary.tsv && grep -q ">${HEADER}" ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa ; then
						NATT=$( grep "^${NONUMB}_${MINUS1}_" ${INPUT_AA%.faa}.DNA_methylases.combined.summary.tsv | awk '{FS="\t"}{OFS="\t"}{ print "DNA Methylase of subtype "$2" motif guess- "$7 }' )
						echo -e "${N_CONTIG}@${START_FNA}-${END1_FNA}\tfind_methylases\tCDS\t${NSTART}\t${NEND}\t.\t${HSTRAND}\t.\tcodon_start=1;transl_table=11;product=${NATT};ID=${HEADER}" >> neighborhood_annotations/${LINE}neighbors.gff
					elif grep -q "^${NONUMB}_${MINUS2}_" ${INPUT_AA%.faa}.DNA_methylases.combined.summary.tsv && grep -q ">${HEADER}" ${INPUT_AA%.faa}.DNA_methylases.final_prediction.faa ; then
						NATT=$( grep "^${NONUMB}_${MINUS2}_" ${INPUT_AA%.faa}.DNA_methylases.combined.summary.tsv | awk '{FS="\t"}{OFS="\t"}{ print "DNA Methylase of subtype: "$2" motif guess: "$7 }' )
						echo -e "${N_CONTIG}@${START_FNA}-${END1_FNA}\tfind_methylases\tCDS\t${NSTART}\t${NEND}\t.\t${HSTRAND}\t.\tcodon_start=1;transl_table=11;product=${NATT};ID=${HEADER}" >> neighborhood_annotations/${LINE}neighbors.gff
					elif grep -q "${HEADER}" neighborhood_annotations/${LINE}neighbors.S_subunit.out ; then
						NATT="Putative Specificity Subunit Gene"
						echo -e "${N_CONTIG}@${START_FNA}-${END1_FNA}\tfind_methylases\tCDS\t${NSTART}\t${NEND}\t.\t${HSTRAND}\t.\tcodon_start=1;transl_table=11;product=${NATT};ID=${HEADER}" >> neighborhood_annotations/${LINE}neighbors.gff
					elif grep -q "${HEADER}" neighborhood_annotations/${LINE}neighbors.RE.out ; then
						NATT="Putative Restrction Endonuclease"
						echo -e "${N_CONTIG}@${START_FNA}-${END1_FNA}\tfind_methylases\tCDS\t${NSTART}\t${NEND}\t.\t${HSTRAND}\t.\tcodon_start=1;transl_table=11;product=${NATT};ID=${HEADER}" >> neighborhood_annotations/${LINE}neighbors.gff
					elif grep -q "${HEADER}" neighborhood_annotations/${LINE}neighbors.combined_HMM.sort.out ; then
						NATT=$( grep "${HEADER}" neighborhood_annotations/${LINE}neighbors.combined_HMM.sort.out | cut -f1 | cut -d "/" -f 2- | sed 's/;/ /g ; s/-/ /g' )
						echo -e "${N_CONTIG}@${START_FNA}-${END1_FNA}\tfind_methylases\tCDS\t${NSTART}\t${NEND}\t.\t${HSTRAND}\t.\tcodon_start=1;transl_table=11;product=${NATT};ID=${HEADER}" >> neighborhood_annotations/${LINE}neighbors.gff
					else
						NATT="Hypothetical protein"
						echo -e "${N_CONTIG}@${START_FNA}-${END1_FNA}\tfind_methylases\tCDS\t${NSTART}\t${NEND}\t.\t${HSTRAND}\t.\tcodon_start=1;transl_table=11;product=${NATT};ID=${HEADER}" >> neighborhood_annotations/${LINE}neighbors.gff
					fi

				done
				if [ -s neighborhood_annotations/${LINE}neighbors.gff ] && [ -s neighborhood_annotations/${LINE}neighbors.fna ] ; then
					python ${FIND_METHYLASE_DIR}/gff_to_genbank2.py neighborhood_annotations/${LINE}neighbors.gff neighborhood_annotations/${LINE}neighbors.fna
					if [ -s neighborhood_annotations/${LINE}neighbors.gb ] ; then
						echo "neighborhood gene map created @ neighborhood_annotations/${LINE}neighbors.gb"
					fi
				fi
			done
		else
			echo "DNA methylase predictions not found"
		fi
	else
		echo "input type was not nucleotide sequence, so neighborhoods were not determined"
	fi
fi

echo "##### Done with all tasks #####"
echo "See output at: $OUTPUT_DIRECTORY"
MDYT=$( date +"%m-%d-%y---%T" )
echo "End time -- $MDYT"


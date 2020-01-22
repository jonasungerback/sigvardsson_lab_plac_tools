#!/bin/bash

#Takes a 2D bed file consisting of 6 columns and a reference genome as two arguments and annotate these with Homer and creates a new table with the Homer annotated anchor-points.

# $1 = 2D-bed file.
# $2 = homer reference genome
LOG_FILE=script.log
exec > >(while read -r line; do printf '%s %s\n' "$(date --rfc-3339=seconds)" "$line" | tee -a $LOG_FILE; done)
exec 2> >(while read -r line; do printf '%s %s\n' "$(date --rfc-3339=seconds)" "$line" | tee -a $LOG_FILE; done >&2)
#Test so the variables are not empty
if [ -z "$1" ] || [ -z "$2" ]; then echo "At least one of the variables is missing." && exit

else

    echo "You are good to go!!!"
    echo ""

    # Take out the columns important for the interactions. The if statement checks if there is a string called start.

    if grep -lq -e start $1
        then
        # code if found. Then remove line 1 to make it Homer compatible.
            echo "############################"
            echo "The word start found in line 1"
            echo "############################"
            awk 'BEGIN {FS=OFS="\t"} ; { if (NR!=1) print $1,$2,$3,$4,$5,$6}' $1 > interaction_file
        else
        # Print as is without removing line 1.
            echo "############################"
            echo "The word start NOT found in line 1"
            echo "############################"
            awk 'BEGIN {FS=OFS="\t"} ; {  print $1,$2,$3,$4,$5,$6}' $1 > interaction_file
        fi

	#These two commands creates two new tab-separated tables with columns appearing in the following order chr start end "anchor#" contact count isPeak (1=yes, 0=no) P-value and Q-value. Must print
	
	mkdir -p $1.$2.annotated_interactions

    awk 'BEGIN {FS=OFS="\t"} ; { print $1,$2,$3,"anchor1"}' interaction_file > $1.$2.annotated_interactions/anchor1.bed
    awk 'BEGIN {FS=OFS="\t"} ; { print $4,$5,$6,"anchor2"}' interaction_file > $1.$2.annotated_interactions/anchor2.bed

    rm interaction_file

	#Annotate the files with homer annotatePeaks
	annotatePeaks.pl $1.$2.annotated_interactions/anchor1.bed $2 > $1.$2.annotated_interactions/anchor1.$2.txt
	annotatePeaks.pl $1.$2.annotated_interactions/anchor2.bed $2 > $1.$2.annotated_interactions/anchor2.$2.txt

	#Remove the first line of the homer output and sort anchor1 output based on columns 1. Hopefully this gives the same order as the anchor file. This sort based on the name of the first column which is anchor-1 .. anchor-n. Resort this and it recreates the order in the original anchor1.bed
    awk 'BEGIN {FS=OFS="\t"} ; { if (NR!=1) print }' $1.$2.annotated_interactions/anchor1.$2.txt|sort -k1,1V > $1.$2.annotated_interactions/anchor1.$2.s.txt
    #Repeat for anchor 2.
    awk 'BEGIN {FS=OFS="\t"} ; { if (NR!=1) print }' $1.$2.annotated_interactions/anchor2.$2.txt|sort -k1,1V > $1.$2.annotated_interactions/anchor2.$2.s.txt


	#Merge the wanted columns from the different data frames to a new file.
	# Peak1 chr_1 s1 e1 strand1 cc_1 isPeak1 Ann1 Dist_TSS_1 refseq1 Ensembl1 gene1 Peak2 chr_2 s2 e2 strand2 cc_2 isPeak2 Ann2 Dist_TSS_2 refseq1 Ensembl2 gene2 P-value Q-value
	paste $1.$2.annotated_interactions/anchor1.$2.s.txt $1.$2.annotated_interactions/anchor2.$2.s.txt $1.$2.annotated_interactions/anchor1.bed $1.$2.annotated_interactions/anchor2.bed |awk 'BEGIN {FS=OFS="\t"};{ print $1,$39, $40, $41, $5, $8, $10,$11, $12, $13, $14, $15, $16, $19, $20, $43, $44, $45, $24, $27, $29, $30, $31, $32, $33, $34, $35 ,$38 }' >  $1.$2.annotated_interactions/$1.$2.annotated.tmp.txt

	touch $1.$2.annotated_interactions/colheader_temp.txt
	touch $1.$2.annotated_interactions/col_header.txt

	{ printf '#Peak1\tchr_1\ts1\te1\tstrand1\tAnn1\tDist_TSS1\tNearest_Promoter_1\tEntrezID_1\tUnigene_1\trefseq_1\tEnsembl_1\tgene_1\tGene_type_1\tPeak2\tchr_2\ts2\te2\tstrand2\tAnn2\tDist_TSS2\tNearest_Promoter_2\tEntrezID_2\tUnigene_2\trefseq_2\tEnsembl_2\tgene_2\tGene_type_2\n'; cat $1.$2.annotated_interactions/colheader_temp.txt; } > $1.$2.annotated_interactions/col_header.txt

	rm $1.$2.annotated_interactions/colheader_temp.txt

	cat $1.$2.annotated_interactions/col_header.txt $1.$2.annotated_interactions/$1.$2.annotated.tmp.txt > $1.$2.annotated_interactions/$1.$2.annotated.txt

    rm $1.$2.annotated_interactions/$1.$2.annotated.tmp.txt $1.$2.annotated_interactions/col_header.txt $1.$2.annotated_interactions/anchor*

    echo "Annotation complete!!!!"

fi


#!/bin/bash

#Takes a 2D bed file [chr1,start1,end1,chr2,start2,end2, numerical] with one extra column and a Homer annotated peak file and pulls out the interactions with a peak in one or the other end.

# Updated to v.2 180822. Now takes a Bias corrected file from the updated FitHiChIP pipeline.
# Updated to v.3 181008. This version does not print the TF peak coordinates but convert them to a two binary columns telling weather an anchor point matched a TF.
#Updated to v.4  181028. Also saving the bias_corrected_cc.

# $1 = 2D-bed file
# $2 = Homer annotated peak file
# $3 = Name of TF as a string

#Example: bash 2Dbed_overlap_Homer_annotated_peaks_v2.sh 2D-BED.bed homer.annotated.txt TF

#Test so the variables are not empty
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then echo "At least one of the variables is missing." && exit

else

    # Turn the interaction file seven column format
    awk 'BEGIN {FS=OFS="\t"} ; {  print $1,$2,$3,$4,$5,$6,$7}' $1 > interaction_file

	#Turn Homer file to a bed file:
	awk 'BEGIN {FS=OFS="\t"} ; { if (NR!=1) print $2,$3,$4,$5}' $2 > $2.bed

	#These two commands creates two new tab-separated tables with columns appearing in the following order chr start end "anchor#" contact count isPeak (1=yes, 0=no) P-value and Q-value
	awk 'BEGIN {FS=OFS="\t"} ; { if (NR!=1) print $1,$2,$3}' $1 > anchor1.tmp1.bed
	awk 'BEGIN {FS=OFS="\t"} ; { if (NR!=1) print $4,$5,$6}' $1 > anchor2.tmp1.bed

    #Add an unique counter to both end points.
    # Anchor1
    awk '{printf "peak_%d\t%s\n", NR,$0}' < anchor1.tmp1.bed |cut -f 1 > anchor1.tmp2.bed
    paste anchor1.tmp1.bed anchor1.tmp2.bed > anchor1.bed

    # Anchor2
    awk '{printf "peak_%d\t%s\n", NR,$0}' < anchor2.tmp1.bed |cut -f 1 > anchor2.tmp2.bed
    paste anchor2.tmp1.bed anchor2.tmp2.bed > anchor2.bed
    rm *tmp*

    #Remove first line from interaction file:
    awk 'BEGIN {FS=OFS="\t"} ; { if (NR!=1) print }' interaction_file > interaction_file.tmp

    #Add unique id to interaction file
    awk '{printf "peak_%d\t%s\n", NR,$0}' < interaction_file.tmp > interaction_file.id

    #Sort the interaction file
    sort interaction_file.id > interaction_file.s.id
    #Make the interaction file with the Peak id as the first column.
    scp interaction_file.s.id $1.$3.id.txt
	rm interaction_file.tmp interaction_file.id

    echo "################################################################################################"
    echo "######## This part calculates the TF overlaps but MAINTAINS the interaction coordinates. #######"
    echo "################################################################################################"
    #Try to make a function of the block below.
	#use bedtools to take out binding interactions and saving the positions for the peak files.
    #Change -wa to -wb to write the TF coordinates to the final file.
	bedtools intersect -wa -a anchor1.bed -b $2.bed > anchor1.peak.bed
	bedtools intersect -wa -a anchor2.bed -b $2.bed > anchor2.peak.bed

	#Resort and only take unique entries
	awk 'BEGIN {OFS="\t"} ; { print $4,$1,$2,$3}' anchor1.peak.bed |sort -u -k1,1 > anchor1.peak.sorted.bed
	awk 'BEGIN {OFS="\t"} ; { print $4,$1,$2,$3}' anchor2.peak.bed |sort -u -k1,1 > anchor2.peak.sorted.bed

	#Match the anchor back to the interaction file for both anchors
	awk 'NR==FNR{a[$1];next} ($1 in a)' anchor1.peak.sorted.bed interaction_file.s.id > $1.s.anchor1.id
	awk 'NR==FNR{a[$1];next} ($1 in a)' anchor2.peak.sorted.bed interaction_file.s.id > $1.s.anchor2.id

    #Replace columne #3 and #4 in $1.s.anchor1.id with column #3 and #4 in anchor1.peak.sorted.bed.
    paste anchor1.peak.sorted.bed $1.s.anchor1.id|awk 'BEGIN {OFS="\t"} ; { print $1,$2,$3,$4,$9,$10,$11,$12}' > tmp1.id
    mv tmp1.id $1.s.anchor1.id

    #Replace columne #6 and #7 in $1.s.anchor2.id with column #3 and #4 in anchor2.peak.sorted.bed.
    paste anchor2.peak.sorted.bed $1.s.anchor2.id|awk 'BEGIN {OFS="\t"} ; { print $5,$6,$7,$8,$2,$3,$4,$12}' > tmp2.id
    mv tmp2.id $1.s.anchor2.id


    #Create a list of all the id's that have a peak in both ends.
    join $1.s.anchor1.id $1.s.anchor2.id -t $'\t'|awk 'BEGIN {OFS="\t"} ; { print $1,$2,$3,$4,$12,$13,$14,$15}'|sort -u -k1,1 > anchor_join

		# Catenate and retrieve the interaction files with peak overlap in one end
	cat $1.s.anchor1.id $1.s.anchor2.id|sort|awk 'BEGIN {OFS="\t"} ; { print $1,$2,$3,$4,$5,$6,$7,$8}'> $1.$3.peak_overlap.tmp.txt

    #Replace the lines in the tmp-file with peaks in both ends and remove double rows.
    awk 'BEGIN {OFS="\t"} ;NR==FNR{a[$1]=$0;next;}a[$1]{$0=a[$1]}1' anchor_join $1.$3.peak_overlap.tmp.txt|sort -u -k1,1 > $1.$3.peak_overlap.tmp2.txt
    mv $1.$3.peak_overlap.tmp2.txt  $1.$3.peak_overlap.tmp.txt
    ####### Save $1.$3.peak_overlap.tmp.txt


    ################################################################################################
    ################################################################################################
    ################################################################################################

    echo ""
    echo "First cleaning up!!!!!!!!"
    echo ""
    rm anchor1.peak.bed anchor1.peak.sorted.bed anchor2.peak.sorted.bed $1.s.anchor1.id $1.s.anchor2.id anchor_join

    echo "################################################################################################"
    echo "########     This part calculates the TF overlaps but CHANGES to TF coordinates.       #########"
    echo "################################################################################################"

    #use bedtools to take out binding interactions and saving the positions for the peak files.
    bedtools intersect -wb -a anchor1.bed -b $2.bed > anchor1.peak.bed
    bedtools intersect -wb -a anchor2.bed -b $2.bed > anchor2.peak.bed

    #Resort and only take unique entries
    awk 'BEGIN {OFS="\t"} ; { print $4,$1,$2,$3}' anchor1.peak.bed |sort -u -k1,1 > anchor1.peak.sorted.bed
    awk 'BEGIN {OFS="\t"} ; { print $4,$1,$2,$3}' anchor2.peak.bed |sort -u -k1,1 > anchor2.peak.sorted.bed

    #Sort the interaction file
    #sort interaction_file.id > interaction_file.s.id

    #Match the anchor back to the interaction file for both anchors
    awk 'NR==FNR{a[$1];next} ($1 in a)' anchor1.peak.sorted.bed interaction_file.s.id > $1.s.anchor1.id
    awk 'NR==FNR{a[$1];next} ($1 in a)' anchor2.peak.sorted.bed interaction_file.s.id > $1.s.anchor2.id

    #Replace columne #3 and #4 in $1.s.anchor1.id with column #3 and #4 in anchor1.peak.sorted.bed.
    paste anchor1.peak.sorted.bed $1.s.anchor1.id|awk 'BEGIN {OFS="\t"} ; { print $1,$2,$3,$4,$9,$10,$11,$12}' > tmp1.id
    mv tmp1.id $1.s.anchor1.id

    #Replace columne #6 and #7 in $1.s.anchor2.id with column #3 and #4 in anchor2.peak.sorted.bed.
    paste anchor2.peak.sorted.bed $1.s.anchor2.id|awk 'BEGIN {OFS="\t"} ; { print $5,$6,$7,$8,$2,$3,$4,$12}' > tmp2.id
    mv tmp2.id $1.s.anchor2.id

    dos2unix $1.s.anchor1.id
    dos2unix $1.s.anchor2.id

    #Create a list of all the id's that have a peak in both ends.
    join $1.s.anchor1.id $1.s.anchor2.id -t $'\t'|awk 'BEGIN {OFS="\t"} ; { print $1,$2,$3,$4,$12,$13,$14,$8}'|sort -u -k1,1 > anchor_join

    # Catenate and retrieve the interaction files with peak overlap in one end
    cat $1.s.anchor1.id $1.s.anchor2.id|sort|awk 'BEGIN {OFS="\t"} ; { print $1,$2,$3,$4,$5,$6,$7,$8}'> $1.$3.peak_overlap.tmp3.txt

    #Replace the lines in the tmp-file with a peaks in both ends and remove double rows.
    awk 'BEGIN {OFS="\t"} ;NR==FNR{a[$1]=$0;next;}a[$1]{$0=a[$1]}1' anchor_join $1.$3.peak_overlap.tmp3.txt|sort -u -k1,1 > TF_anchor_coord

    #Calculate the difference between e1-s1 and e2-s2 for F_anchor_cord
    awk 'BEGIN {OFS="\t"} ; { print $1,$4-$3,$7-$6}' TF_anchor_coord > TF_anchor_coord_diff
    #Make a column where if first coordinate diff is smaller than 5000 print 1 else print0
    awk '$2 < 5000 { print 1 ; } $2 >= 5000 { print 0 ;}' TF_anchor_coord_diff > anchor1_diff
    #Same diff for second coordinate
    awk '$3 < 5000 { print 1 ; } $3 >= 5000 { print 0 ;}' TF_anchor_coord_diff > anchor2_diff
    rm TF_anchor_coord TF_anchor_coord_diff $1.s.anchor1.id $1.s.anchor2.id

    ##############################################################################################
    ##############################################################################################
    ##############################################################################################

    #Merge the diff check onto $1.$3.peak_overlap.tmp.txt
    dos2unix $1.$3.peak_overlap.tmp.txt
    paste $1.$3.peak_overlap.tmp.txt anchor1_diff anchor2_diff > $1.$3.peak_overlap.tmp2.txt

	#Add column names back to interaction file.
        touch colheader_temp.txt
        touch col_header.txt

    { printf 'PeakID\tchr_1\ts1\te1\tchr_2\ts2\te2\tcc\tisTF_left_anchor\tisTF_right_anchor\n'; cat colheader_temp.txt; } > col_header.txt

    cat col_header.txt $1.$3.peak_overlap.tmp2.txt > $1.$3.peak_overlap.tmp3.txt
    awk 'BEGIN {OFS="\t"} ; { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' $1.$3.peak_overlap.tmp3.txt > $1.$3.peak_overlap.txt

    echo "$1.peak_overlap.txt printed to file."
    echo "Final clean-up!"
    rm  $1.$3.peak_overlap.tmp3.txt col_header.txt colheader_temp.txt interaction_file*  $1.$3.peak_overlap.tmp.txt  anchor* $1.$3.peak_overlap.tmp2.txt $2.bed

    #Print number of lines in file to screen.
    echo ""
    echo "The number of rows in the interaction file is: " $(wc -l $1.$3.peak_overlap.txt)
    echo ""

fi


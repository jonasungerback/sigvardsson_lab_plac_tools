# sigvardsson_lab_plac_tools

### 2Dbed_overlap_Homer_annotated_peaks_v2.sh

Takes a 2D bed file [chr1,start1,end1,chr2,start2,end2, numerical] with one extra column and a Homer annotated peak file and pulls out the interactions with a peak in one or the other end.

Example: `bash 2Dbed_overlap_Homer_annotated_peaks_v2.sh 2D-BED.bed homer.annotated.txt TF`

### annotate2Dbed.sh

Takes a 2D bed file consisting of 6 columns and a reference genome as two arguments and annotate these with Homer and creates a new table with the Homer annotated anchor-points.

Example: `bash annotate2Dbed.sh 2D-BED.bed homer_ref_genome`

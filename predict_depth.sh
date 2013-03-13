#!/bin/sh

#written 2-14-13
#if [ $# -eq 0 || $# -gt 3 ] # usage for no or more than 3 arguments                                                                                 
if [ $# -eq 0 -o $# -gt 10 ]
then
    echo "DESCRIPTION: Shell wrapper for predict_depth.r"
    echo 
    echo "USAGE: >predict_depth.sh <file_in> <col_num> <num_reads> <file_out_prefix> <genome_size> <scale_by_unannotated> <read_length> <min_overlap> <produce_pdf> <show>"
    echo   
    
    echo "     <file_in>:         NO DEFAULT; string, name of the input file" 
    echo "                        (two column, tab delimited, annotations then abundances" 
    echo "                        first row, first column (annotation header) is blank, first" 
    echo "                        row, scecond column can have a text header)"

    echo
    echo "     <col_num>:         DEFAULT = 1; 1 based index of column to process from file_in"

    echo
    echo "     <num_reads>        DEFAULT = 1000000; estimated number of reads in wgs sequencing run"

    echo
    echo "     <file_out_prefix>: DEFAULT = \"depth_prediction\"  ;string, prefix for output file(s)" 
    echo "                        *.txt (& *.pdf if produce_fig = 1)"
   
    echo 
    echo "     <genome_size>:     DEFAULT = 4000000; size in bp of genomes" 

    echo
    echo "     <scale_by_unannotated>:   DEFAULT = TRUE; scale calculated depth including unannotated reads"

    echo
    echo "     <read_length>:     DEFAULT = 125; average read length in bp"

    echo
    echo "     <min_overlap>      DEFAULT = 30; min overlap for assembly"
 
    echo
    echo "     <produce_pdf>:     DEFAULT = FALSE;   TRUE|FALSE, produce a pdf visulization of the output"
    
    echo
    echo "     <show>:            DEFAULT = 10; Integer, number of taxa to show in pdf"
    echo
    
    echo "EXAMPLES: predict_depth.sh test_data.txt "
    echo "          predict_depth.sh test_data.txt 1 my_output 4000000 25 125 30 TRUE TRUE 10"
    echo
    echo "NOTES: Annoying -- to change argument 2, you have to supply 1 and 2 etc."
    echo "There are additional arguments in the function file (predict_depth.r)"
    echo
    echo "Prediction values use calculations defined in the following:"
    echo 
    echo "    Lander, E. S. and Waterman, M. S., \"Genomic mapping by fingerprinting"
    echo "    random clones: a mathematical analysis\", Genomics 2, 231-239 (1988)."
    echo


    exit 1                                                                                         
fi

time_stamp=`date +%m-%d-%y_%X`;  # create the time stamp month-day-year_hour:min:sec   # :nanosec removed nano

# grab arguments from prompt
file_in=$1
col_num=$2
num_reads=$3
file_out_prefix=$4
genom_size=$5
scale_by_unannotated=$6
read_length=$7
min_overlap=$8
produce_pdf=$9
show=${10}

# set default values
: ${col_num:=1}
: ${num_reads:=1000000}
: ${file_out_prefix:="depth_prediction"}
: ${genome_size:=4000000}
: ${coverage:=25}
: ${scale_by_unannotated:=0}
: ${read_length:=125}
: ${min_overlap:=30}
: ${produce_pdf:=0}
: ${show:=10}

echo "# shell generated script to run predict_depth.r" > shell.predict_depth.r.$time_stamp.r

echo "source(\"~/predict_depth/predict_depth.watstats.r\")" >> shell.predict_depth.r.$time_stamp.r       

echo "predict_depth.watstats(abundance_matrix=\"$file_in\", col_num=$col_num, num_reads=$num_reads, input_type=\"file\", file_out_prefix = \"$file_out_prefix\", genome_size=$genome_size, read_length=$read_length, min_overlap=$min_overlap, scale_by_unannotated=$scale_by_unannotated, create_figure=$produce_pdf, num_to_show=$show)" >>  shell.predict_depth.r.$time_stamp.r

R --vanilla --slave <  shell.predict_depth.r.$time_stamp.r
rm  shell.predict_depth.r.$time_stamp.r

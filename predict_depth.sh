#!/bin/sh

#written 2-14-13
#if [ $# -eq 0 || $# -gt 3 ] # usage for no or more than 3 arguments                                                                                 
if [ $# -eq 0 -o $# -gt 7 ]
then
    echo "DESCRIPTION: Shell wrapper for predict_depth.r"
    echo 
    echo "USAGE: >predict_depth.sh <file_in> <col_num> <file_out_prefix> <genome_size> <coverage> <produce_pdf> <show>"
    echo   
    
    echo "     <file_in>:         NO DEFAULT; string, name of the input file" 
    echo "                        (two column, tab delimited, annotations then abundances" 
    echo "                        first row, first column (annotation header) is blank, first" 
    echo "                        row, scecond column can have a text header)"

    echo
    echo "     <col_num>:         DEFAULT = 1; 1 based index of column to process from file_in"
    echo


    echo
    echo "     <file_out_prefix>: DEFAULT = \"depth_prediction\"  ;string, prefix for output file(s)" 
    echo "                        *.txt (& *.pdf if produce_fig = 1)"
    echo
   
    echo 
    echo "     <genome_size>:     DEFAULT = 4000000; size in bp of genomes" 
    echo

    echo
    echo "     <coverage>:        DEFAULT = 30; Desired level of coverage"
    echo


    echo
    echo "     <produce_pdf>:     DEFAULT = FALSE;   TRUE|FALSE, produce a pdf visulization of the output"
    echo
    
    echo
    echo "     <show>:            DEFAULT = 10; Integer, number of taxa to show in pdf"
    echo
    
    echo "EXAMPLES: predict_depth.sh test_data.txt "
    echo "          predict_depth.sh test_data.txt 1 my_output 3000000 100 TRUE 20"
    echo
    echo "NOTES: Annoying -- to change argument 2, you have to supply 1 and 2 etc."
    echo "There are additional arguments in the function file (predict_depth.r)"
    echo
    exit 1                                                                                         
fi

time_stamp=`date +%m-%d-%y_%X:%N`;  # create the time stamp month-day-year_hour:min:sec:nanosec

# grab arguments from prompt
file_in=$1
col_num=$2
file_out_prefix=$3
genom_size=$4
coverage=$5
produce_pdf=$6
show=$7

# set default values
: ${col_num:=1}
: ${file_out_prefix:="depth_prediction"}
: ${genome_size:=4000000}
: ${coverage:=30}
: ${produce_pdf:=1}
: ${show:=10}



echo "# shell generated script to run predict_depth.r" > shell.predict_depth.r.$time_stamp.r

echo "source(\"~/predict_depth/predict_depth.r\")" >> shell.predict_depth.r.$time_stamp.r       

echo "predict_depth(abundance_matrix=\"$file_in\", col_num=$col_num, input_type=\"file\", file_out_prefix = \"$file_out_prefix\", genome_size=$genome_size, coverage=$coverage, num_to_show=$show, create_figure=$produce_pdf)" >>  shell.predict_depth.r.$time_stamp.r

R --vanilla --slave <  shell.predict_depth.r.$time_stamp.r
rm  shell.predict_depth.r.$time_stamp.r

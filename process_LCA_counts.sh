#!/bin/sh

#written 3-1-13
#if [ $# -eq 0 || $# -gt 4 ] # usage for no or more than 3 arguments                                                                                 
if [ $# -eq 0 -o $# -gt 4 ]
then
    echo "DESCRIPTION: Shell wrapper for process_LCA_counts.r"
    echo 
    echo "USAGE: > process_LCA_counts.sh <file_in> <include_ambiguous_counts> <tax_level> <relative_abundance>"
    echo   
    
    echo "     <file_in>:         NO DEFAULT; string, name of the input file" 
    echo "                        (two column, tab delimited, annotations then abundances" 
    echo "                        first row, first column (annotation header) is blank, first" 
    echo "                        row, scecond column can have a text header)"

    echo
    echo "     <include_ambiguous_counts>: DEFAULT = FALSE; "
    echo

    echo
    echo "     <tax_level>:       DEFAULT = \"genus\" ;string, tax level, one from the following:" 
    echo "                        \"domain\", \"phylum\", \"class\", \"order\",\"family\", \"genus\", \"species\""
    echo

    echo
    echo "     <relative_abundance>  DEFAULT = TRUE; logical, produce relative abundance output in addition to raw abundance output"
    echo
    
    echo "EXAMPLES:  process_LCA_counts.sh LCA_test_data.txt"
    echo "           process_LCA_counts.sh LCA_test_data.txt TRUE order TRUE"
    echo
    echo "NOTES: Annoying -- to change argument 2, you have to supply 1 and 2 etc."
    echo "There are additional arguments in the function file (process_LCA_counts.r)"
    echo
    exit 1                                                                                         
fi

time_stamp=`date +%m-%d-%y_%X:%N`;  # create the time stamp month-day-year_hour:min:sec:nanosec

# grab arguments from prompt
file_in=$1
include_ambiguous_counts=$2
tax_level=$3

# set default values
: ${col_num:=1}
: ${include_ambiguous_counts:=0}
: ${tax_level:="genus"}
: ${relative_abundance:=1}

echo "# shell generated script to run process_LCA_counts.r" > shell.process_LCA_counts.r.$time_stamp.r

echo "source(\"~/predict_depth/process_LCA_counts.r\")" >> shell.process_LCA_counts.r.$time_stamp.r       

echo "process_LCA_counts(abundance_matrix=\"$file_in\", include_ambiguous_counts=$include_ambiguous_counts, tax_level=\"$tax_level\", relative_abundance=$relative_abundance)" >>  shell.process_LCA_counts.r.$time_stamp.r

R --vanilla --slave <  shell.process_LCA_counts.r.$time_stamp.r
rm  shell.process_LCA_counts.r.$time_stamp.r

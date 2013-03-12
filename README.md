predict_depth
=============
Clone into ~



DESCRIPTION: r code and shell wrappers for Wizard tools
----------------------------------------------------------------------------------------------------\
--------
(process_LCA_counts.sh)
DESCRIPTION: Shell wrapper for process_LCA_counts.r

USAGE: > process_LCA_counts.sh <file_in> <include_ambiguous_counts> <tax_level> <relative_abundance>

     <file_in>:         NO DEFAULT; string, name of the input file
                        (two column, tab delimited, annotations then abundances
                        first row, first column (annotation header) is blank, first
                        row, scecond column can have a text header)

     <include_ambiguous_counts>: DEFAULT = FALSE;

     <tax_level>:       DEFAULT = "genus" ;string, tax level, one from the following:
                        "domain", "phylum", "class", "order","family", "genus", "species"

     <relative_abundance>  DEFAULT = TRUE; logical, produce relative abundance output in addition to raw abundance output


EXAMPLES:  process_LCA_counts.sh LCA_test_data.txt
           process_LCA_counts.sh LCA_test_data.txt TRUE order

NOTES: Annoying -- to change argument 2, you have to supply 1 and 2 etc.
There are additional arguments in the function file (process_LCA_counts.r)

***                                                                      ***
*** When would you use one version of merging the counts over the other? ***
*** (behavior defined by include_ambiguous_counts)                       ***

  If you are interested in a specific organism - you would only want to consider
  the counts that can be unambiguously associated with it, and would choose:
  include_ambiguous_counts = FALSE

  If you are interested in exploring a range of closely related organisms, you
  would want to consider all counts that can (unambiguously or ambiguously)
  associated with it, and would choose:
  include_ambiguous_counts = TRUE

------------------------------------------------------------------------------------------------------------
(predict_depth.sh)
DESCRIPTION: Shell wrapper for predict_depth.r

USAGE: >predict_depth.sh <file_in> <col_num> <file_out_prefix> <genome_size> <coverage> <scale_by_unannotated> <produce_pdf> <show>

     <file_in>:         NO DEFAULT; string, name of the input file
                        (two column, tab delimited, annotations then abundances
                        first row, first column (annotation header) is blank, first
                        row, scecond column can have a text header)

     <col_num>:         DEFAULT = 1; 1 based index of column to process from file_in

     <file_out_prefix>: DEFAULT = "depth_prediction"  ;string, prefix for output file(s)
                        *.txt (& *.pdf if produce_fig = 1)

     <genome_size>:     DEFAULT = 4000000; size in bp of genomes

     <coverage>:        DEFAULT = 30; Desired level of coverage
     
     <scale_by_unannotated>:   DEFAULT = TRUE; scale calculated depth including unannotated reads

     <read_length>:     DEFAULT = 125; average read length in bp

     <min_overlap>      DEFAULT = 30; min overlap for assembly

     <produce_pdf>:     DEFAULT = FALSE;   TRUE|FALSE, produce a pdf visulization of the output

     <show>:            DEFAULT = 10; Integer, number of taxa to show in pdf

EXAMPLES: predict_depth.sh test_data.txt 
          predict_depth.sh test_data.txt 1 my_output 4000000 25 125 30 TRUE TRUE 10

NOTES: Annoying -- to change argument 2, you have to supply 1 and 2 etc.
There are additional arguments in the function file (predict_depth.r)

Prediction values use calculations defined in the following:

    Lander, E. S. and Waterman, M. S., "Genomic mapping by fingerprinting
    random clones: a mathematical analysis", Genomics 2, 231-239 (1988).

------------------------------------------------------------------------------------------------------------

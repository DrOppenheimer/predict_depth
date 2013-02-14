predict_depth
=============
Clone into ~



DESCRIPTION: Shell wrapper for predict_depth.r

USAGE: >predict_depth.sh <file_in> <file_out_prefix> <produce_pdf> <show>

     <file_in>:         NO DEFAULT: string, name of the input file
                        (two column, tab delimited, annotations then abundances
                        first row, first column (annotation header) is blank, first
                        row, scecond column can have a text header)

     <file_out_prefix>: Default = "depth_prediction"  :string, prefix for output file(s)
                        *.txt (& *.pdf if produce_fig = 1)

     <produce_pdf>:     Default = TRUE:   TRUE|FALSE, produce a pdf visulization of the output

     <show>:            Default = 10: Integer, number of taxa to show in pdf

EXAMPLES: predict_depth.sh test_data.txt 
          predict_depth.sh test_data.txt "depth_prediction" TRUE 20

Annoying: to change argument 2, you have to supply 1 and 2 etc.
There are additional arguments in the function file (predict_depth.r)

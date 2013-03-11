process_LCA_counts<- function(
                              abundance_matrix,
                              include_ambiguous_counts = FALSE,
                              input_type = c("object","file"),
                              tax_level = c("domain", "phylum", "class", "order","family", "genus", "species"),
                              output_type = c("r.matrix", "file"),
                              relative_abundance = TRUE,
                              ambig_count_file = FALSE,
                              debug = FALSE
                              )
{

  # Print usage if no args are supplied
  if (nargs() == 0){ print_usage() }
  
  # supply default values if needed
  if(length(input_type)==2){
    input_type <- "file"
  }

  if(length(tax_level)==7){
    tax_level <- "genus"
  }

  if(length(output_type)==2){
    output_type <- "file"
  }

  # load packages
  require(hash)

  # name file out
  file_out <- gsub(" ", "", paste(abundance_matrix, ".LCA_processed_counts"))
 
  # get the index of the selected taxonomic level - return usge if argument is not valid
  if( length(tax_level) > 1){
    stop("you can only enter one tax level")
  }
  tax_index <<- grep(tax_level, c("domain", "phylum", "class", "order","family", "genus", "species"))
  if( tax_index == 0 ){ stop("you have to select a single tax level") }

  # load the data into a matrix
  if (input_type=="file"){
    my_data.matrix <- data.matrix(read.table(abundance_matrix, row.names=1, check.names=FALSE, header=TRUE, sep="\t", comment.char="", quote=""))
  }else{
    my_data.matrix <- data.matrix(abundance_matrix)
  }

  num_taxa_start <- dim(my_data.matrix)[[1]]
  #if(debug==TRUE){print(paste("Num start: ", num_taxa_start))}
  
  num_samples <-  dim(my_data.matrix)[[2]]
  tax_hash <- hash() # hash to fold counts that are defined (have numerical value) at the selected level
  tax_hash.ambig <- hash() # hash to contain values that are ambiguous  ("-") at the selected level
  
  for (j in 1:num_taxa_start){

    # get annotation string and counts
    line_annotation.character <- dimnames(my_data.matrix)[[1]][j]
    line_annotation.list <-  strsplit(line_annotation.character,";")
    line_counts.numeric <- as.numeric(my_data.matrix[j,])
    
    if ( identical( line_annotation.list[[1]][tax_index], "-") ){
    # process counts that are not defined at the selected level
      
      # create annotation string - is also used as the key for the counts
      for (m in 1:length(line_annotation.list[[1]])){
        if (m==1){
          if ( identical( line_annotation.list[[1]][m], "-") ){
          }else{
            line_annotation_ambig.character <- line_annotation.list[[1]][m]
          }     
        }else{
          if ( identical( line_annotation.list[[1]][m], "-") ){
          }else{
            line_annotation_ambig.character <- gsub(" ", "",paste(line_annotation_ambig.character, line_annotation.list[[1]][m], sep=";"))
          }
        }
      }

      # add counts to the hash for ambiguous counts
      if ( has.key(line_annotation_ambig.character, tax_hash.ambig)==TRUE ){
        tax_hash.ambig[ line_annotation_ambig.character ] <- tax_hash.ambig[[ line_annotation_ambig.character ]] + line_counts.numeric
      }else{
        tax_hash.ambig[ line_annotation_ambig.character ] <- line_counts.numeric
      }    
      
    }else{
    # process counts that are defined at the selected level

      # create annotation string - is also used as the key for the counts
      for (l in 1:tax_index){
        if (l==1){
          line_annotation_new.character <- line_annotation.list[[1]][l]
        }else{
          line_annotation_new.character <- gsub(" ", "",paste(line_annotation_new.character, line_annotation.list[[1]][l], sep=";"))
        }
      }
        
      # hash values lower order values into selected higher level order
      if ( has.key(line_annotation_new.character, tax_hash)==TRUE ){
        tax_hash[ line_annotation_new.character ] <- tax_hash[[ line_annotation_new.character ]] + line_counts.numeric
      }else{
        tax_hash[ line_annotation_new.character ] <- line_counts.numeric
      }
      
    }
    
  }
  
  # add ambiguous counts if that option is selected
  # note that the ambiguous counts are added to all taxa
  # that are a substring match (staring from beginning of strings)

  if( include_ambiguous_counts==TRUE ){
    #if(debug==TRUE){print("HELLO.2")}

    for (n in 1: length(keys(tax_hash.ambig)) ){
      ambig_key <<- keys(tax_hash.ambig)[n]

      for (o in 1:length(keys(tax_hash)) ){
        #if(debug==TRUE){print("HELLO.4")}
        tax_key <<- keys(tax_hash)[o]

        if( grepl(ambig_key,tax_key, fixed=TRUE)==TRUE ){
          tax_hash[ tax_key ] <- tax_hash[[ tax_key ]] + tax_hash.ambig[[ ambig_key ]]
        }
        
      }

    }

  }

  #if(debug){ tax_keys <<- keys(tax_hash); ambig_keys <<- keys(tax_hash.ambig)  }
  

  # place data in matrix for easy writing
  num_taxa <- length(keys(tax_hash))
  output.matrix <- matrix(0, num_taxa, num_samples)
  
  dimnames(output.matrix)[[2]] <- dimnames(my_data.matrix)[[2]]
  dimnames(output.matrix)[[1]] <- c(rep("",length(keys(tax_hash))))
  
  for (k in 1:length(keys(tax_hash))){
    dimnames(output.matrix)[[1]][k] <- keys(tax_hash)[k]
    output.matrix[k,] <- tax_hash[[ keys(tax_hash)[k] ]]
  }

  if ( output_type=="r.matrix"){
    return(output.matrix)
  }else{
    write.table(output.matrix, file = file_out, col.names=NA, row.names = TRUE, sep="\t", quote=FALSE)
  }

  # print a file that contains the relative abundances (for each metagenome), if the option is selection
  if( relative_abundance==TRUE ){
    
    relative.output.matrix <- matrix(0, num_taxa, num_samples)
    dimnames(relative.output.matrix)[[2]] <- dimnames(my_data.matrix)[[2]]
    dimnames(relative.output.matrix)[[1]] <- c(rep("",length(keys(tax_hash))))

    for (k in 1:length(keys(tax_hash))){
      dimnames(output.matrix)[[1]][k] <- keys(tax_hash)[k]
      output.matrix[k,] <- tax_hash[[ keys(tax_hash)[k] ]]
    }
    
    for( my_col in 1:dim(output.matrix)[2] ) {
      for( my_row in 1:dim(output.matrix)[1] ) {    
        relative.output.matrix[my_row,my_col] = ((( output.matrix[my_row,my_col] - min(output.matrix[,my_col]))*(1))/(max(output.matrix[,my_col]) - min(output.matrix[,my_col])))
        
        #NewValue = (((OldValue - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
        #Or a little more readable:
        #OldRange = (OldMax - OldMin)  
        #NewRange = (NewMax - NewMin)  
        #NewValue = (((OldValue - OldMin) * NewRange) / OldRange) + NewMin
        
      }
    }
    
    write.table(relative.output.matrix, file =  (paste(file_out,".scaled",sep="")), col.names=NA, row.names = TRUE, sep="\t", quote=FALSE)
    
  }
  
  # print file that has the ambiguous counts, if option is selected
  if( ambig_count_file==TRUE ){

    num_taxa.ambig <- length(keys(tax_hash.ambig))
    output.matrix.ambig <- matrix(0, num_taxa.ambig, num_samples)
  
    dimnames(output.matrix.ambig)[[2]] <- dimnames(my_data.matrix)[[2]]
    dimnames(output.matrix.ambig)[[1]] <- c(rep("",length(keys(tax_hash.ambig))))
  
  for (z in 1:length(keys(tax_hash.ambig))){
    dimnames(output.matrix.ambig)[[1]][z] <- keys(tax_hash.ambig)[z]
    output.matrix.ambig[z,] <- tax_hash.ambig[[ keys(tax_hash.ambig)[z] ]]
  }

    write.table(output.matrix.ambig, file = paste(file_out,".ambig", sep=""), col.names=NA, row.names = TRUE, sep="\t", quote=FALSE)

  }

 
}
    

   print_usage <- function() {
  writeLines(
 "  ------------------------------------------------------------------------------
  process_LCA_counts.r                     Kevin P. Keegan, kkeegan@anl.gov  Feb 2013
  ------------------------------------------------------------------------------
  DESCRIPTION:
  Script to process LCA counts for predict_depth.r

  USAGE:
  process_LCA_counts(abundance_matrix,
                              input_type = c(\"object\",\"file\"),   # default = file
                              include_ambiguous_counts = TRUE|FALSE, # default = FALSE
                              tax_level = c(\"domain\", \"phylum\", \"class\", \"order\",\"family\", \"species\"), # default=\"genus\"
                              output_type = c(\"r.matrix\", \"file\") #default = \"r.matrix\",
                              ambig_count_file = TRUE|FALSE, # default = FALSE
                              debug = TRUE|FALSE # default = FALSE

  NOTES:
  This script can process a table of LCA or Best hit counts to generate counts
  needed for predict_depth.  Counts are merged in one of two ways, LCA counts
  at higher level are either left out ( include_ambiguous_counts = FALSE ) or
  included ( include_ambiguous_counts = FALSE ).  If they are included, counts
  higher order counts are added to each lower order that has them as a subtring.
  As an example:

  If genus is the selected level, counts from the following two higher order LCAs: 

       order;family;genus

       domain;phylum;class;order;family

  Would be added to each genus they match.

  *** When would you use one version of merging the counts over the other? ***

  If you are interested in a specific organism - you would only want to consider
  the counts that can be unambiguously associated with it, and would choose:
  include_ambiguous_counts = FALSE

  If you are interested in exploring a range of closely related organisms, you
  would want to consider all counts that can (unambiguously or ambiguously)
  associated with it, and would choose:
  include_ambiguous_counts = TRUE

  ------------------------------------------------------------------------------"
             )
  stop("you did not enter the correct args -- see above")
}

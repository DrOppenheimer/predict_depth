process_LCA_counts<- function(
                              abundance_matrix,
                              input_type = c("object","file"),
                              tax_level = c("domain", "phylum", "class", "order","family", "genus", "species"),
                              output_type = c("r.matrix", "file"),
                              debug = TRUE
                              )
{

  # Print usage if no args are supplied
  if (nargs() == 0){ print_usage() }
  
  # load packages
  require(hash)

  if(debug==TRUE){print(paste("abundance_matrix: ", abundance_matrix))}
  
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
    my_data.matrix <<- data.matrix(read.table(abundance_matrix, row.names=1, check.names=FALSE, header=TRUE, sep="\t", comment.char="", quote=""))
  }else{
    my_data.matrix <<- data.matrix(abundance_matrix)
  }

  num_taxa_start <- dim(my_data.matrix)[[1]]
  if(debug==TRUE){print(paste("Num start: ", num_taxa_start))}
  
  num_samples <-  dim(my_data.matrix)[[2]]
  tax_hash <<- hash()

  #num_taxa_end <- length(grep("_", my_data.matrix[,tax_index]))

  #if(debug==TRUE){print(paste("Num end:   ", num_taxa_end))}
  
  #output.matrix <- matrix(NA, num_taxa.end, num_samples)
  
  for (j in 1:num_taxa_start){

    if(debug==TRUE){print("HELLO.1")}
    if(debug==TRUE){print(paste("j: ", j))}
    # replace taxa string with string that ends at selected level - new tax string is used as hash key to combine counts
    line_annotation.character <- dimnames(my_data.matrix)[[1]][j]
    ###if(debug==TRUE){print(paste("OG annotation: ", line_annotation.character))}

    line_annotation.list <<-  strsplit(line_annotation.character,";")

    #if(debug==TRUE){ print(paste("- check", grep( "-",line_annotation.list[[1]][tax_index] ))) }
    if(debug==TRUE){ print(paste("annotation entry [",tax_index,"] :", line_annotation.list[[1]][tax_index] )) }
    ###if(debug==TRUE){print("GOT HERE")}

    if ( identical( line_annotation.list[[1]][tax_index], "-") ){

      if(debug==TRUE){print("HELLO.1-5")}
    # drop line if counts at selected level are "_" -- i.e. are not defined
    # if you want to propegate counts down - this is where you'd start  
    }else{

      if(debug==TRUE){print("HELLO.2")}
      # get the counts for the current taxa(line)
      #line_counts.character <<- as.character(my_data.matrix[j,])
      line_counts.numeric <- as.numeric(my_data.matrix[j,])
      if(debug==TRUE){print(paste("Line counts: ", line_counts.numeric))}
      
      for (l in 1:tax_index){


        if(debug==TRUE){print("HELLO.3")}
        if (l==1){

          if(debug==TRUE){print("HELLO.4")}
          line_annotation_new.character <- line_annotation.list[[1]][l]
        }else{

          if(debug==TRUE){print("HELLO.5")}
          line_annotation_new.character <- gsub(" ", "",paste(line_annotation_new.character, line_annotation.list[[1]][l], sep=";"))
        }

      }
      
      if (debug==TRUE){print(paste("line_annotation_list         : ", line_annotation.list))}
        
    # hash values lower order values into selected higher level order

      if (debug==TRUE){print(paste("line_annotation_new.character: ",line_annotation_new.character))}
      
      if ( has.key(line_annotation_new.character, tax_hash)==TRUE ){

        if(debug==TRUE){print("HELLO.6")}
        tax_hash[ line_annotation_new.character ] <<- tax_hash[[ line_annotation_new.character ]] + line_counts.numeric
      }else{

        if(debug==TRUE){print("HELLO.7")}
        tax_hash[ line_annotation_new.character ] <<- line_counts.numeric
      }

      if(debug==TRUE){print("HELLO.8")}
      
    }

  }
  
  if(debug==TRUE){print("HELLO.9")}
  # place data in matrix for easy writing
  num_taxa <<- length(keys(tax_hash))
  output.matrix <<- matrix(0, num_taxa, num_samples)
  
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

  #}

}
    

    
    ## for (k in 1:length(keys(my_hash))){
    ##   print(paste("key  : ",keys(my_hash)[k]))
    ##   print(paste("value: ",my_hash[[ keys(my_hash)[k] ]]))
    ## }

    ## # This works
    ## for (k in 1:length(keys(my_hash))){
    ##   print(paste("key  : ",keys(my_hash)[k]))
    ##   print(paste("value: ",my_hash[[ keys(my_hash)[k] ]]))
    ## }

    
##     writeLines(  )
    
##  keys(my_hash)
    

##         tax_hash[  ]
    
## sum(as.numeric(my_test[1:2]))

##   }
    
##   }

## my_hash[ test ] <- (my_hash[[ test ]] + 3)

## my_hash[ test ] <- 3
  

## has.key signature(key = "ANY", hash = "hash"): Test for existence of key

  
## value , <-  my_hash[[key]]
  
## Examples
## h <- hash( keys=letters, values=1:26 )
## h <- hash( letters, 1:26 )
## h$a # 1
## h$foo <- "bar"
## h[ "foo" ]
## h[[ "foo" ]]
## clear(h)
## rm(h)

  
  
  
  
## # convert tax level to int
## if ( ! identical(class(tax_level), "numeric" ){
  


## }


   
# Check to see if abundance_matrix is a file - if so, load it - if not, assume it is an R object and copy it
 


   print_usage <- function() {
  writeLines(
 "  ------------------------------------------------------------------------------
  process_LCA_counts.r                     Kevin P. Keegan, kkeegan@anl.gov  Feb 2013
  ------------------------------------------------------------------------------
  DESCRIPTION:
  Script to process LCA counts for predict_depth.r

  USAGE:
  process_LCA_counts(abundance_matrix,
                              input_type = c(\"object\",\"file\"),
                              tax_level = c(\"domain\", \"phylum\", \"class\", \"order\",\"family\", \"species\"),
                              output_type = c(\"r.matrix\", \"file\")

  NOTES:

  ------------------------------------------------------------------------------"
             )
  stop("you did not enter the correct args -- see above")
}

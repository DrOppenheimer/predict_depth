dist_hist <- function(distance_matrix_file, groups_list = "EHFI_groups", debug = 0){

POO
# dist_hist <- function(distance_matrix_file, file_path = "./", groups_list = "groups", perm_dists_path = "./DISTs", debug = 0){ 

  legend_text <- c("within_group","between_group")
  legend_colors <- c("red", "green")

# First, deal with the original file

  distance_matrix <<- data.matrix(read.table(distance_matrix_file, row.names=1, header=TRUE, check.names=FALSE, sep="\t", comment.char="", quote=""))
  
  groups_dataframe <<- read.table(groups_list, header=FALSE, check.names=FALSE, sep = ",", comment.char="", quote="", fill=TRUE, blank.lines.skip=TRUE)
  #groups_characterlist <<- readLines(groups_list)

  groups_in <<- readLines(groups_list)
  num_groups <<- dim(data.frame(groups_in))[1]

  within_group_distances <- matrix()
  between_group_distances <- matrix()
  
  for (i in 1:num_groups){
    
    if ( debug>0 ) { print(paste("Group:", i)) }
    
    group_samples <- strsplit(groups_in[i], ",")
    num_samples <- dim(data.frame(group_samples))[1]
    group_distances <<- matrix (NA, size(combn(num_samples,2))[2], 1) #####
    gd_index <- 1 
    
    for (m in 2:num_samples){ # get all of the unique non-redundant pairs
      
      for (n in 1:(m-1)){
        
        if (
            identical( toString(charmatch( group_samples[[1]][m], dimnames(distance_matrix)[[1]])), "NA" ) ||
            identical( toString(charmatch(   group_samples[[1]][n], dimnames(distance_matrix)[[1]])), "NA" )
            )
          {
            
            stop (paste("sample names in the groups file:\n", groups_list,"\n (", noquote(group_samples[[1]][m]), " or ", noquote(group_samples[[1]][n]), ")\n are not in the distance_matrix:\n",distance_matrix_file, "\n",
                        "check the names in the groups and distance_matrix files to make sure that they match\n\n"))    
          }
        
        if (debug>0){ print(paste("gd_index =", gd_index, "distance:", distance_matrix[ noquote(group_samples[[1]][m]), noquote(group_samples[[1]][n]) ] )) } #####
        group_distances[gd_index,1] <<- distance_matrix[ noquote(group_samples[[1]][m]), noquote(group_samples[[1]][n]) ] #####
        gd_index <- gd_index + 1
        
      }
      
    }
    
    if ( i == 1 ){
      within_group_distances <- group_distances
    }else{
      within_group_distances <- rbind(within_group_distances, group_distances)
    }
    
  }

  
  #  
  for (p in 2:num_groups){ # get all of the unique non-redundant pairs
    
    for (q in 1:(p-1)){
      
      alpha_samples <- strsplit(groups_in[p], ",")
      num_alpha_samples <- dim(data.frame(alpha_samples))[1]
                                        #alpha_distances <- matrix(NA, size(combn(num_alpha_samples,2))[2], 1) ### ###  
      
      beta_samples <- strsplit(groups_in[q], ",")
      num_beta_samples <- dim(data.frame(beta_samples))[1]
                                        #beta_distances <- matrix(NA, size(combn(num_beta_samples,2))[2], 1) ### ###
      
      alpha_beta_distances <<- matrix(NA, (num_alpha_samples*num_beta_samples), 1) # alpha and beta should be the same size
      dist_index <- 1 ### ###
      
      write(gsub(" ", "", (paste(">>Group(", p, ")::Group(", q, ")"))), file=output_file, append=TRUE)
      
      for (na in 1:num_alpha_samples){
        for (nb in 1:num_beta_samples){
          
          alpha_beta_distances[dist_index,1] <<- distance_matrix[ noquote(alpha_samples[[1]][na]), noquote(beta_samples[[1]][nb]) ] ### ###
          dist_index <- dist_index + 1 ### ###
      
        }
      }
      
    }

    if( p == 1 ){
      between_group_distances <- alpha_beta_distances
    }else{
      between_group_distances <- rbind(between_group_distances, alpha_beta_distances)
    }
    
    
  }   
  
  
within_group.hist <-  hist(within_group_distances, plot=FALSE, breaks=20)
between_group.hist <-  hist(between_group_distances, plot=FALSE, breaks=20)

plot(x = within_group.hist$breaks, y = c(within_group.hist$counts, 0), ylab="frequency", xlab="breaks", type="l", col=legend_colors[1])
lines(x = between_group.hist$breaks, y = c(between_group.hist$counts, 0), ylab="frequency", xlab="breaks", type="l", col=legend_colors[2])
legend("topright", legend = legend_text, col = legend_colors)
  
}































  
  








dist_hist <- function(original_dist,
                          #input_type = "file",
                          #original_dist,
                          perm_dists_path = "./DISTs",
                          file_out_prefix = "my_dist_hist",
                          create_figure=FALSE){

# Print usage if
   #if (nargs() == 0){print_usage()}
   
# mport the original dist file
original.dist <- data.matrix(read.table(original_dist, row.names=1, check.names=FALSE, header=TRUE, sep="\t", comment.char="", quote=""))

original_hist <- hist(as.dist(original.dist), plot=FALSE, breaks=20)

sum.dist <- matrix()

perm_dists <- dir(perm_dists_path)

for (i in 1:length(perm_dists)){

  perm.dist <- data.matrix(read.table(  paste( perm_dists_path,"/",perm_dists[i], sep="") , row.names=1, check.names=FALSE, header=TRUE, sep="\t", comment.char="", quote=""))

  # perm_hist <- hist(perm.dist, plot=FALSE, breaks=20)

  if ( i==1 ){
    sum.dist = perm.dist
  }else{
    sum.dist = sum.dist + perm.dist
  }
 
}

avg.perm.dist = ( sum.dist / length(perm_dists) )

#perm_hist <- hist(avg.perm.dist, plot=FALSE, breaks=20)
perm_hist <- hist(as.dist(avg.perm.dist), plot=FALSE, breaks=20)

print(attributes(as.dist(original.dist)))
print(attributes(as.dist(avg.perm.dist)))

  #num_breaks <- length(my_hist$breaks)
  
  #print(paste(my_hist$counts, 0))
  


  

  #print(my_hist$breaks)
  #my_hist <- hist(temp.matrix)

}
























## dist_hist <- function(original_dist,
##                           #input_type = "file",
##                           #original_dist,
##                           perm_dists_path = "./DISTs",
##                           file_out_prefix = "my_dist_hist",
##                           create_figure=FALSE){

## # Print usage if
##    #if (nargs() == 0){print_usage()}
   
## # mport the original dist file
## original.dist <- data.matrix(read.table(original_dist, row.names=1, check.names=FALSE, header=TRUE, sep="\t", comment.char="", quote=""))

## original_hist <- hist(as.dist(original.dist), plot=FALSE, breaks=20)

## sum.dist <- matrix()

## perm_dists <- dir(perm_dists_path)

## for (i in 1:length(perm_dists)){

##   perm.dist <- data.matrix(read.table(  paste( perm_dists_path,"/",perm_dists[i], sep="") , row.names=1, check.names=FALSE, header=TRUE, sep="\t", comment.char="", quote=""))

##   # perm_hist <- hist(perm.dist, plot=FALSE, breaks=20)

##   if ( i==1 ){
##     sum.dist = perm.dist
##   }else{
##     sum.dist = sum.dist + perm.dist
##   }
 
## }

## avg.perm.dist = ( sum.dist / length(perm_dists) )

## #perm_hist <- hist(avg.perm.dist, plot=FALSE, breaks=20)
## perm_hist <- hist(as.dist(avg.perm.dist), plot=FALSE, breaks=20)

## print(attributes(as.dist(original.dist)))
## print(attributes(as.dist(avg.perm.dist)))

##   #num_breaks <- length(my_hist$breaks)
  
##   #print(paste(my_hist$counts, 0))
  
## plot(x = original_hist$breaks, y = c(original_hist$counts, 0), ylab="frequency", xlab="breaks", type="l", col="blue")
## lines(x = perm_hist$breaks, y = c(perm_hist$counts, 0), ylab="frequency", xlab="breaks", type="l", col="red")


  

##   #print(my_hist$breaks)
##   #my_hist <- hist(temp.matrix)

## }
  

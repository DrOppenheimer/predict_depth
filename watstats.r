watstats <- function (
                      num_reads = 1000000,
                      percent_data = 20, # use to get num reads
                      genome_length = 4000000,
                      read_length = 125,
                      min_overlap = 30,
                      verbose=TRUE
                      )
  {

     taxa_num_reads = ( num_reads * ( percent_data/100 ) ) # 
     
     alpha    <<- ( taxa_num_reads/genome_length ) # $alpha=$N/$GM; # $GM = $G*1000; (input in original was in KB)
     theta    <<- ( min_overlap/read_length )      # $theta=$T/$L;
     sigma    <<- ( 1-theta )                      # $sigma=1-$theta;

     coverage_redundancy <<- ( (read_length*taxa_num_reads)/genome_length )          # $c=$L*$N/$GM;

     num_contigs      <<- taxa_num_reads*exp(-coverage_redundancy*sigma)    # $i  =$N*exp(-$c*$sigma); 
     if ( num_contigs < 1 ){ num_contigs <<- 1 }                            # $i=1    if $i   < 1;
     if ( num_contigs > num_reads ){ num_contigs <<- num_reads }            # $i  =$N if $i   > $N;
     
     seqs_per_contig  <<- exp(coverage_redundancy*sigma)                             # exp($c*$sigma);
     if ( seqs_per_contig > num_reads ){ seqs_per_contig <<- num_reads }             # $iii=$N if $iii > $N;
     if ( seqs_per_contig < 1 ){ seqs_per_contig <<- 1 }                             # $iii=1  if $iii < 1;

     contig_length    <<- read_length*(((exp(coverage_redundancy*sigma)-1)/coverage_redundancy)+(1-sigma)) # $iv=int($L*(((exp($c*$sigma)-1)/$c)+(1-$sigma)));
     if ( contig_length > genome_length ){ contig_length <<- genome_length }         # $iv=$GM if $iv  > $GM;
  
     percent_coverage <<- 100*num_contigs*contig_length/genome_length                # $compl=int(100*$i*$iv/$GM);
     if ( percent_coverage > 100 ){ percent_coverage <- 100 }
     
     if( verbose==TRUE ){
       print("INPUT")
       print(paste("num_reads           :", round(taxa_num_reads, digits=0)))
       print(paste("percent_data        :", round(percent_data, digits=1)))
       print(paste("genome_length       :", round(genome_length, digits=0)))
       print(paste("read_length         :", round(read_length,digits=0)))
       print(paste("min_overlap         :", round(min_overlap, digits=0)))
       print("")
       print("OUTPUT")
       print(paste("coverage_redundancy :", round(coverage_redundancy, digits=1)))
       print(paste("num_contigs         :", round(num_contigs, digits=1)))
       print(paste("seqs_per_contig     :", round(seqs_per_contig, digits=1)))
       print(paste("contig_length       :", round(contig_length, digits=1)))
       print(paste("percent_coverage    :", round(percent_coverage, digits=1), "%"))
       print('------------------------------------------------------')
     }
     
     watstats_command <- paste("watstats.pl -g", genome_length, "-n", taxa_num_reads, "-l", read_length, "-t",  min_overlap, "intern=TRUE")
     
     my_watstats <<- system(watstats_command)

    #watstats_command <- paste("watstats.pl -g", genome_size, "-n", taxa_num_reads, "-l", read_length, "-t",  min_overlap, "intern=TRUE")

    # system2(watstats_command, stdout=my_watstats.2)
    
    #my_watstats.2 <<- capture.output(system(watstats_command), file = NULL, append = FALSE)
    
    my_date <<- system("date", intern=TRUE)


  }
# watstats.pl -g 5000 -n 1000 -l 480 -t 30

# watstats 1.1 ; Based on
#    Lander, E. S. and Waterman, M. S., "Genomic mapping by fingerprinting
#    random clones: a mathematical analysis", Genomics 2, 231-239 (1988).

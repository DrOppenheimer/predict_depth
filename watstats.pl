#!/usr/bin/env  perl
# compute some statistics for shotgun projects
# $Id: watstats,v 1.2 2002/12/06 11:21:49 fm Exp fm $
# Programming by Pierre Rioux, OGMP, 1994
# Fixing by Folker Meyer, Bielefeld University, 2001

use Getopt::Std;
              
################################################################
# read arguments
getopts('hvg:n:l:t:');

if (defined $opt_h) {
	&usage;
	die "\n";;
}

if (defined $opt_v){
 $debug=1;
}

# number of reads (default 500)
if (defined $opt_n){
  if ($debug) {
     print "$0 Using N:$opt_n \n";
  }
  $N=$opt_n;
}
else {
 $N=500;
}

# genome size in kb (default = 70)
if (defined $opt_g){
  if ($debug) {
     print "$0 Using G:$opt_g \n";
  }
  $G=$opt_g;
}
else {
 $G=70
}

# avg read length (default = 480)
if (defined $opt_l){
  if ($debug) {
     print "$0 Using L:$opt_l \n";
  }
  $L=$opt_l;
}
else {
 $L=480
}

# overlap length (default = 30)
if (defined $opt_t){
  if ($debug) {
     print "$0 Using T:$opt_t \n";
  }
  $T=$opt_t;
}
else {
 $T=30
}
# parsing done
##########################################

$O=200;

$TABLE=0;
$OCEAN=0;

    $VV="";
    foreach $n ("G","L","N","T","O") {
        eval "(\$x,\$y,\$z) = (\$$n =~ /^\\D*(\\d+)\\D+(\\d+)\\D+(\\d+)/);";
        next if $x eq "";
        if ($VV) {
            print "Error: you can't make more than one variable parameter.\n";
            $VV="ERROR";
            last;
            }
        ($from,$to,$step) = ($x,$y,$z);
        $VV=$n;
        eval "\$$VV=$from;";
        $TABLE=1;
        }
    next if ($VV eq "ERROR");

    if (!$VV) {
        $VV="G";
        $from=$G;$to=$G;$step=1;
        $TABLE=0;
        }

    next if ($G <= 0 || $T <= 0 || $L <= 0 || $N <= 0 || $O <= 0);

    if ($TABLE) {
        eval "\$$VV=\"From $from to $to, step of $step\"";
        }

    print "Length of genome            G = ${G}000 bp\n";
    print "Number of sequences         N = $N\n";
    print "Average length of sequences L = $L bp\n";
    print "Minimum overlap             T = $T bp\n";
    print "Ocean threshold             O = $O bp\n" if $OCEAN;

    if ($TABLE) {
        print "\n";
        print "$VV\tCover\t#Contig\t#Seq/C\tContLen\t%Compl";
        print "\tOcean>$O" if $OCEAN;
        print "\n-\t-----\t-------\t------\t-------\t------";
        print "\t---------" if $OCEAN;
        print "\n";
        }

    for ($range=$from;$range<=$to;$range+=$step) {
        eval "\$$VV=\$range;";

        $GM = $G*1000;
        $alpha=$N/$GM;
        $theta=$T/$L;
        $sigma=1-$theta;
        $c    =$L*$N/$GM;

        #$sigma=1;  # Override: means overlap ALWAYS detected.

        $cover=$c;
        $i  =$N*exp(-$c*$sigma);                            $i=1    if $i   < 1; $i  =$N if $i   > $N;
        $iii=exp($c*$sigma);                                $iii=$N if $iii > $N;$iii=1  if $iii < 1;
        $iv=int($L*(((exp($c*$sigma)-1)/$c)+(1-$sigma)));   $iv=$GM if $iv  > $GM;
        $vi=int(100*exp(-$c*($O/$L+$theta)));
        $compl=int(100*$i*$iv/$GM);

        &AdjustValue(*i);
        &AdjustValue(*iii);
        #&AdjustValue(*iv);
        &AdjustValue(*cover);

        if (!$TABLE) {
            print "\n";
            print "Redundancy of coverage                   = $cover\n";
            print "Expected number of contigs               = $i\n";
            print "The expected number of seqs in a contig  = $iii\n";
            print "The expected length of each contig       = $iv\n";
            print "Percentage of coverage of genome         = $compl%\n";
            printf"The prob of having an ocean > %s bp     = %-3.1f%%\n",$O,$vi
                if $OCEAN;
            }
        else {
            if ($OCEAN) {
                print "$range\t$cover\t$i\t$iii\t$iv\t$compl%\t";printf "%-3.1f\n",$vi;
                }
            else {
                print "$range\t$cover\t$i\t$iii\t$iv\t$compl%\n";
                }
            }
        }
        if ($TABLE) {
            eval "\$$VV=\"From $from to $to, step of $step\"";
            }


#########################################################
sub usage{
print <<BANNER;
USAGE: watstats -g GENOMESIZE in KB 
                -n # of reads (either number or from-to,step)
                -l avg. read length
                -t min. overlap needed

example:  watstats -g 5000 -n 1000-70000,500 -l 480 -t 30

watstats 1.1 ; Based on

    Lander, E. S. and Waterman, M. S., "Genomic mapping by fingerprinting
    random clones: a mathematical analysis", Genomics 2, 231-239 (1988).

Programming by Pierre Rioux, OGMP, 1994
Fixing by Folker Meyer, ZfGB, 2001
BANNER
}
##################################
sub AdjustValue {
    local(*VarName) = @_;
    $VarName=sprintf("%4.1f",$VarName);
}

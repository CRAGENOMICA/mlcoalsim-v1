#!/usr/local/bin/perl

use Getopt::Long;
use File::Find;
use strict;

my (%options,$flag,$flag2,$flag3,@headers,$line,$percent,$i,$j,$k,$count,$rowp,$rowm,$value,$floor,$filenameout,$filename,@files);
my @null_intervals;
my (%parameter_null,%parameter_altern);

%parameter_null = (
    N_iterations=>"",
    Seed1=>"",
    N_loci=>"",
    Sfix_allthetas=>"",
    Range_thetant=>"",
    N_samples=>"",
    N_sites=>"",
    Recombination=>"",
    Mutations=>"",
    Thetaw=>"",
    Nintn=>"",
    Nrec=>"",
    Npast=>"",
    Tpast=>"",
    Npop=>"",
    Mig_rate=>"",
    Ifselection=>"",
    Pop_size=>"",
    Pop_sel=>"",
    Sel_nt=>"",
    Sinit=>""
);
 %parameter_altern = %parameter_null;

#options
GetOptions(
    'obs=s'  => \$options{'obs'}, #path to the file with the observed values (MANDATORY: header EQUAL to simulated file: data separated by tabs, include blanks)
    'sim=s' => \$options{'sim'},    #path to the file with mlcoal simulations
);

#Usage
if($options{'obs'} eq undef || $options{'sim'} eq undef) {
   print "\nUsage: \nperl probstats.pl -obs (file with observed values) -sim (file with simulated values)\n\n";
   exit;
}

#start with the null model
open(INPUTNULL,"$options{'obs'}") or die;

#reading the file
$flag=0;
while(<INPUTNULL>){
    #if(($_ =~/[A-Z]/ && $_!~/\tna/) || ($_=~/[a-z]/ && $_!~/\tna/)) 
    if($_=~/(TD)/ || $_=~/\[0\]/)
        {$flag=1;}
    if($flag==2) {
        @null_intervals = split(/\t/,$_);
    }
    if($flag==1) {
        @headers = split(/\t/,$_);
        $flag=2;
    }
}
close(INPUTNULL);

#set na, nan or "" to 1234567890
for($i=0;$i<=$#headers;$i++) {
    if($null_intervals[$i] eq "na" || $null_intervals[$i] eq "nan" || $null_intervals[$i] eq "") {
         $null_intervals[$i] = 1234567890;
    }
}


#print this results in an output file
#open the output file
($options{'obs'})=~/(\w+)\./;
$filenameout = ">"."$1"."_probstats".".txt";
open(OUTPUT,$filenameout) or die;
#print quantiles for std model
print OUTPUT "Probability to have a value in the simulated model smaller than the observed value\n\n";

print OUTPUT "Name of the observed data file: $options{'obs'}\n\n";
print OUTPUT "Observed data values:\n";
for($i=0;$i<=$#headers;$i++) {
    print OUTPUT "\t$headers[$i]";
}
for($i=0;$i<=$#headers;$i++) {
    if($null_intervals[$i] != 1234567890) {print OUTPUT "\t$null_intervals[$i]";}
    else {print OUTPUT "\tna";}
}
print OUTPUT "\n";
#then calculate the probability given the simulated model file
#print label
print OUTPUT "\nfile\t";
print OUTPUT "N_iterations\t";
print OUTPUT "Seed1\t";
print OUTPUT "N_loci\t";
print OUTPUT "Sfix_allthetas\t";
print OUTPUT "Range_thetant\t";
print OUTPUT "N_samples\t";
print OUTPUT "N_sites\t";
print OUTPUT "Recombination\t";
print OUTPUT "Mutations\t";
print OUTPUT "Thetaw\t";
print OUTPUT "Nintn\t";
print OUTPUT "Nrec\t";
print OUTPUT "Npast\t";
print OUTPUT "Tpast\t";
print OUTPUT "Npop\t";
print OUTPUT "Mig_rate\t";
print OUTPUT "Ifselection\t";
print OUTPUT "Pop_size\t";
print OUTPUT "Pop_sel\t";
print OUTPUT "Sel_nt\t";
print OUTPUT "Sinit\t";

print OUTPUT "Probability";

for($i=0;$i<=$#headers;$i++) {
    print OUTPUT "\t$headers[$i]";
}

#read simulated file
open(INPUTALTERN,"$options{'sim'}") or die;
#introduce data in a matrix
my @matrix_altern;
my @headers_altern;
$flag=0;
$flag2=0;
while(<INPUTALTERN>){
    if($_=~/(TD)/ || $_=~/\[0\]/) 
        {$flag=1;}
    if($flag==2) {
        push(@matrix_altern,[split(/\t/,$_)]);
    }
    if($flag==1) {
        @headers_altern = split(/\t/,$_);
        $flag=2;
    }
    if($flag==0) {
        #I also want to take the parameters values...
        if($_=~/N_iterations:\s+(\d+)/) {$parameter_altern{N_iterations} = $1;}
        if($_=~/Seed1:\s+(\d+)/) {$parameter_altern{Seed1} = $1;}

        if($_=~/N_loci:\s+(\d+)/) {$parameter_altern{N_loci} = $1;}
        if($_=~/Sfix_allthetas:\s+(\d+)/) {$parameter_altern{Sfix_allthetas} = $1;}
        if($_=~/Range_thetant:\s+(\d+)/) {$parameter_altern{Range_thetant} = $1;}
        
        if($_=~/N_samples:\s+(.+)/) {$parameter_altern{N_samples} = $1;}
        if($_=~/N_sites:\s+(.+)/) {$parameter_altern{N_sites} = $1;}
        if($_=~/Recombination:\s+(.+)/) {$parameter_altern{Recombination} = $1;}
        if($_=~/Mutations:\s+(.+)/) {$parameter_altern{Mutations} = $1;}
        if($_=~/Thetaw:\s+(.+)/) {$parameter_altern{Thetaw} = $1;}
        
        if($_=~/Nintn:\s+(\d+)/) {$parameter_altern{Nintn} = $1;}
        if($_=~/Nrec:\s+(.+)/) {$parameter_altern{Nrec} = $1;}
        if($_=~/Npast:\s+(.+)/) {$parameter_altern{Npast} = $1;}
        if($_=~/Tpast:\s+(.+)/) {$parameter_altern{Tpast} = $1;}
        
        if($_=~/Npop:\s+(\d+)/) {$parameter_altern{Npop} = $1;}
        if($_=~/Mig_rate:\s+(\d+)/) {$parameter_null{Mig_rate} = $1;}
        
        if($_=~/Ifselection:\s+(.+)/) {$parameter_altern{Ifselection} = $1;}
        if($_=~/Pop_size:\s+(\d+)/) {$parameter_altern{Pop_size} = $1;}
        if($_=~/Pop_sel:\s+(.+)/) {$parameter_altern{Pop_sel} = $1;}
        if($_=~/Sel_nt:\s+(.+)/) {$parameter_altern{Sel_nt} = $1;}
        if($_=~/Sinit:\s+(.+)/) {$parameter_altern{Sinit} = $1;}
    }
}
close(INPUTALTERN);
#check for headers
if($#headers != $#headers_altern) {
    print "\nError: Not the same headers in obs and simulated values.\n";
}
else {
    for($i=0;$i<=$#headers;$i++) {
        if($headers[$i] ne $headers_altern[$i]) {
          print "\nError: Not the same headers in obs and simulated values.\n";
          die;
        }
    }
    #for each column
    my @alternative_probp;
    my @alternative_probm;
    for($i=0;$i<=$#headers_altern;$i++) {
        #count the number of values and how many values are below and above the quantiles
        $count = 0;
        $rowp = 0;
        $rowm = 0;
        for($j=0;$j<=$#matrix_altern;$j++) {
            if($matrix_altern[$j][$i] ne "na" && $matrix_altern[$j][$i] ne "nan" && $matrix_altern[$j][$i] ne "" && $null_intervals[$i] != 1234567890) {
                $count++;
                if($matrix_altern[$j][$i] <= $null_intervals[$i]) {$rowp = $rowp + 1;}
                if($matrix_altern[$j][$i] >= $null_intervals[$i]) {$rowm = $rowm + 1;}
            }
        }
        #calculate proportion
        if($count>0) {
            $rowp = $rowp/$count;
            $rowm = $rowm/$count;
        }
        else {$rowp = $rowm = "na";}
        #keep in a matrix
        $alternative_probp[$i] = $rowp;
        $alternative_probm[$i] = $rowm;
    }
    #print probabilities for the simulated model

    print OUTPUT "$filename\t";
    
    print OUTPUT "$parameter_altern{N_iterations}\t";
    print OUTPUT "$parameter_altern{Seed1}\t";
    print OUTPUT "$parameter_altern{N_loci}\t";
    print OUTPUT "$parameter_altern{Sfix_allthetas}\t";
    print OUTPUT "$parameter_altern{Range_thetant}\t";
    print OUTPUT "$parameter_altern{N_samples}\t";
    print OUTPUT "$parameter_altern{N_sites}\t";
    print OUTPUT "$parameter_altern{Recombination}\t";
    print OUTPUT "$parameter_altern{Mutations}\t";
    print OUTPUT "$parameter_altern{Thetaw}\t";
    print OUTPUT "$parameter_altern{Nintn}\t";
    print OUTPUT "$parameter_altern{Nrec}\t";
    print OUTPUT "$parameter_altern{Npast}\t";
    print OUTPUT "$parameter_altern{Tpast}\t";
    print OUTPUT "$parameter_altern{Npop}\t";
    print OUTPUT "$parameter_altern{Mig_rate}\t";
    print OUTPUT "$parameter_altern{Ifselection}\t";
    print OUTPUT "$parameter_altern{Pop_size}\t";
    print OUTPUT "$parameter_altern{Pop_sel}\t";
    print OUTPUT "$parameter_altern{Sel_nt}\t";
    print OUTPUT "$parameter_altern{Sinit}\t";
    
    print OUTPUT "Prob sim <= obs";

    for($i=0;$i<=$#headers_altern;$i++) {
        print OUTPUT "\t$alternative_probp[$i]";
    }
    print OUTPUT "\n";
    print OUTPUT "$filename\t";
    
    print OUTPUT "$parameter_altern{N_iterations}\t";
    print OUTPUT "$parameter_altern{Seed1}\t";
    print OUTPUT "$parameter_altern{N_loci}\t";
    print OUTPUT "$parameter_altern{Sfix_allthetas}\t";
    print OUTPUT "$parameter_altern{Range_thetant}\t";
    print OUTPUT "$parameter_altern{N_samples}\t";
    print OUTPUT "$parameter_altern{N_sites}\t";
    print OUTPUT "$parameter_altern{Recombination}\t";
    print OUTPUT "$parameter_altern{Mutations}\t";
    print OUTPUT "$parameter_altern{Thetaw}\t";
    print OUTPUT "$parameter_altern{Nintn}\t";
    print OUTPUT "$parameter_altern{Nrec}\t";
    print OUTPUT "$parameter_altern{Npast}\t";
    print OUTPUT "$parameter_altern{Tpast}\t";
    print OUTPUT "$parameter_altern{Npop}\t";
    print OUTPUT "$parameter_altern{Mig_rate}\t";
    print OUTPUT "$parameter_altern{Ifselection}\t";
    print OUTPUT "$parameter_altern{Pop_size}\t";
    print OUTPUT "$parameter_altern{Pop_sel}\t";
    print OUTPUT "$parameter_altern{Sel_nt}\t";
    print OUTPUT "$parameter_altern{Sinit}\t";

    print OUTPUT "Prob sim >= obs";

    for($i=0;$i<=$#headers_altern;$i++) {
        print OUTPUT "\t$alternative_probm[$i]";
    }
    print OUTPUT "\n";
}
close(OUTPUT);


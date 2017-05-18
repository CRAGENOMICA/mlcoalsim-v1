#!/usr/local/bin/perl

use Getopt::Long;
use File::Find;
use strict;

my (%options,$flag,$flag2,$flag3,@headers,$line,$percent,$i,$j,$k,$l,$count,@row,$value,$floor,$filenameout,$filename,@files);
my (@matrix_null,@inputnull_array,@matrix_null_sort,@null_intervals);
my (@matrix_altern,@headers_altern,@matrix_altern_array,@alternative_prob);
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
    'std=s'  => \$options{'std'}, #path to the file with the null model to compare
    'dir=s' => \$options{'dir'},    #path to the directory with alternative model files
    'perc0=s' => \$options{'perc0'},    #proportion of the distribution to the calculate power (1)
    'perc1=s' => \$options{'perc1'},    #proportion of the distribution to the calculate power (2)
);

#Usage
if($options{'std'} eq undef) {
   print "\nUsage: \nperl calc_power.pl -std (file with null model) -dir (directory with alternative files.out) -perc0 (proportion to calculate power) -perc1 (proportion to calculate power)\n\n";
   exit;
}

#start with the null model
open(INPUTNULL,"$options{'std'}") or die;

#reading the file
$flag=0;
while(<INPUTNULL>){
    if($_=~/(TD)/ || $_=~/\[0\]/) {$flag=1;}
    if($flag==2) {push(@inputnull_array,$_);}
    if($flag==1) {
        @headers = split(/\t/,$_);
        $flag=2;
    }
    if($flag==0) {
        #I also want to take the parameters values...
        if($_=~/N_iterations:\s+(\d+)/) {$parameter_null{N_iterations} = $1;}
        if($_=~/Seed1:\s+(\d+)/) {$parameter_null{Seed1} = $1;}

        if($_=~/N_loci:\s+(\d+)/) {$parameter_null{N_loci} = $1;}
        if($_=~/Sfix_allthetas:\s+(\d+)/) {$parameter_null{Sfix_allthetas} = $1;}
        if($_=~/Range_thetant:\s+(\d+)/) {$parameter_null{Range_thetant} = $1;}
        
        if($_=~/N_samples:\s+(.+)/) {$parameter_null{N_samples} = $1;}
        if($_=~/N_sites:\s+(.+)/) {$parameter_null{N_sites} = $1;}
        if($_=~/Recombination:\s+(.+)/) {$parameter_null{Recombination} = $1;}
        if($_=~/Mutations:\s+(.+)/) {$parameter_null{Mutations} = $1;}
        if($_=~/Thetaw:\s+(.+)/) {$parameter_null{Thetaw} = $1;}
        
        if($_=~/Nintn:\s+(\d+)/) {$parameter_null{Nintn} = $1;}
        if($_=~/Nrec:\s+(.+)/) {$parameter_null{Nrec} = $1;}
        if($_=~/Npast:\s+(.+)/) {$parameter_null{Npast} = $1;}
        if($_=~/Tpast:\s+(.+)/) {$parameter_null{Tpast} = $1;}
        
        if($_=~/Npop:\s+(\d+)/) {$parameter_null{Npop} = $1;}
        if($_=~/Mig_rate:\s+(\d+)/) {$parameter_null{Mig_rate} = $1;}
        
        if($_=~/Ifselection:\s+(.+)/) {$parameter_null{Ifselection} = $1;}
        if($_=~/Pop_size:\s+(\d+)/) {$parameter_null{Pop_size} = $1;}
        if($_=~/Pop_sel:\s+(.+)/) {$parameter_null{Pop_sel} = $1;}
        if($_=~/Sel_nt:\s+(.+)/) {$parameter_null{Sel_nt} = $1;}
        if($_=~/Sinit:\s+(.+)/) {$parameter_null{Sinit} = $1;}
    }
}
close(INPUTNULL);

#define columns in a new array
foreach $line (@inputnull_array) {
    my @vector_null = split(/\t/,$line);
    push(@matrix_null,[@vector_null]);
}

#set na, nan or "" to 1234567890
for($i=0;$i<=$#headers;$i++) {
    for($j=0;$j<=$#matrix_null;$j++) {
        if($matrix_null[$j][$i] eq "na" || $matrix_null[$j][$i] eq "nan" || $matrix_null[$j][$i] eq "") {
             $matrix_null[$j][$i] = 1234567890;
        }
    }
}

#sort each column and take the value at perc0%, perc1%, 1-perc0% and 1-perc1% 
for($i=0;$i<=$#headers;$i++) {
    #sort column 'i'
    @matrix_null_sort = sort {$a->[$i] <=> $b->[$i]} @matrix_null;
    $count = 0;
    for($j=0;$j<=$#matrix_null_sort;$j++) {
        if($matrix_null_sort[$j][$i] < 1234567890) {
             $count++;
        }
    }
    #calculate quantiles
    for($k=0;$k<4;$k++) {
        $row[$k] = "na";
        if($k==0) {$value = 1 + ($count-1) * $options{'perc0'};}
        if($k==1) {$value = 1 + ($count-1) * $options{'perc1'};}
        if($k==2) {$value = 1 + ($count-1) * (1.0 - $options{'perc1'});}
        if($k==3) {$value = 1 + ($count-1) * (1.0 - $options{'perc0'});}
        ($floor) = ($value=~/(\d+)/);
        if(($matrix_null_sort[$floor][$i] < 1234567890) || ($matrix_null_sort[$floor+1][$i]  < 1234567890)) {
            if($floor == $value) {
                if($k<2) {$row[$k] = $matrix_null_sort[$floor-1][$i];} #the matrix start from zero
                else {$row[$k] = $matrix_null_sort[$floor+0][$i];}
            }
            else {
                if($k<2) {$row[$k] = ($matrix_null_sort[$floor-1][$i] + $matrix_null_sort[$floor-2][$i])/2;}
                else {$row[$k] = ($matrix_null_sort[$floor+0][$i] + $matrix_null_sort[$floor+1][$i])/2};
            }
        }
    }
    #keep all quantilles in a matrix
    push(@null_intervals,[@row]);
}

#print this results in an output file

#open the output file
($options{'std'})=~/(\w+)\./;
$filenameout = ">"."$1"."_results".".txt";
open(OUTPUT,$filenameout) or die;

#print quantiles for std model
print OUTPUT "Name of the standard model input file: $options{'std'}\n\n";

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
print OUTPUT "Sinit\n";

print OUTPUT "$parameter_null{N_iterations}\t";
print OUTPUT "$parameter_null{Seed1}\t";
print OUTPUT "$parameter_null{N_loci}\t";
print OUTPUT "$parameter_null{Sfix_allthetas}\t";
print OUTPUT "$parameter_null{Range_thetant}\t";

if($parameter_null{N_loci} eq '1') {
    print OUTPUT "$parameter_null{N_samples}\t";
    print OUTPUT "$parameter_null{N_sites}\t";
    print OUTPUT "$parameter_null{Recombination}\t";
    print OUTPUT "$parameter_null{Mutations}\t";
    print OUTPUT "$parameter_null{Thetaw}\t";
    print OUTPUT "$parameter_null{Nintn}\t";
    print OUTPUT "$parameter_null{Nrec}\t";
    print OUTPUT "$parameter_null{Npast}\t";
    print OUTPUT "$parameter_null{Tpast}\t";
    print OUTPUT "$parameter_null{Npop}\t";
    print OUTPUT "$parameter_null{Mig_rate}\t";
    print OUTPUT "$parameter_null{Ifselection}\t";
    print OUTPUT "$parameter_null{Pop_size}\t";
    print OUTPUT "$parameter_null{Pop_sel}\t";
    print OUTPUT "$parameter_null{Sel_nt}\t";
    print OUTPUT "$parameter_null{Sinit}\n";
}
else {
    print OUTPUT "not shown\t";
    print OUTPUT "not shown\t";
    print OUTPUT "not shown\t";
    print OUTPUT "not shown\t";
    print OUTPUT "not shown\t";
    print OUTPUT "not shown\t";
    print OUTPUT "not shown\t";
    print OUTPUT "not shown\t";
    print OUTPUT "not shown\t";
    print OUTPUT "not shown\t";
    print OUTPUT "not shown\t";
    print OUTPUT "not shown\t";
    print OUTPUT "not shown\t";
    print OUTPUT "not shown\t";
    print OUTPUT "not shown\t";
    print OUTPUT "not shown\t";
}

print OUTPUT "\n";
print OUTPUT "Quantiles:";
for($i=0;$i<=$#headers;$i++) {
    print OUTPUT "\t$headers[$i]";
}
#print OUTPUT "\n";
for($k=0;$k<4;$k++) {
    if($k==0) {$percent = $options{'perc0'}*100; print OUTPUT 'Q['.$percent.'%]:'."\t";}
    if($k==1) {$percent = $options{'perc1'}*100; print OUTPUT 'Q['.$percent.'%]:'."\t";}
    if($k==2) {$percent = (1.0-$options{'perc1'})*100; print OUTPUT 'Q['.$percent.'%]:'."\t";}
    if($k==3) {$percent = (1.0-$options{'perc0'})*100; print OUTPUT 'Q['.$percent.'%]:'."\t";}
    for($i=0;$i<=$#headers;$i++) {
        print OUTPUT "$null_intervals[$i][$k]\t";
    }
    print OUTPUT "\n";
}
print OUTPUT "\n";

#then calculate the power for the alternative models

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
#take all files.out from the list directory (including subdirectories)
find( sub{ m/\.out$/i and unshift(@files,$File::Find::name);},$options{'dir'}, );
# @files is a vector with the name of each file.out in the columns

foreach $filename (@files) {
    #read each file
    open(INPUTALTERN,"$filename") or die;
    #introduce data in a matrix
    @matrix_altern = undef;
    unshift(@matrix_altern);
    @headers_altern = undef;
    unshift(@headers_altern);
    $flag=0;
    $flag2=0;
    while(<INPUTALTERN>){
        if($_=~/(TD)/ || $_=~/\[0\]/) {$flag=1;}
        if($flag==2) {
            if($flag2==0) {
                @matrix_altern[0] = $_;
                $flag2=1;
            }
            else {push(@matrix_altern,$_);}
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
        print "\nError: Not the same headers in null and alternative models.\n Next file...\n";
    }
    else {
        $l=0;
        for($i=0;$i<=$#headers;$i++) {
            if($headers[$i] ne $headers_altern[$i]) {
                print "\nError: Not the same headers in null and alternative models.\n Next file...\n";
                $l=1;
            }
        }
        if($l==0) {
            #define columns in a new array
            @matrix_altern_array = undef;
            unshift(@matrix_altern_array);
            $flag = 0;
            foreach $line (@matrix_altern) {
                my @vector_altern = split(/\t/,$line);
                if($flag == 0) {
                    @matrix_altern_array[0] = [@vector_altern];
                    $flag = 1;
                }
                else {push(@matrix_altern_array,[@vector_altern]);}
            }
            #for each column
            @alternative_prob = undef;
            unshift(@alternative_prob);
            $flag3 = 0;
            for($i=0;$i<=$#headers_altern;$i++) {
                #count the number of values and how many values are below and above the quantiles
                $count = 0;
                $row[0] = $row[1] = $row[2] = $row[3] = 0;
                for($j=0;$j<=$#matrix_altern_array;$j++) {
                    if($matrix_altern_array[$j][$i] ne "na" && $matrix_altern_array[$j][$i] ne "nan" && $matrix_altern_array[$j][$i] ne "") {
                        $count++;
                        if($matrix_altern_array[$j][$i] < $null_intervals[$i][0]) {$row[0] = $row[0] + 1;}
                        if($matrix_altern_array[$j][$i] < $null_intervals[$i][1]) {$row[1] = $row[1] + 1;}
                        if($matrix_altern_array[$j][$i] > $null_intervals[$i][2]) {$row[2] = $row[2] + 1;}
                        if($matrix_altern_array[$j][$i] > $null_intervals[$i][3]) {$row[3] = $row[3] + 1;}
                    }
                }
                #calculate proportion
                for($k=0;$k<4;$k++) {
                    if($count>0) {$row[$k] = $row[$k]/$count;}
                    else {$row[$k] = "na";}
                }
                #keep in a matrix
                if($flag3==0) {
                    @alternative_prob[0] = [@row];
                    $flag3 = 1;
                }
                else {push(@alternative_prob,[@row]);}
            }
            #print probabilities for alternative models

            for($k=0;$k<4;$k++) {
                print OUTPUT "$filename\t";
                
                print OUTPUT "$parameter_altern{N_iterations}\t";
                print OUTPUT "$parameter_altern{Seed1}\t";
                print OUTPUT "$parameter_altern{N_loci}\t";
                print OUTPUT "$parameter_altern{Sfix_allthetas}\t";
                print OUTPUT "$parameter_altern{Range_thetant}\t";
                
                if($parameter_null{N_loci} eq '1') { 
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
                }
                else {
                    print OUTPUT "not shown\t";
                    print OUTPUT "not shown\t";
                    print OUTPUT "not shown\t";
                    print OUTPUT "not shown\t";
                    print OUTPUT "not shown\t";
                    print OUTPUT "not shown\t";
                    print OUTPUT "not shown\t";
                    print OUTPUT "not shown\t";
                    print OUTPUT "not shown\t";
                    print OUTPUT "not shown\t";
                    print OUTPUT "not shown\t";
                    print OUTPUT "not shown\t";
                    print OUTPUT "not shown\t";
                    print OUTPUT "not shown\t";
                    print OUTPUT "not shown\t";
                    print OUTPUT "not shown\t";
                }
                if($k==0) {$percent = $options{'perc0'}*100; print OUTPUT "P(x\<".$percent."\%)";}
                if($k==1) {$percent = $options{'perc1'}*100; print OUTPUT "P(x\<".$percent."\%)";}
                if($k==2) {$percent = (1.0-$options{'perc1'})*100; print OUTPUT "P(x\>".$percent."\%)";}
                if($k==3) {$percent = (1.0-$options{'perc0'})*100; print OUTPUT "P(x\>".$percent."\%)";}
                for($i=0;$i<=$#headers;$i++) {
                    print OUTPUT "\t$alternative_prob[$i][$k]";
                }
                print OUTPUT "\n";
            }
        }
    }
}
close(OUTPUT);



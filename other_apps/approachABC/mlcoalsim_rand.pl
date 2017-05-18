#!/usr/local/bin/perl

use Getopt::Long;
use File::Basename; #parse file path, names and suffix
use strict;

my(%options,@initfile);
my($i,$j,$k,$flag,$it,$niter);
my(@par_name,@par_col,@par_min,@par_dist,@par_max,@par_num,$nr,@par_value,@linefields,@stdinarray);
my($fileline,$fileline2,$filenameout,$parline,$stdinline,$fileout);
my($filenameoutp,$path,$suffix);
my(@changep,@changev,$npop,$np,@vector,$p);

#options
GetOptions(
    'infile=s' => \$options{'infile'}, #path to the initial mlcoalsim input file
    'num_par=s' => \$options{'num_par'}, #number of parameters to modify
    'outfile=s' => \$options{'outfile'}, #name to the output results file
    'parfile=s' => \$options{'parfile'}, #name to the parameters output file
    'iter=s' => \$options{'iter'}, #number of iterations    
    'os=s' => \$options{'os'}, #PC, linux or mac
);

#Usage
if($options{'infile'} eq undef || $options{'num_par'} eq undef || $options{'outfile'} eq undef || $options{'parfile'} eq undef || $options{'os'} eq undef || $options{'iter'} eq undef) {
   print "\nUsage: \nperl mlcoalsim_grid.pl -infile (path to the initial mlcoalsim input file, with NO comments) -num_par (# of parameters to modify) -iter (# of iterations) -outfile (name of the output file) -parfile (name of the parameters output file -os (pc, linux, mac)\n\n";
   exit;
   }

#open the initial file
open(INPUTINIT,"$options{'infile'}") or die;
#reading the file and keep it in a matrix
while(<INPUTINIT>) {
    push(@initfile,$_);
}
close(INPUTINIT);

#Questions about the parameters to modify
#if(($stdinline = <>) !=~/\w/) {
    for($i=0;$i<$options{'num_par'};$i++) {
        printf "\nEnter the name of the parameter[%d]: ",$i+1;
        $par_name[$i] = <>;
        chomp($par_name[$i]);
        printf "Enter the column position of the parameter[%d] (the parameter is the column 0): ",$i+1;
        $par_col[$i] = <>;
        chomp($par_col[$i]);        
        printf "Enter the distribution for the parameter[%d] ((u)uniform, (l)og-uniform): ",$i+1;
        $par_dist[$i] = <>;
        chomp($par_dist[$i]);
        printf "Enter the minimum value of the parameter[%d]: ",$i+1;
        $par_min[$i] = <>;
        chomp($par_min[$i]);
        printf "Enter the maximum value of the parameter[%d]: ",$i+1;
        $par_max[$i] = <>;
        chomp($par_max[$i]);
    }

    print "\n\nPARAMETERS:\n";
    for($i=0;$i<$options{'num_par'};$i++) {
        printf "Name of the parameter[%d]: %s\n",$i+1,$par_name[$i];
        printf "Column position of the parameter[%d] (the parameter is the column 0): %d\n",$i+1,$par_col[$i];
        printf "Distribution for the parameter[%d] ((u)niform, (l)og-uniform): %s\n",$i+1,$par_dist[$i];
        printf "Minimum value of the parameter[%d]: %g\n",$i+1,$par_min[$i];
        printf "Maximum value of the parameter[%d]: %g\n",$i+1,$par_max[$i];
    }

#change thetant_max, theta, recnt_max, recombination, thetant_min, recnt_min in case npop is modified
$npop = undef;
for($i=0;$i<$options{'num_par'};$i++) {
    if($par_name[$i] =~/npop/) {
        $changep[0] = "thetant_min";
        $changep[1] = "thetant_max";
        $changep[2] = "thetaw";
        $changep[3] = "recnt_min";
        $changep[4] = "recnt_max";
        $changep[5] = "recombination";
        #$changep[6] = "dist_out";
        $np = $i;
        foreach $fileline (@initfile) {
            if($fileline =~/"/) {next;}
            if($fileline =~/npop\s/) {
                @linefields = split(/\s+/,$fileline);
                $npop = $linefields[1];
            }
        }
    }
}
if($npop) {
    for($i=0;$i<6;$i++) {
        foreach $fileline (@initfile) {
            if($fileline =~/"/) {next;}
            if($fileline =~/$changep[$i]/) {
                @vector = split(/\s+/,$fileline);
                for($j=0;$j<=$#vector;$j++) {
                    $changev[$i][$j] = $vector[$j];
                }
            }
        }
    }
}

#open the parameters output file
open(OUTPAR,">$options{'parfile'}") or die;
for($i=0;$i<$options{'num_par'};$i++) {print OUTPAR "$par_name[$i]\t";}
print OUTPAR "\n";

#Generate input files, run files
$it = 1;
$niter = $options{'iter'};
print "\nGenerating $niter runs with mlcoalsim...\n";
#random numbers
srand(123456);

while($it <= $niter) {
    #look for next parameters
    $flag = 0;
    for($i=$options{'num_par'}-1;$i>=0;$i--) {
        $nr = rand;
        if($par_dist[$i] =~/u/) {
            $par_value[$i] = $par_min[$i] + ($par_max[$i] - $par_min[$i]) * $nr;
        }
        if($par_dist[$i] =~/l/) {
            $par_value[$i] = exp(log($par_min[$i]) + (log($par_max[$i]) - log($par_min[$i])) * $nr);
        }
        if($par_name[$i] =~/npop/) {$par_value[$i] = int($par_value[$i]);}
    }
    #create input mlcoalsim file
    $filenameout = "inputfile_"."$it"."_"."$niter".".txt";
    open(OUTPUT,">$filenameout") or die;
    #include all parameters
    foreach $fileline (@initfile) {
        chomp($fileline);
        if($fileline =~/"/) {next;}
        $fileline2 = $fileline;
        if($fileline =~/seed1\s/) {
            @linefields = split(/\s+/,$fileline);
            $linefields[1] = int(rand 1e8);
            $fileline2 = join "\t", @linefields;
        }
        else {
            for($i=0;$i<$options{'num_par'};$i++) {
                if($fileline =~/$par_name[$i]\s/) {
                    @linefields = split(/\s+/,$fileline);
                    $linefields[$par_col[$i]] = $par_value[$i];
                    $fileline2 = join "\t", @linefields;
                 }
            }
            if($npop) {
                for($k=0;$k<6;$k++) {
                    if($fileline =~/$changep[$k]\s/) {
                        @linefields = split(/\s+/,$fileline);
                        for($j=1;$j<=$#linefields;$j++) {
                            $p = undef;
                            for($i=0;$i<$options{'num_par'};$i++) {
                                if($changep[$k] eq $par_name[$i] && $j == $par_col[$i]) {
                                    $p = $par_value[$i];
                                }
                            }
                            if($p == undef) {$p = $changev[$k][$j];}
                            #if($i < 6) {
                                $linefields[$j] = $npop / $par_value[$np] * $p;
                            #}
                            #else {$linefields[$j] = $par_value[$np] / $npop * $p;}
                        }
                        $fileline2 = join "\t", @linefields;
                    }
                }
            }
        }
        print OUTPUT $fileline2."\n";
    }
    close(OUTPUT);
    
    #run and erase input file
    $fileout = $options{'outfile'};
    if($options{'os'} eq "mac") {system("./mlcoalsim_osx $filenameout $fileout > ./res.txt");}
    if($options{'os'} eq "linux") {system("./mlcoalsim $filenameout $fileout > ./res.txt");}
    if($options{'os'} eq "pc") {system("mlcoalsim $filenameout $fileout > ./res.txt");}
    unlink "$filenameout";
    
   #printing parameters
   open(READOUT,"./res.txt") or die;
   $flag=0;
   while(<READOUT>) {
      if($_ =~/succesful/) {$flag=1;}
   }   
   close(READOUT);
   if($flag==1) {
     for($i=0;$i<$options{'num_par'};$i++) {print OUTPAR "$par_value[$i]\t";}
     print OUTPAR "\n";
     $it += 1;
  }
   else {
      print"\nRUN #$it FAILED: \n";
      for($i=0;$i<$options{'num_par'};$i++) {print "$par_name[$i]\t";}
      print "\n";
      for($i=0;$i<$options{'num_par'};$i++) {print "$par_value[$i]\t";} 
      print "\n";
      open(READOUT,"./res.txt") or die;
      while(<READOUT>) {
         print"$_";
      }     
      close(READOUT);
      print"\n";
   }
   unlink"./res.txt";  
}

close(OUTPAR);
($filenameoutp,$path,$suffix) = fileparse($options{'outfile'},qr/\.[^.]*/);    
$fileout = "$path"."$filenameoutp"."_1".".out";
unlink $fileout;
print "\nmlcoalsim_grid.pl process finished.\n";
exit();


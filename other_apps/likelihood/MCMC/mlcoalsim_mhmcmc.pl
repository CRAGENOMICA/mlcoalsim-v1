#!/usr/local/bin/perl

use Getopt::Long;
use File::Basename; #parse file path, names and suffix
use strict;

my(%options,@initfile);
my($i,$j,$k,$p,$it,$niter);
my(@par_name,@par_col,@par_dist,@par_delta,@par_value,@par_min,@par_max);
my(@linefields,@stdinarray);
my($fileline,$fileline2,$filenameout,$parline,$stdinline,$fileout);
my(@likelihoods,@likelihoodsc,$likelihoodline,$likelihood,$likelihoodc);
my(@candidate,$headerline,$filenameoutmc);
my($path,$filenameoutp,$suffix);
my(@changep,@changev,$npop,$np,@vector,$burnin,$ittot,$rn,$flag);

#options
GetOptions(
    'infile=s' => \$options{'infile'}, #path to the initial mlcoalsim input file
    'numpar=s' => \$options{'num_par'}, #number of parameters to modify
    'burniter=s' => \$options{'burn_iter'}, #number of burn-in iterations in the chain
    'mcmciter=s' => \$options{'mcmc_iter'}, #number of iterations in the chain
    'mlcout=s' => \$options{'mlcoalsimout'}, #mlcoalsim output file name
    'mcmcout=s' => \$options{'mcmc_out'}, #MHMCMC output file name
    'parfile=s' => \$options{'parfile'}, #name to the parameters output file
    'os=s' => \$options{'os'}, #PC, linux or mac
);

#Usage
if($options{'infile'} eq undef || $options{'num_par'} eq undef || $options{'burn_iter'} eq undef || $options{'mcmc_iter'} eq undef || $options{'mlcoalsimout'} eq undef || $options{'mcmc_out'} eq undef  || $options{'parfile'} eq undef || $options{'os'} eq undef) {
   print "\nUsage: \nperl mlcoalsim_mhmcmc.pl -infile (path to the initial mlcoalsim input file, with NO COMMENTS) -numpar (number of parameters to modify) -burniter (burn-in iterations) -mcmciter (iterations in the chain) -mlcoalsimout (mlcoalsim output file name) -mcmcout (MHMCMC output file name) -parfile (name of the parameters output file -os (pc, linux, mac)\n\n";
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
        printf "Enter the steps for the parameter[%d] (aritmethic[a],logarithmic[l]): ",$i+1;
        $par_dist[$i] = <>;
        chomp($par_dist[$i]);
        printf "Enter the increasing delta value for MHMCMC of the parameter[%d]: ",$i+1;
        $par_delta[$i] = <>;
        chomp($par_delta[$i]);
        printf "Enter the minimum value of the parameter[%d]: ",$i+1;
        $par_min[$i] = <>;
        chomp($par_min[$i]);
        printf "Enter the maximum value of the parameter[%d]: ",$i+1;
        $par_max[$i] = <>;
        chomp($par_max[$i]);
        printf "Enter the initial value of the parameter[%d]: ",$i+1;
        $par_value[$i] = <>;
        chomp($par_value[$i]);
    }
    print "\nPARAMETERS:\n";
    for($i=0;$i<$options{'num_par'};$i++) {
        printf "\nName of the parameter[%d]: %s\n",$i+1,$par_name[$i];
        printf "Column position of the parameter[%d] (the parameter is the column 0): %d\n",$i+1,$par_col[$i];
        printf "Steps for the parameter[%d] (aritmethic[a],logarithmic[l]): %s\n",$i+1,$par_dist[$i];
        printf "Increasing delta value for MHMCMC of the parameter[%d]: %g\n",$i+1,$par_delta[$i];
        printf "Minimum value of the parameter[%d]: %g\n",$i+1,$par_min[$i];
        printf "Maximum value of the parameter[%d]: %g\n",$i+1,$par_max[$i];
        printf "Initial value of the parameter[%d]: %g\n",$i+1,$par_value[$i];
    }
#}
#else {
#    @stdinarray = split(/\s+/,$stdinline);
#    print "\nPARAMETERS:\n";
#    for($i=0,$j=5;$i<$options{'num_par'};$i++) {
#        printf "Name of the parameter[%d]: %s\n",$i+1,$stdinarray[$j++];
#        $par_name[$i] = $stdinarray[$j];
#        printf "Column position of the parameter[%d] (the parameter is the column 0): %d\n",$i+1,$stdinarray[$j++];
#        $par_col[$i] = $stdinarray[$j];
#        printf "Steps for the parameter[%d] (aritmethic[a],logarithmic[l]): %s",$i+1,$stdinarray[$j++];
#        $par_dist[$i] = $stdinarray[$j];
#        printf "Increasing delta value for MHMCMC of the parameter[%d]: %g\n",$i+1,$stdinarray[$j++];
#        $par_delta[$i] = $stdinarray[$j];
#        printf "Initial value of the parameter[%d]: %g\n",$i+1,$stdinarray[$j++];
#        $par_value[$i] = $stdinarray[$j];
#    }
#}
#change thetant_max, theta, dist_out, recnt_max, recombination, thetant_min, recnt_min in case npop is modified
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
$it = 0;
$niter = $options{'mcmc_iter'};
$burnin = $options{'burn_iter'};
$ittot = $niter + $burnin;
print "\nGenerating $niter+$burnin runs with mlcoalsim...\n";
#random numbers
srand(123456);

#FIRST ITERATION
while($it==0) {
    #create file for keeping the distribution of MHMCMC
    $filenameoutmc = $options{'mcmc_out'};
    open(OUTPUTLK,">$filenameoutmc") or die;
    #create input mlcoalsim file
    $filenameout = "inputfile_"."$it"."_"."$ittot".".txt";
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
    ($filenameoutp,$path,$suffix) = fileparse($options{'mlcoalsimout'},qr/\.[^.]*/);    
    $fileout = "$path"."$filenameoutp"."$suffix";
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
      #if($it >= $burnin) {
      #  for($i=0;$i<$options{'num_par'};$i++) {print OUTPAR "$par_value[$i]\t";}
      #  print OUTPAR "\n";
      #}
    
     #Read header and print initial values
     open(OUTIN,"$fileout") or die;
     $headerline = <OUTIN>;
     chomp($headerline);
     print OUTPUTLK "$headerline\n";
     close(OUTIN);
     #Read the last line of the file likelihoods.txt
     open(OUTLV,"$path"."$filenameoutp"."_1.out") or die;
     #reading the file and keep it in a matrix
     $likelihoodline = <OUTLV>;
     chomp($likelihoodline);
     @likelihoods = split(/\s+/,$likelihoodline);
     $likelihood = $likelihoods[0];
     #for($i=0;$i<=$#likelihoods;$i++) {print OUTPUTLK "$likelihoods[$i]\t";}
     #print OUTPUTLK "\n";
     close(OUTLV);
     unlink("$path"."$filenameoutp"."_1.out");    
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

#LOOP
while($it < $niter + $burnin) {
    #look for next parameters random walk metropolis
    for($i=0;$i<$options{'num_par'};$i++) {
       do {
            $rn = rand;
            #arithmetic
            if($par_dist[$i] =~/a/) { 
                $candidate[$i] = ($par_value[$i] - $par_delta[$i]) + ($rn * 2 * $par_delta[$i]);
            }
            #logarithmic
            if($par_dist[$i] =~/l/) { 
                $candidate[$i] = exp((log($par_value[$i]) - log($par_delta[$i])) + ($rn * 2 * log($par_delta[$i])));
            }
        } while($candidate[$i] < $par_min[$i] || $candidate[$i] > $par_max[$i]);
        if($par_name[$i] =~/npop/) {$candidate[$i] = int($candidate[$i]);
        }
    }
    #probabilities
    #create input mlcoalsim file
    $filenameout = "inputfile_"."$it"."_"."$ittot".".txt";
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
                    $linefields[$par_col[$i]] = $candidate[$i];
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
                                    $p = $candidate[$i];
                                }
                            }
                            if($p == undef) {$p = $changev[$k][$j];}
                            #if($i < 6) {
                                $linefields[$j] = $npop / $candidate[$np] * $p;
                            #}
                            #else {$linefields[$j] = $candidate[$np] / $npop * $p;}
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
    ($filenameoutp,$path,$suffix) = fileparse($options{'mlcoalsimout'},qr/\.[^.]*/);    
    $fileout = "$path"."$filenameoutp"."$suffix";
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
     if($it >= $burnin) {
        for($i=0;$i<$options{'num_par'};$i++) {print OUTPAR "$par_value[$i]\t";}
        print OUTPAR "\n";
     }
    
     #Read the last line of the file likelihoods.txt
     open(OUTLV,"$path"."$filenameoutp"."_1.out") or die;
     $likelihoodline = <OUTLV>;
     chomp($likelihoodline);
     @likelihoodsc = split(/\s+/,$likelihoodline);
     $likelihoodc = $likelihoodsc[0];
     unlink("$path"."$filenameoutp"."_1.out");    
     close(OUTLV);
    
     #accept/reject step    
     if(2*log(rand) <= $likelihoodc - $likelihood) {
        for($i=0;$i<$options{'num_par'};$i++) {
            $par_value[$i] = $candidate[$i];
        }
        @likelihoods = @likelihoodsc;
        $likelihood = $likelihoodc;
     }
     if($it >= $burnin) {
        for($i=0;$i<=$#likelihoods;$i++) {print OUTPUTLK "$likelihoods[$i]\t";}
        print OUTPUTLK "\n";
     }
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
close(OUTPUTLK);
print "mlcoalsim_mhmcmc.pl process finished.\n";
exit();


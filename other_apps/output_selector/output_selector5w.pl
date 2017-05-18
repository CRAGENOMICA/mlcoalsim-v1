###output_selector.pl for mlcoalsim09v3##
use strict;
use Getopt::Long;
use Statistics::Descriptive;
use Statistics::Distributions;

#Opciones
my %options;
GetOptions(
'cols=s'=>\$options{'cols'},
'in=s'=>\$options{'in'},
'obs=s'=>\$options{'obs'},
'out=s'=>\$options{'out'},
'emp=s'=>\$options{'emp'},
'tail=s'=>\$options{'tail'},
'wt=s'=>\$options{'wt'}
);

unless($options{'cols'} or $options{'in'} or $options{'obs'} or $options{'out'} or $options{'tail'}){
print "\nUsage:\n";
print "perl output-selector5.pl -in inputfile(mlcoalsim output) -obs (obs_freq_file) -cols selected_cols_from_inputfile(ej 1:5:9) -out outputfile -emp (calculate empirical P: y/n) -tail (left/right/both) -wt (file with weight values for each iteration)\n\n";
exit;
}

#Abre fichero para la salida
open(OUTPUT,">$options{'out'}") or die "Error happens trying to open $options{'out'}\n";

#Declaraciones varias
my ($b2,$t,$size,@size2);
my (@fila,@filaw,@columna,@stat,@cols,@array_percentile_prev,@sorteds_array,@Sum);

if ($options{'cols'}){
  @cols=split(/:/,$options{'cols'});
}
if ($options{'emp'} eq undef) {$options{'emp'} = "y";}
elsif($options{'emp'} =~/y/) {$options{'emp'} = "y";}
    elsif($options{'emp'} =~/n/) {$options{'emp'} = "n";}

if ($options{'tail'} eq undef) {$options{'tail'} = "both";}
if ($options{'tail'} ne "both" and $options{'tail'} ne "left" and $options{'tail'} ne "right") {$options{'tail'} = "both";}

print "\nversion:\toutput_selector5\nOptions:\ncols\t$options{'cols'}\nin\t$options{'in'}\nobs\t$options{'obs'}\nout\t$options{'out'}\nweight\t$options{'wt'}\ntail\t$options{'tail'}\nemp\t$options{'emp'}\n\n";
print "Reading input files...";
print OUTPUT "\nversion:\toutput_selector5\nOptions:\ncols\t$options{'cols'}\nin\t$options{'in'}\nobs\t$options{'obs'}\nout\t$options{'out'}\nweight\t$options{'wt'}\ntail\t$options{'tail'}\nemp\t$options{'emp'}\n\n";
print OUTPUT "COMPUTING THE COMPOSITE PROBABILITY:\n\n";
#Lee los datos del fichero con los datos de lo observado
my (@obs,@obs_pre,@erase_col);
open(OBS,"$options{'obs'}") or die "Can not open observed file\n";
while(<OBS>){
    chomp($_);
    if($_ !~/(\(|\[|\()/ || $_=~/\bna\b/) {
        if ($_=~/\d+/) {
            if (length($_)-rindex($_,"\t")==1){chop($_);}
            @obs_pre=split(/\t/,$_);
            my $cnt_cols=0;
            foreach my $element (@cols){
                if(($obs_pre[$element-1] eq "na") || ($obs_pre[$element-1] == -10000)){
                    push(@erase_col,$cnt_cols);
                    $cnt_cols++;
                    print "\nObserved column $element does not exist, eliminated from analysis...";
                    print OUTPUT "\nObserved column $element does not exist, eliminated from analysis...";
                    next;
                }
                if ($element > scalar(@obs_pre)) {
                    print "\nObserved column ($element) doesn't exist, eliminated from analysis...";
                    print OUTPUT "\nObserved column ($element) doesn't exist, eliminated from analysis...";
                    push(@erase_col,$cnt_cols);
                    $cnt_cols++;
                    next;
                }
                push (@obs,$obs_pre[$element-1]);
                $cnt_cols++;
            }
        }
    }
}
for(my $j=$#erase_col;$j>=0;$j--){
    my $element = $erase_col[$j];
    splice(@cols,$element,1);
}
close (OBS);

#Abre fichero del output del mlcoalsim y lo lee eliminando la cabecera
open (INPUT,"$options{'in'}") or die "Can' t open simulated file\n";
$b2="F";
my $s=0;
while(<INPUT>){
    if ($_=~/value/ or $_=~/\[/ or $_=~/avg/ or $_=~/\(/ ){
        $b2="T";
    }
    if ($b2 eq "T"){
        chomp($_);
        @columna=split(/\t/,$_);
        my @columna2;
        my $j=0;
        foreach my $col (@cols){
            push(@columna2,$columna[$col-1]);
            $stat[$j] = Statistics::Descriptive::Full->new() if ($s==0);
            if($columna[$col-1] ne "na") {
                $stat[$j]->add_data($columna[$col-1])  if ($s!=0) ;
            }
            $j++;
        }
        push(@fila,[@columna2]);
        @columna2= ();
        $s++;
    }
}
close(INPUT);

#Abre fichero de pesos (weights) de ABC y lo lee eliminando la cabecera
if($options{'wt'}) {
  open (WEIGHT,"$options{'wt'}") or die "Can' t open weight file\n";
  my $s2=0;
  my $b2="F";
  while(<WEIGHT>){
    if ($_=~/Weight/ ){
        $b2="T";
    }
    if ($b2 eq "T"){
      chomp($_);
      push(@filaw,$_);
      $s2++;
    }
  }
  close(WEIGHT);
  if($s2 != $s) {
    print "The weight file and the input file have not the same number of rows. Exiting..";
    exit;
  }
}
else {
  #put he same weight to all values (1)
  for(my $xc=0; $xc<$s; $xc++) {
    push(@filaw,1);
  }
}


my $j=0;
foreach my $element (@cols){
    if($stat[$j]->count <= 0) {
        print "\nSimulated column $element can not be calculated, not included in analysis...";
    }
    $j++;
}
print "\nCalculating composite probabilities from empirical simulations...";

# Calcula las probabilidades de los datos y las probabilidades de lo observado en los datos
my @pro_obs;
push(@pro_obs,"P1");
my @percentile;
my @perc_allstats;
$j=1;
foreach my $element (@cols){
    if(($stat[$j-1]->count) > 0) {
        my (@array_percentile,@prob_col);
        my $array_per;
        #include all data (for each column) in a single array named @array_percentile
        $size = 0; #total size ($#fila)
        $size2[$j-1]= 0; #size2 is the total sum of the weight
        $percentile[0] = 0;
        for my $p (1 .. $#fila){
            my @temp;
            $temp[0]=$p;
            $temp[1]=$fila[$p][$j-1];
            $temp[2]=$filaw[$p];
            push (@array_percentile,[@temp]);
            $size += $filaw[$p];
            if($temp[1] ne "na") {
              $size2[$j-1] += $filaw[$p];
              $percentile[0] += $temp[2] * $temp[1]; #median calculated
            }
       }
        
        my @sorted_array_per=sort{@$a[1]<=>@$b[1]} @array_percentile; #keep in 0 the position and in 1 the value (sorted by the value). In 2 is the weight
        #array_percentile is the column desired. sorted_array_per is this column but sorted    
 
        if($options{'emp'} eq "y" || $options{'emp'} eq "Y") {              
            my %util_hash;#the position (weight) of each different and sorted element. This is a hash indicating the number of the cummulated weight where each value is (for repeated only the first).
            my $element2_ant = -100000;
            my $pos=0; #pos is the cummulative probability (weight) of each element in inverse order
            for(my $p=$#fila-1;$p>=0;$p--) {
                my $element2 = $sorted_array_per[$p][1];
                if($element2 ne "na") {
                    if($element2_ant!=$element2) {
                        $util_hash{$element2} = $pos;
                    }
                    $element2_ant=$element2;
                    $pos += $sorted_array_per[$p][2];                    
                }
            }
            
            $element2_ant = -100000;
            my $cnt_hash_array=0;
            my $cnt_na=0;
            my $cnt=0;my $cnt2=0;
            my $cnt_hash_array_post=0;
            my $count_na=0;
            my $count_value=0;
            my $sumper=0; my $psumper; #sum of all probabilities to have the percentiles
        
            for(my $p=0;$p<$#fila;$p++) {
                my $element2;
                my $add_equal;
                $element2 = $sorted_array_per[$p][1];
                if($element2 ne "na") {
                    my $probab1;my $probab2=1;my $pro_obs;
                    if($element2!=$element2_ant) {#it's faster when there are many repeated elements
                        $probab1=($cnt_hash_array)/($size2[$j-1]-$sorted_array_per[$p][2]);
                        if($probab1 < 1e-6) {$probab1 = 0;}
                        $probab2 = ($util_hash{$element2})/($size2[$j-1]-$sorted_array_per[$p][2]);
                        if($probab2 < 1e-6) {$probab2 = 0;}
                        $add_equal = (($size2[$j-1]-$sorted_array_per[$p][2]) - ($cnt_hash_array+$util_hash{$element2}))/($size2[$j-1]-$sorted_array_per[$p][2]);
                        if($add_equal < 1e-6) {$add_equal = 0;}
                        #$pro_obs=(0.5-abs(0.5-($probab1)))+(0.5-abs(0.5-($probab2))) + $add_equal;
                        $pro_obs=0;
                        if($options{'tail'} eq "both") {
                            $pro_obs += (0.5-abs(0.5-$probab1)) + (0.5-abs(0.5-$probab2));
                        }
                        if($options{'tail'} eq "left") {
                            $pro_obs += $probab1;
                        }
                        if($options{'tail'} eq "right") {
                            $pro_obs += $probab2;
                        }
                        $pro_obs += $add_equal;
                        if($pro_obs < 1e-6) {$pro_obs = 0;}
                        if($pro_obs > 1-1e-6) {$pro_obs = 1;}
                        #debugging
                        if($pro_obs < 0 || $pro_obs > 1) {
                            print "\nError in script: ";
                            print "p1:$probab1\tp2:$probab2\tp:$p\tj:$j\tn:$util_hash{$element2_ant}";
                            print "\telement2:$element2\telement2_ant:$element2_ant\tsize2:$size2[$j-1]\tadd:$add_equal\n";
                            exit;
                        }
                        
                        #redo matrix with the same order than original
                        $cnt_hash_array_post = $count_value + $count_na + 1;#$cnt_hash_array + $cnt_na + 1;
                        while($cnt_hash_array_post < $size &&
                              ($sorted_array_per[$p][1] == $sorted_array_per[$cnt_hash_array_post][1] ||
                              $sorted_array_per[$cnt_hash_array_post][1] eq "na")) {
                                $cnt_hash_array_post++;
                        }
                        for (my $m=$count_value + $count_na;$m<$cnt_hash_array_post;$m++){
                            $sorteds_array[$sorted_array_per[$m][0]][$j-1]=$pro_obs;
                        }
                    }
                    if($element2 < $obs[$j-1]) {$cnt += $sorted_array_per[$p][2];}
                    elsif($element2 > $obs[$j-1]) {$cnt2 += $sorted_array_per[$p][2];}
                    $cnt_hash_array += $sorted_array_per[$p][2];
                    $count_value++;
                    $element2_ant=$element2;
 
                    #percentile calculations
                    $psumper = $sorted_array_per[$p][2]/$size2[$j-1];
                    $sumper += $psumper; #percentile
                    
                    if($sumper >= 0.01 && ($sumper - $psumper) < 0.01 ) {
                      $percentile[1] = $sorted_array_per[$p][1];
                    }
                    if($sumper >= 0.025 && ($sumper - $psumper) < 0.025 ) {
                      $percentile[2] = $sorted_array_per[$p][1];
                    }
                    if($sumper >= 0.05 && ($sumper - $psumper) < 0.05 ) {
                      $percentile[3] = $sorted_array_per[$p][1];
                    }
                    if($sumper >= 0.25 && ($sumper - $psumper) < 0.25 ) {
                      $percentile[4] = $sorted_array_per[$p][1];
                    }
                    if($sumper >= 0.50 && ($sumper - $psumper) < 0.50 ) {
                      $percentile[5] = $sorted_array_per[$p][1];
                    }
                    if($sumper >= 0.75 && ($sumper - $psumper) < 0.75 ) {
                      $percentile[6] = $sorted_array_per[$p][1];
                    }
                    if($sumper >= 0.95 && ($sumper - $psumper) < 0.95 ) {
                      $percentile[7] = $sorted_array_per[$p][1];
                    }
                    if($sumper >= 0.975 && ($sumper - $psumper) < 0.975 ) {
                      $percentile[8] = $sorted_array_per[$p][1];
                    }
                    if($sumper >= 0.99 && ($sumper - $psumper) < 0.99 ) {
                      $percentile[9] = $sorted_array_per[$p][1];
                    }
               }
                else {
                    $sorteds_array[$sorted_array_per[$count_value+$count_na][0]][$j-1]="na";
                    $cnt_na += $sorted_array_per[$p][2];
                    $count_na++;
                }
            }
            my $probab1=($cnt/$size2[$j-1]);
            my $probab2=($cnt2/$size2[$j-1]);
            my $add_equal = ($size2[$j-1] -($cnt+$cnt2))/$size2[$j-1];

            my $pro_obs=0;
            if($options{'tail'} eq "both") {
                $pro_obs += (0.5-abs(0.5-$probab1)) + (0.5-abs(0.5-$probab2));
            }
            if($options{'tail'} eq "left") {
                $pro_obs += $probab1;
            }
            if($options{'tail'} eq "right") {
                $pro_obs += $probab2;
            }
            $pro_obs += $add_equal;

            if($pro_obs <= 0) {$pro_obs = 1.0/$size2[$j-1];}
            push (@pro_obs, $pro_obs);
        }
        else {
            my $cnt=0;my $cnt2=0;
            my $sumper=0; my $psumper; #sum of all probabilities to have the percentiles
            for(my $p=0;$p<$#fila;$p++) {
                my $element2;
                $element2 = $array_percentile[$p][1];
                if($element2 ne "na") {
                    if($element2 < $obs[$j-1]) {$cnt += $sorted_array_per[$p][2];}
                    elsif($element2 > $obs[$j-1]) {$cnt2 += $sorted_array_per[$p][2];}

                    #percentiles calculations
                    $psumper = $sorted_array_per[$p][2]/$size2[$j-1];
                    $sumper += $psumper; #percentile
                    
                    if($sumper >= 0.01 && ($sumper - $psumper) < 0.01 ) {
                      $percentile[1] = $sorted_array_per[$p][1];
                    }
                    if($sumper >= 0.025 && ($sumper - $psumper) < 0.025 ) {
                      $percentile[2] = $sorted_array_per[$p][1];
                    }
                    if($sumper >= 0.05 && ($sumper - $psumper) < 0.05 ) {
                      $percentile[3] = $sorted_array_per[$p][1];
                    }
                    if($sumper >= 0.25 && ($sumper - $psumper) < 0.25 ) {
                      $percentile[4] = $sorted_array_per[$p][1];
                    }
                    if($sumper >= 0.50 && ($sumper - $psumper) < 0.50 ) {
                      $percentile[5] = $sorted_array_per[$p][1];
                    }
                    if($sumper >= 0.75 && ($sumper - $psumper) < 0.75 ) {
                      $percentile[6] = $sorted_array_per[$p][1];
                    }
                    if($sumper >= 0.95 && ($sumper - $psumper) < 0.95 ) {
                      $percentile[7] = $sorted_array_per[$p][1];
                    }
                    if($sumper >= 0.975 && ($sumper - $psumper) < 0.975 ) {
                      $percentile[8] = $sorted_array_per[$p][1];
                    }
                    if($sumper >= 0.99 && ($sumper - $psumper) < 0.99 ) {
                      $percentile[9] = $sorted_array_per[$p][1];
                    }
                }
            }
            my $probab1=($cnt/$size2[$j-1]);
            my $probab2=($cnt2/$size2[$j-1]);
            my $add_equal = ($size2[$j-1] -($cnt+$cnt2))/$size2[$j-1];
            
            my $pro_obs=0;
            if($options{'tail'} eq "both") {
                $pro_obs += (0.5-abs(0.5-$probab1)) + (0.5-abs(0.5-$probab2));
            }
            if($options{'tail'} eq "left") {
                $pro_obs += $probab1;
            }
            if($options{'tail'} eq "right") {
                $pro_obs += $probab2;
            }
            $pro_obs += $add_equal;
            
            if($pro_obs <= 0) {$pro_obs = 1.0/$size2[$j-1];}
            push (@pro_obs, $pro_obs);
        }
        push(@perc_allstats,[@percentile]);
    }
    $j++;
}


#Impresi—n final del output
print OUTPUT "\n\nStatistic\tmean\t1%\t2.5%\t5%\t25%\t50%\t75%\t95%\t97.5%\t99%\tvalid_iter\n";
$j=0;
foreach my $element (@cols){
    if($stat[$j]->count > 0) {
      printf OUTPUT "%s",$fila[0][$j];
            
      printf OUTPUT "\t%.6f",$perc_allstats[$j][0]/$size2[$j];
      printf OUTPUT "\t%.6f",$perc_allstats[$j][1];
      printf OUTPUT "\t%.6f",$perc_allstats[$j][2];
      printf OUTPUT "\t%.6f",$perc_allstats[$j][3];
      printf OUTPUT "\t%.6f",$perc_allstats[$j][4];
      printf OUTPUT "\t%.6f",$perc_allstats[$j][5];
      printf OUTPUT "\t%.6f",$perc_allstats[$j][6];
      printf OUTPUT "\t%.6f",$perc_allstats[$j][7];
      printf OUTPUT "\t%.6f",$perc_allstats[$j][8];
      printf OUTPUT "\t%.6f",$perc_allstats[$j][9];
      
      printf OUTPUT "\t%ld" ,$stat[$j]->count;
      print  OUTPUT "\n";
      
    }
    else {printf OUTPUT "Simulated column %d can not be calculated\n",$element;}
    $j++;
}

#print "\n\nStatistic\tOBSERVED\tP1\n";
print OUTPUT "\n\nStatistic\tOBSERVED\tP1\n";

$j=0;
foreach my $element (@cols){
    if($stat[$j]->count > 0) {
        printf OUTPUT "%s",$fila[0][$j];
        printf OUTPUT "\t%.6f",$obs[$j];
        printf OUTPUT "\t%.6f",$pro_obs[$j+1];
        print  OUTPUT "\n";
    }
    $j++;
}

print OUTPUT "\n";

#print OUTPUT2 (all probabilities for Jacknife) AND CALCULATE THE SUM OF THE 2LOG P (I think the values have not to be weighted (inside) until counting the Sum)
@Sum = ();
$j=0;
my $i;
if($options{'emp'} eq "y" || $options{'emp'} eq "Y") {
    my $output2 = "$options{'out'}"."_distPC.out";
    open(OUTPUT2,">$output2") or die "Error happens trying to create the file $output2\n";
    
    $j=0;
    foreach my $element (@cols) {
        if($stat[$j]->count > 0) {
            printf OUTPUT2 "%s\t",$fila[0][$j];
        }
        $j++;
    }
    print OUTPUT2 "\tP1\n";
    
    for my $i (1 .. $#sorteds_array){
        $Sum[$i]=0;
        my $value_ant=0;
        my $cna = 0;
        for my $j (0 .. $#{$sorteds_array[$i]}){
            if($sorteds_array[$i][$j] ne "na") {printf OUTPUT2 "%.6f\t",$sorteds_array[$i][$j];}
            else {print OUTPUT2 "na\t";}
            if($sorteds_array[$i][$j] ne "na") {
                if($sorteds_array[$i][$j] != 0) {
                    $Sum[$i] += -2*log($sorteds_array[$i][$j]);
                }
                else {$Sum[$i] += -2*log(1/($size+1));}
            }
            else {$cna += 1;}
        }
        if($cna == $#{$sorteds_array[$i]} + 1) {
            $Sum[$i] = "na";
        }
        $value_ant=$sorteds_array[$i][$j];
        if($Sum[$i] ne "na") {printf OUTPUT2 "\t%.6f\n",$Sum[$i];}
        else {print OUTPUT2 "na\n";}
    }
    close(OUTPUT2);
}

my $ndgf = 2*($#cols+1);
my $sumlog_obs=0;
my $flag0 = 0;
foreach my $pro_obs (@pro_obs){
  if($pro_obs ne "P1") {
    if($pro_obs > 0) {$sumlog_obs += -2*log($pro_obs); $flag0=1;}
    else {$flag0=0;last;}
  }
}

if($flag0==1) {
    my $sumprob=0;
    my $sumprobp=0;
    my $sumprobn=0;
    my $sumprobe=0;
    if($options{'emp'} eq "y" || $options{'emp'} eq "Y") {
        my $cna = 0; my $sizes;
        for my $i (1 .. $#Sum) {
          if($Sum[$i] ne "na") {
            if($sumlog_obs > $Sum[$i]) {$sumprobp += $filaw[$i];}
            if($sumlog_obs < $Sum[$i]) {$sumprobn += $filaw[$i];}
            if($sumlog_obs == $Sum[$i]) {$sumprobe += $filaw[$i];}
            $sizes += $filaw[$i];
          }
          else {
            $cna += 1;
          }
        }
        $sumprobp /= ($sizes);
        if($sumprobp < 0) {$sumprobp = 0;}
        $sumprobn /= ($sizes);
        if($sumprobn < 0) {$sumprobn = 0;}
        $sumprobe /= ($sizes);
        if($sumprobe < 0) {$sumprobe = 0;}
        #$sumprob = (0.5-abs(0.5-$sumprobp))+(0.5-abs(0.5-$sumprobn));
        #$sumprob += $sumprobe;
        $sumprob = 1-$sumprobp;
    }
    my $chi2=0;
    #calculate chi2 with $#cols df...
    $chi2=Statistics::Distributions::chisqrprob (2*($#cols+1),$sumlog_obs);
    
    print  "\n\ COMPOSITE PROBABILITY:";
    print  "\nP1:\tCst\tProb(Chi2[$ndgf df])";
    if($options{'emp'} eq "y" || $options{'emp'} eq "Y") {print "\tProb(empirical)";}
    print "\n";
    printf "RESULT:\t%.6f\t",$sumlog_obs;
    printf "%.6f\t",$chi2;
    if($options{'emp'} eq "y" || $options{'emp'} eq "Y") {printf "%.6f\t",$sumprob;}

    print  OUTPUT "\nP1:\tCst\tProb(Chi2[$ndgf df])";
    if($options{'emp'} eq "y" || $options{'emp'} eq "Y") {print OUTPUT "\tProb(empirical)";}
    print  OUTPUT "\n";
    printf OUTPUT "RESULT:\t%.6f\t",$sumlog_obs;
    printf OUTPUT "%.6f\t",$chi2;
    if($options{'emp'} eq "y" || $options{'emp'} eq "Y") {printf OUTPUT "%.6f\t",$sumprob;}
}
else {
    print "RESULT:\tinf\t0.000000";
    if($options{'emp'} eq "y" || $options{'emp'} eq "Y") {print "\t0.000000";}
    print OUTPUT "RESULT:\tinf\t0.000000";
    if($options{'emp'} eq "y" || $options{'emp'} eq "Y") {print "\t0.000000";}
}
print "\n";
print OUTPUT "\n";

close(OUTPUT);

print "\nDone\n";

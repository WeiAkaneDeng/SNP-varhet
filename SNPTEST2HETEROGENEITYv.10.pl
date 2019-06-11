###############################################
#
#	SNPTEST2HETEROGENEITY.pl     v.10
#
###############################################
#RESOLVE COMMAND LINE OPTIONS
#############################
while (my $arg = shift) 
{
  if ($arg eq "--sample" || $arg eq "-sample" || $arg eq "--s" || $arg eq "-s"){$fileSample=shift;}
  elsif ($arg eq "--gen" || $arg eq "-gen" || $arg eq "--g" || $arg eq "-g"){$fileGenotype=shift;}
  elsif ($arg eq "--out" || $arg eq "-out" || $arg eq "--o" || $arg eq "-o"){$fileOutput=shift;}
  elsif ($arg eq "--excl" || $arg eq "-excl" || $arg eq "--e" || $arg eq "-e"){$fileExclusion=shift;}
  elsif ($arg eq "--header" || $arg eq "-header" || $arg eq "--h" || $arg eq "-h"){$headerResid=shift;}
  elsif ($arg eq "--strand" || $arg eq "-strand"){$fileStrand=shift;}
  elsif ($arg eq "--build" || $arg eq "-build"){$build=shift;}
  elsif ($arg eq "--help" || $arg eq "--h" || $arg eq "-h"|| $arg eq "-help"){&printHelp();}
  elsif ($arg eq "--chr" || $arg eq "-chr"){$forceChr=shift;}
  else {die "Unknown command line argument $arg. Please type: perl SNPTEST2HETEROGENEITY.pl --help to see command line options\n";}
}

###########################
#CHECK COMMAND LINE OPTIONS
###########################
if ($fileSample eq ""){die "Sample file was not defined as a command line option. Please type: perl SNPTEST2HETEROGENEITY.pl --help to see command line options\n";}
if ($fileGenotype eq ""){die "Genotype file was not defined as a command line option. Please type: perl SNPTEST2HETEROGENEITY.pl --help to see command line options\n";}
if ($fileOutput eq ""){die "Output file was not defined as a command line option. Please type: perl SNPTEST2HETEROGENEITY.pl --help to see command line options\n";}
if ($build eq ""){die "Build of genome was not defined as a command line option. Please type: perl SNPTEST2HETEROGENEITY.pl --help to see command line options\n";}
open S, "$fileSample" or die "Cannot open sample file: $fileSample. Please make sure that the file exists! Exit program!\n";
if ($file =~ /\.gz$/) { open(G, "gunzip -c $fileGenotype |") || die "Cannot open genotype file: $fileGenotype. Please make sure that the file exists! Exit program!\n";}
else { open(G, $fileGenotype) || die "Cannot open genotype file: $fileGenotype. Please make sure that the file exists! Exit program!\n";}
open O, ">$fileOutput" or die "Cannot open output file: $fileOutput. Please make sure that you have writing permissions to current folder! Exit program!\n";
open LOG, ">$fileOutput" . ".log" or die "Cannot open log file! Exit program!\n";
if ($fileStrand ne ""){open STRAND, "$fileStrand" or die "Cannot open strand file: $fileStrand. Please make sure that the file exists! Exit program!\n";}
else{print "No strand file given. All markers are expected to be from positive strand.\n";}
if ($fileExclusion ne ""){open EXCL, "$fileExclusion" or die "Cannot open sample exclusion file: $fileExclusion. Please make sure that the file exists! Exit program!\n";}
else{print "No sample exclusion file given. All samples used in analysis (and expected to have same gender).\n";}
if ($headerResid eq ""){die "No resid_pheno header name defined. Please type: perl SNPTEST2HETEROGENEITY.pl --help to see command line options\n";}
if ($forceChr eq ""){print "No chromosome name given. If first column doesnt contain CHR data, it wil be set missing.\n";}
else {print "Chromosome is set to: $forceChr\n";}

#############################
#PRINT HEADER FOR OUTPUT FILE
#############################
print O "MarkerName\tSTRAND\tBUILD\tCHR\tPOS\tMINOR_ALLELE\tMAJOR_ALLELE\tN\tN_0\tN_1\tN_2\tMAF\t";
print O "RESID_MEAN_PHENO_0\tRESID_MEAN_PHENO_1\tRESID_MEAN_PHENO_2\tRESID_SD_PHENO_0\tRESID_SD_PHENO_1\tRESID_SD_PHENO_2\t";
print O "ABSLT_MEAN_PHENO_0\tABSLT_MEAN_PHENO_1\tABSLT_MEAN_PHENO_2\tABSLT_SD_PHENO_0\tABSLT_SD_PHENO_1\tABSLT_SD_PHENO_2\t";
print O "HWE_P\tCALL_RATE\tIMPUTED\tINFO_TYPE\tINFO\tBAYES_FACTOR\n";

####################
#READ EXCLUSION LIST
####################
if ($fileExclusion ne "")
{
	$i=0;
	while(<EXCL>)
	{
		chomp;
		$excluded{$_}=1;
		$i++;
	}
	close EXCL;
	print "Overall $i samples found from the exclusion list.\n";
	print LOG "Overall $i samples found from the exclusion list.\n";
}

#################
#READ STRAND FILE
#################
if ($fileStrand ne "")
{
        $i=0;
        while(<STRAND>)
        {
		s/\r\n/\n/g;    #this will hopefully get rid of dos line endings
	        s/\r/\n/g;      #this will hopefully get rid of mac line endings
                chomp;
		($a,$b)= split(/\s+/);
                if ($b eq "-" || $b eq ""){$strandSwitch{$a}=1;$i++;}
                $i++;
        }
        close STRAND;
        print "Overall $i marker found from negative strand list.\n";
	print LOG "Overall $i marker found from negative strand list.\n";
}

##################
#READ SAMPLE FILE
###################
$i=0;
$colResid=-1;
$countNA=$countExcl=$countOK=0;
@array_for_se;
while(<S>)
{
	chomp;	
	@data = split(/\s/);
	if ($i==0)
	{
		for ($j=2; $j<scalar(@data);$j++)
		{
			if ($data[$j] eq $headerResid){$colResid=$j;}
		}
		if ($colResid==-1){print LOG "Cannot find $headerResid from sample file. Exit program.\n"; die "Cannot find $headerResid from sample file. Exit program.\n";}
	}	
	elsif($i>1)
	{
		$ind[$i-2] = $data[0];
		$valResid[$i-2] = $data[$colResid];
		if ($valResid[$i-2] eq "NA")
		{
			$included[$i-2]=0;
			$countNA++;
		}
		elsif ($excluded{$data[0]}==1)
		{
			$included[$i-2]=0;
			$valResid[$i-2] = "NA";$valZscore[$i-2] = "NA";$valAbslt[$i-2] = "NA";
			$countExcl++;
		}
		else {$countOK++;push(@array_for_se, $valResid[$i-2]);push(@statResid, $data[$colResid]); $included[$i-2]=1;}

	}
	$i++;
}

		@valResid = winsor(\@valResid); 

print "Sample file read.\n";
print "Samples with phenotype NA: $countNA\n";
print "Samples excluded: $countExcl\n";
print "Samples in analysis: $countOK\n";
print LOG "Sample file read.\n";
print LOG "Samples with phenotype NA: $countNA\n";
print LOG "Samples excluded: $countExcl\n";
print LOG "Samples in analysis: $countOK\n";
close S;

####################
# Read genotype file
####################
while(<G>)
{
	chomp;
	@data = split(/\s/);
	if ($forceChr ne ""){$chr = $forceChr;}
	else {$chr = $data[0];}
	$pos = $data[2];
	$snpid = $data[1];
	$a1 = $data[3];
	$a2 = $data[4];
	if ($strandSwitch{$snpid}==1){$strand="-";}
	else {$strand="+";}
	@a11_resid=();@a12_resid=();@a22_resid=();
	@a11_abslt=();@a12_abslt=();@a22_abslt=();
	$ind=$missing=$j=0;
	if ($curMarker%1000==0){print "Reading marker: $curMarker\r";}
	for ($i=5; $i<scalar(@data); $i+=3)
	{
		if ($data[$i]>=0.95 && $included[$j]==1)
		{
			if ($valResid[$j] ne "NA")
				{push(@a11_resid,$valResid[$j]);}
		}
		elsif ($data[$i+1]>=0.95 && $included[$j]==1)
		{
			if ($valResid[$j] ne "NA")
                {push(@a12_resid,$valResid[$j]);}
		}
		elsif ($data[$i+2]>=0.95 && $included[$j]==1)
		{
			if ($valResid[$j] ne "NA" && $valZscore[$j] ne "NA" && $valAbslt[$j] ne "NA")
                {push(@a22_resid,$valResid[$j]);}
		}
		elsif($included[$j]==1) {$missing++;}
		$j++;
	}
	$mean11_resid  = $mean11_abslt = $mean12_resid  = $mean12_abslt = $mean22_resid = $mean22_abslt = ".";
	$sd11_resid = $sd11_abslt = $sd12_resid  = $sd12_abslt = $sd22_resid = $sd22_abslt = ".";

    if (scalar(@a11_resid)>0)
	{
		$mean11_resid = mean(\@a11_resid);$sd11_resid = sd(\@a11_resid);
        if ($sd11_resid>0)
        {
            for ($i=0; $i<scalar(@a11_resid); $i++)
            {
                push(@a11_abslt, abs(($a11_resid[$i]-$mean11_resid)/$sd11_resid));
            }
            $mean11_abslt = mean(\@a11_abslt);$sd11_abslt = sd(\@a11_abslt);
        }
	}
	if (scalar(@a12_resid)>0)
    {
		$mean12_resid = mean(\@a12_resid);$sd12_resid = sd(\@a12_resid);
        if ($sd12_resid>0)
        {
            for ($i=0; $i<scalar(@a12_resid); $i++)
            {
                push(@a12_abslt, abs(($a12_resid[$i]-$mean12_resid)/$sd12_resid));
            }
            $mean12_abslt = mean(\@a12_abslt);$sd12_abslt = sd(\@a12_abslt);
        }
	}
	if (scalar(@a22_resid)>0)
    {
		$mean22_resid = mean(\@a22_resid);$sd22_resid = sd(\@a22_resid);
        if ($sd22_resid>0)
        {
            for ($i=0; $i<scalar(@a22_resid); $i++)
            {
                push(@a22_abslt, abs(($a22_resid[$i]-$mean22_resid)/$sd22_resid));
            }
            $mean22_abslt = mean(\@a22_abslt);$sd22_abslt = sd(\@a22_abslt);
        }
	}

        
    
	print O "$snpid\t$strand\t$build\t$chr\t$pos\t";
	if (scalar(@a11_abslt)>scalar(@a22_abslt))
	{
		print O "$a2\t$a1\t";
		printf O "%d\t", scalar(@a11_abslt) + scalar(@a12_abslt) + scalar(@a22_abslt);
		printf O "%d\t%d\t%d\t", scalar(@a11_abslt), scalar(@a12_abslt), scalar(@a22_abslt);
		printf O "%G\t", (2*scalar(@a22_abslt)+scalar(@a12_abslt))/(2*(scalar(@a11_abslt) + scalar(@a12_abslt) + scalar(@a22_abslt)));
		printf O "%G\t%G\t%G\t", $mean11_resid, $mean12_resid, $mean22_resid;
		printf O "%G\t%G\t%G\t", $sd11_resid, $sd12_resid, $sd22_resid;
		printf O "%G\t%G\t%G\t", $mean11_abslt, $mean12_abslt, $mean22_abslt;
		printf O "%G\t%G\t%G\t", $sd11_abslt, $sd12_abslt, $sd22_abslt;
	}
	else
	{
		print O "$a1\t$a2\t";
                printf O "%d\t", scalar(@a22_abslt) + scalar(@a12_abslt) + scalar(@a11_abslt);
                printf O "%d\t%d\t%d\t", scalar(@a22_abslt), scalar(@a12_abslt), scalar(@a11_abslt);
                if (scalar(@a11_abslt) + scalar(@a12_abslt) + scalar(@a22_abslt)>0){
		    printf O "%G\t", (2*scalar(@a11_abslt)+scalar(@a12_abslt))/(2*(scalar(@a11_abslt) + scalar(@a12_abslt) + scalar(@a22_abslt)));}
	        else {print O ".\t";}
                printf O "%G\t%G\t%G\t", $mean22_resid, $mean12_resid, $mean11_resid;
                printf O "%G\t%G\t%G\t", $sd22_resid, $sd12_resid, $sd11_resid;
                printf O "%G\t%G\t%G\t", $mean22_abslt, $mean12_abslt, $mean11_abslt;
                printf O "%G\t%G\t%G\t", $sd22_abslt, $sd12_abslt, $sd11_abslt;
	}
	if ($data[0] eq "---")
	{
		print O ".\t.\t1\t.\t.\t.\n";
	}
	else
	{
#		printf O "%.4E\t", hwe(scalar(@a11_abslt), scalar(@a12_abslt), scalar(@a22_abslt));
		print O ".\t";
		printf O "%.3f\t", ((scalar(@a22_abslt) + scalar(@a12_abslt) + scalar(@a11_abslt))/(scalar(@a22_abslt) + scalar(@a12_abslt) + scalar(@a11_abslt) + $missing));
		print O "0\t0\t.\t.\n";
	}
	$curMarker++;
}


print LOG "\nMean\t", mean(\@statResid);
print LOG "\nSD\t", sd(\@statResid);
print LOG "\nSkewness\t", skewness(\@statResid);
print LOG "\nKurtosis\t", kurtosis(\@statResid);

@statResid = winsor(\@statResid);
print LOG "\nTransformed Mean\t", mean(\@statResid);
print LOG "\nTransformed SD\t", sd(\@statResid);
print LOG "\nSkewness\t", skewness(\@statResid);
print LOG "\nKurtosis\t", kurtosis(\@statResid);

print "\nAnalysis finished\n";
print LOG "\nAnalysis finished\n";


#sub 
#zScore
#{
#
#$meanZ = mean(\@array_for_se);
#$sdZ = sd(\@array_for_se);
#for ($i=0; $i<scalar(@valResid); $i++)
#{
#        $valZscore[$i] = ($valResid[$i]-$meanZ)/$sdZ;
#        $valAbslt[$i] = abs($valZscore[$i]);
# }
#}

###############
#HELP PRINT SUB
###############
sub 
printHelp
{
	print "\nScript SNPTEST2HETEROGENEITY.pl has following command line options:\n";
	print "\n\t--sample\tSample file name (mandatory)\n";
	print "\t--gen\tGenotype file name (mandatory)\n";
	print "\t--out\tOutput file name (mandatory)\n";
	print "\t--strand\tStrand file name (optional, by default all markers are + strand). \n\t\tFile has either two columns: markername strand; or one column, if only \n\t\tnegative trand markers are listed. Each marker, not present in strand file will be set to + strand\n";
	print "\t--build\tBuild of the genome (mandatory)\n";
        print "\t--excl\tExclusion file name (optional, defines the sample ID's, which will \n\t\tbe removed from analysis. File has one column with sample ID's)\n";
	print "\t--header\tColumn header in sample file (mandatory, define header name for residual of phenotype)\n";
	print "\t--chr\tChromosome number (optional, if no chromosome number defined, chromosome \n\t\tnumber is searched from first column and if it is \"---\" then left missing)\n\n";
	exit(0);
}

##############
# MEAN SUB
##############
sub
mean
{
  @_ == 1 or die ('Sub usage: $mean = mean(\@array);');
  my ($array_ref) = @_;
  my $count = scalar @$array_ref;
  my @array = @$array_ref;	
  my $sum = 0;
  for (my $x = 0; $x<$count; $x++)
  {
    $sum+=$array[$x];
  }
  return ($sum/$count);
}

#############
#  SD SUB
############
sub
sd
{
  @_ == 1 or die ('Sub usage: $sd = sd(\@array);');
  my ($array_ref) = @_;
  my $count = scalar @$array_ref;
  my @array = @$array_ref;
  my $sum = 0;
  for (my $x = 0; $x<$count; $x++)
  {
    $sum+=$array[$x];
  }
  my $mean = ($sum/$count);
  my $sum = 0;
  for (my $x = 0; $x<$count; $x++)
  {
    $sum += (($mean - $array[$x])*($mean - $array[$x]));
  }
    if ($count>1){return sqrt($sum/($count-1));}
    return 0;
}


#############
#  skewness SUB
############
sub
skewness
{
    @_ == 1 or die ('Sub usage: $skewness = skewness(\@array);');
	my ($array_ref) = @_;
	my $count = scalar @$array_ref;
	my @array = @$array_ref;
	my $sum = 0;
	for (my $x = 0; $x<$count; $x++)
  {
    $sum+=$array[$x];
  }
  my $mean = ($sum/$count);
  my $sum1 = 0;
  my $sum2 = 0;
  for (my $x = 0; $x<$count; $x++)
  {
    $sum1 += (($array[$x]-$mean)**2);
	$sum2 += (($array[$x]-$mean)**3);

  }
    if ($count>1){return ($sum2/(($count-1)*((sqrt($sum1/($count-1)))**3)));}
    return 0;
}
 
 
#############
#  kurtosis SUB
############

sub kurtosis {
    @_ == 1 or die ('Sub usage: $kurtosis = kurtosis(\@array);');
	my ($array_ref) = @_;
	my $count = scalar @$array_ref;
	my @array = @$array_ref;
	my $sum = 0;
	for (my $x = 0; $x<$count; $x++)
  {
    $sum+=$array[$x];
  }
  my $mean = ($sum/$count);
  my $sum3 = 0;
  my $sum4 = 0;
  for (my $x = 0; $x<$count; $x++)
  {
    $sum3 += (($array[$x]-$mean)**2);
	$sum4 += (($array[$x]-$mean)**4);
	
  }	
    if ($count>1){return (($count-1)*$sum4/(($sum3)**2));}
    return 0;
}

###########
##  winsor SUB
###########

sub
winsor()
{
 my ($y) = @_;
 my @a = @$y;
 my $n = scalar @a;
 my $mean1 = mean(\@a);
 my $sd1 = sd(\@a);
 my $winsor_top = $mean1 + 3*$sd1;
 my $winsor_bot = $mean1 - 3*$sd1;
 
 for(my $i=0;$i<$n;$i++)
 {
  if ( $a[$i] eq "NA") 
  {}
  elsif ( $a[$i] >= $winsor_top )
  {
   $a[$i] = $winsor_top ;
  }
  elsif ( $a[$i] <=  $winsor_bot )
  {
   $a[$i] = $winsor_bot ;
  }
  }
 return @a;
}

# SNP-varhet
Updated August 2018

###### Wei Q. Deng (<deng@utstat.toronto.edu>) and Reedik Magi (<reedikm@gmail.com>)

The software includes perl scripts to perform genome-wide analysis and meta analysis of SNP exhibiting heterogeneity in phenotypic variances in the current release.  The perl scripts are used to obtain summary statistics in individual study centres while the R scripts are used for meta-analyzing the results centrally. Acceptable data formats include TPED and SNPTEST, which could be generated using PLINK. The publication with derivation can be found here (<https://www.nature.com/articles/ejhg2013166>).

#### For individual participating cohorts in

**TPED format**, please use: 
TPED2HETEROGENEITYv.4.pl

**SNPTEST format**, please use: 
SNPTEST2HETEROGENEITYv.10.pl

The inputs required are the phenotype file, indication of genotype file format, necessary genotype file in either TPED or SNPTEST format, the genome build, the name for the trait of interest from the phenotype file, and a name for the output file.

#### List of Options:

> --sample
Sample file name (mandatory)
Sample file must be in SNPTESTv.2 format with missing data coded as NA (http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format.html ).


> --chr      
Chromosome number (optional, if no chromosome number was defined, chromosome number is searched from first column and if it is "---" then left missing)

>--gen                
Genotype file name (mandatory).
Genotype file in SNPTESTv.2 as described in link above.

>--out
Output file name (mandatory)
See Data File General Formatting Recommendations section above for detailed instructions on files names.

> --strand
Strand file (optional, if this option is not included in the analysis, by default all markers are assigned + strand. Please flip all markers to + strand before using script in the software of your chose (PLINK, SNPTEST, etc.). 
File has either two columns: markername strand; or one column, listing only negative strand markers, where markers absent in the strand file will be set to + strand.


> --snp
Select specific SNPs for analysis

>--pheno
A phenotype file with FID, IID columns and the phenotype variable of interest. The file should be tab delimited. 

>--build
Build of the genome (mandatory)

>--excl
Exclusion file name (optional, defines the sample ID's, which will be removed from analysis. File has one column with sample ID's). This option is useful if data contain both male and female samples and you want to exclude all samples of one sex for sex-specific analysis.


>--header    
Column header in sample file (mandatory, define header name for residual of phenotype). This option defines the name of analyzed phenotype as it is given in the first row of sample file. !!!Missing data must be coded as “NA”!!!


>--help
Invoke the help file




#### An example of typical use for analysis of variance heterogeneity for BMI:

```
Software name \
--pheno Phenotype_file \
--tped Genotype_file_TPED_format \
--tfam FAM_file \
--build 36 \
--header BMI \
--out output
```

Note that we apply a filtering in the analysis to ensure all possible genotype groups are represented and the variance heterogeneity results are valid: lowest number of genotype count has to be greater than 20.


In the output (as txt.gz files), the following are returned:
```
MarkerName: rs SNP names as in provided by the genotype file
CHR
POS
MINOR_ALLELE
N
N_0
N_1
N_2
MAF 
RESID_SD_PHENO_0
RESID_SD_PHENO_1
RESID_SD_PHENO_2
Z_0   =  ABSLT_MEAN_PHENO_0* RESID_SD_PHENO_0
Z_1
Z_2
Z_VAR_0 =  (ABSLT_SD_PHENO_0* RESID_SD_PHENO_0)^2
Z_VAR_1
Z_VAR_2
```
In addition, variance heterogeneity test specific for sex chromosome is currently under development and will be available in the next release. The raw R script is available in an R package here <https://github.com/WeiAkaneDeng/Xvarhet>. 








### Meta-analysis centre:

The outputs from the individual participating cohort can be combined directly under the MLevene program.

1. Make an input directory and an input text file (meta_input.txt) that contains names of the output files from cohorts:

e.g. they might look something like this:

BMI_outfile_cohort1.txt.gz
BMI_outfile_cohort2.txt.gz
BMI_outfile_cohort3.txt.gz

2. Run the MLevene function to get results

```
MLevene_C="./MLEVINE"
$MLevene_C -i meta_input.txt -o outfile.txt &
```

The output text file will look like this, with column Levene_MetaPval denoting the meta-analyzed Levene’s test p-value.

- Name: marker name

- EA: Effective allele	

- NEA: Non-effective allele	

- N: total sample size	

- N_0: aggregated major allele homozygote count	

- N_1: aggregated heterozygote count	

- N_2: aggregated minor allele homozygote count	

- CohortCount: total number of cohort sites included	

- MetaStatN: numerator

- MetaStatD: denominator	

- MetaLevene: statistics	

- P: meta Levene’s p-value	

- Pressence: indicating the SNP is present in which site with respect to orders in the input file.





#### Toy data example:

The toy data in TPED file format is included in *./Stage1TestingFiles/* (with 6 study sites providing individual-level data).

In bash, the following can be used to obtain site-specific analysis results:

```
for study in {1..6}
do
perl TPED2HETEROGENEITYv.4.pl --pheno ./Stage1TestingFiles/Study${study}/BMI_${study}.txt  --tped ./Stage1TestingFiles/Study${study}/marker_${study}.tped --tfam ./Stage1TestingFiles/Study${study}/tped_test${study}.tfam --header bmi_lg --build 36 --out ./Stage1TestingOutputs/BMI_STUDY${study}.txt
done
```

The results (summary statistics) should be in the corresponding *./Stage1TestingOutputs/* folder with a log file summarizing the characteristics of the quantitative trait from that site for any anomalies (severe skewness or kurtosis). 

From this folder, find all files to be meta-analyzed and output to a text file in the same directory as the c software:

```
ls Stage1TestingOutputs/*.txt > BMI_Levene_input.txt

MLevene_C="./mlevine_1.0.2/MLEVINE" 
$MLevene_C –i BMI_Levene_input.txt -o BMI_Levene_out.txt
```

The final outputs are in *BMI_Levene_out.txt*, with the p-value.

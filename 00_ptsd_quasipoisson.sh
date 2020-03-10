 awk '{print $1,$2,$16,$17,$18,$19,$20,$21,$22,$24,$25,$14}' UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno > UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.allcov
 awk '{print $1,$2,$16,$17,$18,$19,$20,$21,$22,$25}' UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno > UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.cov
  awk '{print $1,$2,$17,$18,$19,$20,$21,$22}' UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno > UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.cov_red
  
 
 
 awk '{if(NR==1) print $1,$2,"PCL_zero"; else print $1,$2,$8 - 6}'  UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno  >  UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.m0
 
 
for chr in {1..22}
do 
 cp ukbb_ptsd_bgen_eur_unrel_10.bim ukbb_ptsd_bgen_eur_unrel_10.bk
 Rscript rename_duplicates.r ukbb_ptsd_bgen_eur_unrel_"$chr".bim
 done
 
 Rscript  Rserve.R
 
 #Make sure that the LAST column of the covariates is the LT count!!!!!
 regression_file=quasipoisson_rserve.r
 
 for chr in {1..22}
 do
   awk '{print $2,$6}' /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bim_rename > /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bim.a2
  done
   
 counter=1
 for chr in {1..22}
 do
 echo  $counter $chr

#Make a file for just ukb snps..
for chr in {1..22}
do
./plink --bed /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bed \
 --bim /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bim_rename \
 --fam /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".fam \
 --extract  w_hm3.noMHC.snplist  \
 --make-bed --out plink_ldsc/ldscsnps_"$chr"  \
 --keep UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.m0 
done

#Split them into segments of..?


for chr in {1..22}
do
 time ./plink --bed /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bed \
 --bim /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bim_rename \
 --fam /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".fam \
 --extract  w_hm3.noMHC.snplist  \
 --geno 0.05 \
  --a2-allele /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bim.a2  --out plink_quasipoisson/ldscsnps_"$chr"  \
  --pheno UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.m0 \
  --covar UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.cov_red --R $regression_file &
  sleep 1200 #sleep between each initiation to load the data into ram and calculate allele frequencies
  
  counter=$((counter+1))
  if [ $counter -eq 7 ]
  then
   echo "Waiting for existing jobs to finish"
   wait
   
   counter=1
  fi
 done
 
 #chr 22 test
 cat plink_ldsc/ldscsnps_"$chr".bim | awk '{print $2}' > chr22.snplist 
  split -n 6 chr22.snplist -d 
  
  for xi in {1..5}
  do
./plink --bfile plink_ldsc/ldscsnps_"$chr"  --geno 0.05 --extract x"$xi" \
 --pheno UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.m0 --covar UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.cov \
 --a2-allele /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bim.a2  \
 --R quasipoisson_rserve.r --out plink_quasipoisson/ldscsnps_x0"$xi"_"$chr" &
 
 
 --exclude plink_test/dupvars_22.dupvar
 for chr in  {1..22}
 do
  LC_ALL=C join -1 1 -2 2 <(awk '{print $2,$6}' bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bim | LC_ALL=C sort -k1b,1) <(LC_ALL=C sort -k2b,2 plink_test/ldscsnps_"$chr".auto.R) | awk '{if (NR==1) print "SNP","A2","CHR","BP","A1","B","SE","Z","P","BW2","SEW2","ZW2","PW2"; else print}' > plink_test/jointest_chr"$chr".rm
 done
 cat plink_test/jointest_chr*.rm | awk '{if (NR==1 || $1 != "SNP") print}' > plink_test/jointest_allchr.rm
 
 conda activate ldsc

    #Weight 1 
    /mnt/ukbb/royce/ldsc/munge_sumstats.py --sumstats  plink_test/jointest_allchr.rm  --signed-sumstats B,0 --ignore Z --N 132988 --p P --out  jointest_chrw1.rm.munge.gz 
    #Weight 2 
    /mnt/ukbb/royce/ldsc/munge_sumstats.py --sumstats  plink_test/jointest_allchr.rm  --signed-sumstats BW2,0 --ignore Z,P --N 132988 --p PW2 --out  jointest_chrw2.rm.munge.gz 

    #Noweight
    #/mnt/ukbb/royce/ldsc/munge_sumstats.py --sumstats  plink_test/jointest_allchr.rm  --signed-sumstats BRreg,0 --ignore Z,P --N 135000 --p Preg --out  jointest_chrnow.rm.munge.gz 
    
     /mnt/ukbb/royce/ldsc/ldsc.py \
    --h2  jointest_chrw1.rm.munge.gz.sumstats.gz \
    --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
    --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
    --out  jointest_chrw1.rmldsc
    
      /mnt/ukbb/royce/ldsc/ldsc.py \
    --h2  jointest_chrw2.rm.munge.gz.sumstats.gz \
    --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
    --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
    --out  jointest_chrw2.rmldsc

##Make the regular stuff to compare to
cut -d " " -f1-10  UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.cov >   UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.cov2
awk '{print $1,$2,$11}' UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.cov | grep -v NA > UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.cov.keep
 for chr in {1..22}
 do
 echo  $counter $chr
#need the covar file without LT in it!! At the same time, subjects without LT need to be removed!!

 time ./plink2 --bed /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bed \
 --bim /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bim_rename \
 --fam /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".fam \
  --extract  w_hm3.noMHC.snplist  --exclude plink_test/dupvars_22.dupvar --out plink_test/ldscsnps_regular_"$chr"_  \
  --pheno UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.m0 \
  --covar UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.cov2 \
  --keep UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.cov.keep \
  --glm cols=chrom,pos,ref,alt,ax,a1freq,firth,test,nobs,beta,se,ci,tz,p hide-covar \
  --memory 4000 
  
  done
  
  cat plink_test/ldscsnps_regular_*_.PCL_sum_imp.glm.linear | sed 's/#//g' | awk '{if (NR==1 || $1!= "CHROM" ) print}' > plink_test/ldscsnps_allchr_regular.PCL_sum_imp.glm.linear
  
  conda activate ldsc

    #Weight 1 
    /mnt/ukbb/royce/ldsc/munge_sumstats.py --sumstats  plink_test/ldscsnps_allchr_regular.PCL_sum_imp.glm.linear --a1 A1 --a2 AX --snp ID --signed-sumstats BETA,0 --ignore T_STAT --N 132988 --p P --out  jointest_chruw.rm.munge.gz 
   
  
     /mnt/ukbb/royce/ldsc/ldsc.py \
    --h2  jointest_chruw.rm.munge.gz.sumstats.gz \
    --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
    --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
    --out  jointest_chruw.rmldsc
    
  
         /mnt/ukbb/royce/ldsc/ldsc.py \
    --rg  jointest_chrw1.rm.munge.gz.sumstats.gz,jointest_chrnow.rm.munge.gz.sumstats.gz \
    --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
    --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
    --out  jointest_w1_wnow.rmldsc     
   98% correlation between analysis...
   
   
    #using the correct one:
Total Observed scale h2: 0.1328 (0.0376)
Lambda GC: 1.1747
Mean Chi^2: 1.2335
Intercept: 0.9198 (0.0476)

#Using the basic one:

    
 library(data.table)
 errorbars=1
 data <- fread("plink_test/jointest_allchr.rm",data.table=F)
unadj_filtered <- sort(data$P)
UNADJ <- -log(unadj_filtered,10)
QQ <- -log(ppoints(length(UNADJ)),10)
GCfactor= round(median(qchisq(unadj_filtered,1,lower.tail=F),na.rm=T)/.455,3)

par(bty='l')
plot(c(0,max(QQ)), c(0,max(UNADJ)+ 1), xlab='Expected -log10(p)', ylab='Observed -log10(p)', col='white', cex=1.3, cex.axis=1.2, cex.lab=1.5,pch=20)
if(errorbars == 1)
{
    #code for error bars
    ranks <- c(1:length(QQ))
    CIlower <- qbeta(.025, ranks, length(QQ)-ranks +1)
    CIupper <- qbeta(.975, ranks, length(QQ)-ranks +1)
    plotCIlower <- -log(CIlower,10)
    plotCIupper <- -log(CIupper,10)
    segments(x0=QQ,x1=QQ, y0=plotCIlower,y1=plotCIupper,lwd=2,col='grey')
}
abline(0,1,col='red', lwd=2)
points(QQ, UNADJ, ,pch=20,col='blue')
legend('topleft', paste('GC Lambda =', GCfactor),  bty='n', cex=1.5, xjust=1)
 
 unadj_filtered <- sort(data$PW2)
 UNADJ <- -log(unadj_filtered,10)
 QQ <- -log(ppoints(length(UNADJ)),10)
 points(QQ, UNADJ, ,pch=20,col='red')
round(median(qchisq(unadj_filtered,1,lower.tail=F),na.rm=T)/.455,3)

 correlate ranks of significant snps
 dx <- subset(data,P < 5e-3)
 cor.test(dx$P,dx$Preg)
 #I have some weird intercept effect? A set of 100% correlated SNPs then a set of deflated snps. Apparent in whole data
 #Its very important that I note this, its very weird!!
 plot(-log10(dx$P),-log10(dx$Preg))
 #It is inherent in weighting, its not my correction factor
 plot(-log10(data$wrongp),-log10(data$Preg))
  plot(-log10(data$P),-log10(data$Preg))
 
 plot(-log10(data$pw2),-log10(data$Preg))
  
 plot(-log10(data$P),-log10(data$PW2))
 
  ./plink --bfile /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_$chr  --out plink_test/dupvars_"$chr"  --list-duplicate-vars suppress-first

 ./plink --bed /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bed --bim /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bim_rename --fam /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".fam --extract  w_hm3.noMHC.snplist --a2-allele /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bim.a2 --out plink_test/ldscsnps_"$chr"  \
  --pheno UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.m0 \
  --covar UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.cov --R $regression_file


 
 plink_assoc_"$abbr"_"$phenotype"/"$abbr"_"$phenotype"_"$file_name".regression
 
 "$bgn_folder"/"$file_name" 
 
 #CHANGE OUT THE .a2 files on the LISA SERVER!!
 
 #Get the span of all chromosomes in order to get the number of commands i should run...

 
 #15 per node, running 400 snps/min each... 6000 snps min, 360k per hour per node... for 8.5 million snps...24... 24 nodes is what i need to finish this in one hour..
 
 for chr in {1..22}
 do
  start=$(head -n1  bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bim | awk '{print $4}')
  stop=$(tail -n1 bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bim | awk '{print $4}')
  length=$(wc -l bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bim | awk '{print $1}')
  echo $chr $start $stop $length >> chromosomelengths.txt
 done 
 
 ..cat all bim files and look at the span for 360k snps...
  
  #do i want to write chromosomal spans or just write in 
  
  
  
 #adapt for an m6 analysis
 counter=1
 for chr in {1..22}
 do
 echo  $counter $chr

 time ./plink --bed /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bed \
 --bim /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bim_rename \
 --fam /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".fam \
  --extract  w_hm3.noMHC.snplist --a2-allele /mnt/ukbb/adam/ptsd/bgeneur/ukbb_ptsd_bgen_eur_unrel_"$chr".bim.a2 --exclude plink_test/dupvars_22.dupvar --out plink_test/ldscsnps_"$chr"  \
  --pheno UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.m0 \
  --covar UKB_ptsd_eur_unrelated_allm_oct11_2019.pheno.cov --R $regression_file &
  sleep 1500 #sleep between each initiation to load the data into ram and calculate allele frequencies
  
  counter=$((counter+1))
  if [ $counter -eq 7 ]
  then
   echo "Waiting for existing jobs to finish"
   wait
   
   counter=1
  fi
 done

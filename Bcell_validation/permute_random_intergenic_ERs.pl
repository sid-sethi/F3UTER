#!/usr/bin/perl
# wrapper to produce random regions
# Author: Sid

my $outPath = "/Bcells/Intergenic_3" ;
my $tissue = "GCB1" ; # replicate 1
my $erPath = "/GCB1/ERs/Intergenic" ;


if(! -d "$outPath/Permutation") {system "mkdir -m a=rwx $outPath/Permutation" ;}


# produce the regions to exclude: file "regions_to_exclude.txt"
system "R_v3.6.2/bin/Rscript regions_to_exclude_for_permutation.R $erPath $outPath/Permutation $tissue" ;

# convert ER file to bed
system "cut -f1-3 $erPath/3_prime/$tissue.txt | sed 1d > $erPath/3_prime/$tissue.bed" ;

if(! -d "$outPath/Permutation/$tissue\_ShuffledRegions") {system "mkdir -m a=rwx $outPath/Permutation/$tissue\_ShuffledRegions" ;}


for(my $i=1; $i <=10000; $i++){

print "shuffling iteration $i\r" ;

system "/bedtools/bin/bedtools shuffle -i $erPath/3_prime/$tissue.bed -g chromInfo_hg38_v2.txt -chrom -excl $outPath/Permutation/$tissue\_regions_to_exclude.txt -seed $i > $outPath/Permutation/$tissue\_ShuffledRegions/sim.$i.txt" ;

}

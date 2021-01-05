#!/usr/bin/perl
#
#Author: Sid

my $i;

# output directory
my $dir = "" ;


system "mkdir -m a=rwx $dir/Grid_logs" ;

for($i = 1; $i<21; $i++){

  if(! -d "$dir/Rep_$i"){
    system "mkdir -m a=rwx $dir/Rep_$i" ;
  }

  system "qsub -q long -o $dir/Grid_logs/rep_$i.out -e $dir/Grid_logs/rep_$i.log classification_models.pbs -F \"$dir/Rep_$i rep_$i\"" ;


}

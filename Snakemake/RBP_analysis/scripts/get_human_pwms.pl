#!/usr/bin/perl
#
#Author: Sid

use Getopt::Long;
use Cwd ;

my ($path,$help,$pwm,$db) = "" ;
GetOptions(
    'help'   => \$help,
    'pwm=s' => \$pwm ,
    'db=s' => \$db ,
    'path=s' => \$path ,

) or die "\n**********  Incorrect usage!  ***********\nrun with -help option to see the useage\n ";

sub useage { die(qq/
	USAGE : perl <script> <arguments>
	ARGUMENTS :
                    REQUIRED
                    -pwm -> input attract pwm file
                    -db -> attract human db file
                    OPTIONAL
                    -path -> OUTPUT path where the files need to be saved
                    -help -> prints this help message
               \n/);
}

if($help) { &useage ;}
#if() { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &useage ;}
#unless( ) { print "\n \n" ; &useage ;}
if( !$path || $path =~ /\./ ) { $cwd = getcwd ; $path = "$cwd" ; }

open (hm, $db)or die "Cannot open file $db \n";
my %hash ;
foreach(<hm>){
	$_ =~ s/\n|\r|\s+$//g ;
  if($_ =~ /^Gene_name/){next ;}
	my(@cols) = split("\t",$_) ;
  my $motif_id = $cols[11] ; my $gene_name = $cols[0] ; my $gene_id = $cols[1] ;
	$hash{$motif_id} = "$gene_name:$gene_id:$motif_id" ;
}
close hm ;


$/ = "\n>" ;
open (in, $pwm)or die "Cannot open file $pwm \n";
open(out, ">$path/pwm_human.txt");

foreach(<in>) {
	$_ =~ s/\r|\s+$//g ;
  my ($id,$matrix) = split("\n",$_, 2) ;
  $id =~ s/>//;
  $matrix =~ s/\n>//;
  my($motif_id, $motif_len) = split("\t", $id) ;

	if(exists $hash{$motif_id}) {
	 	print out ">$hash{$motif_id}\t$motif_id\n$matrix\n" ;
	 }
}

close in ; close out ;

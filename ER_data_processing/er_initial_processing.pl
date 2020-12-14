#Author: Sid
# Initial processing of ER files
# Adding gene strand and biotype info

use Getopt::Long;
use Cwd ;

my ($path,$help,$in, $dir) = "" ;
GetOptions(
    'help'   => \$help,
    'i=s' => \$in ,
    'p=s' => \$path ,
    'd=s' => \$dir ,

) or die "\n**********  Incorrect usage!  ***********\nrun with -help option to see the useage\n ";

sub useage { die(qq/
	USAGE : perl <script> <arguments>
	ARGUMENTS :
                    REQUIRED
                    i -> input ER files per tissue
                    OR
                    d -> directory containing all ER files
                    OPTIONAL
                    -p -> OUTPUT path where the files need to be saved [default = cwd ]
                    -help -> prints this help message
               \n/);
}

if($help) { &useage ;}
#if() { print "\n MISSING ARGUMENTS : Give all the required options\n" ; &useage ;}
#unless( ) { print "\n \n" ; &useage ;}
if( !$path || $path =~ /\./ ) { $cwd = getcwd ; $path = "$cwd" ; }
if( $dir =~ /\./ ) { $cwd = getcwd ; $dir = "$cwd" ; }

my %geneStrand ; my %geneType ;
open (info,"</tmp_mnt/filer1/bioinformatics/Sid/Annotation/gene_info_v92.txt")or die "Cannot open file /tmp_mnt/filer1/bioinformatics/Sid/Annotation/gene_info_v92.txt \n";
foreach(<info>){
  $_ =~ s/\n|\r|\s+$//g ;
  if($_ =~ /^Gene.*/){next;}
  my($id,$name,$chr,$start,$end,$str,$type) = split("\t",$_);
  my $strand ;
  if($str == 1){$strand = "+" ;}
  else{$strand = "-" ;}
  $geneStrand{$id} = $strand ;
  $geneType{$id} = $type ;
}
close info ;

if($in) {

  $out = $in ;
  $out =~ s/\.txt/_processed.txt/ ;

  open (in, $in)or die "Cannot open file $in \n";
  open(RES, ">$path/$out");
  foreach my $line(<in>){
    $line =~ s/\n|\r|\s+$//g ;
    if($line =~ /^seqname.*/){ print RES "$line\toverlap_any_gene_v92_strand\toverlap_any_gene_v92_biotype\tuniq_genes_split_read_annot_strand\tuniq_genes_split_read_annot_biotype\tnearest_any_gene_v92_name_strand\tnearest_any_gene_v92_name_biotype\n"; next;}
    my @elements  = split("\t",$line) ;
    my $gene_name = $elements[9]; my $genes_sr = $elements[28] ; # genes_name = overlapping genes, genes_sr = associated via split read
    my $gene_nearest = $elements[11] ; # nearest_any_gene_v92_name

    my ($out1,$out2) = &extract_gene_info($gene_name) ;
    my ($out3,$out4) = &extract_gene_info($genes_sr) ;
    my ($out5,$out6) = &extract_gene_info($gene_nearest) ;

    print RES "$line\t$out1\t$out2\t$out3\t$out4\t$out5\t$out6\n" ;
  }
  close RES;
}


if($dir) {

  opendir(DIR, $dir) or die $!;
  my @files = grep { $_ ne '.' && $_ ne '..' } readdir DIR;

  foreach my $in(@files){

    $out = $in ;
    $out =~ s/\.txt/_processed.txt/ ;

    print "processing $in\n" ;

    open (in, "$dir/$in")or die "Cannot open file $dir/$in \n";
    open(RES, ">$path/$out");
    foreach my $line(<in>){
      $line =~ s/\n|\r|\s+$//g ;
      if($line =~ /^seqname.*/){ print RES "$line\toverlap_any_gene_v92_strand\toverlap_any_gene_v92_biotype\tuniq_genes_split_read_annot_strand\tuniq_genes_split_read_annot_biotype\tnearest_any_gene_v92_name_strand\tnearest_any_gene_v92_name_biotype\n"; next;}
      my @elements  = split("\t",$line) ;
      my $gene_name = $elements[9]; my $genes_sr = $elements[28] ; # genes_name = overlapping genes, genes_sr = associated via split read
      my $gene_nearest = $elements[11] ; # nearest_any_gene_v92_name

      my ($out1,$out2) = &extract_gene_info($gene_name) ;
      my ($out3,$out4) = &extract_gene_info($genes_sr) ;
      my ($out5,$out6) = &extract_gene_info($gene_nearest) ;

      print RES "$line\t$out1\t$out2\t$out3\t$out4\t$out5\t$out6\n" ;
    }
    close RES;
  }
}




#########################################################
## function to extract info from the annotation hash ###
sub extract_gene_info {
  my $gene_string = shift ;

  my @output1; my @output2 ;
  my @ids = split(",",$gene_string) ;

  foreach my $id(@ids){
    my $strandValue; my $typeValue ;
    $id =~ s/^\s+|\n|\r|\s+$//g ;
    if($id eq "NA"){$strandValue = "NA" ; $typeValue = "NA" ;}
    else{
      if(exists $geneStrand{$id}){
        $strandValue = $geneStrand{$id} ;
        $typeValue = $geneType{$id} ;
      }
      else {$strandValue = "NOT_FOUND" ; $typeValue = "NOT_FOUND" ; }
    }
    push(@output1, $strandValue) ;
    push(@output2, $typeValue) ;
  }

  my $out1 = join(",", @output1);
  my $out2 = join(",", @output2);

  return ($out1, $out2) ;
}
################################################

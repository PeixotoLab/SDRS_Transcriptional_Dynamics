use strict;

my $Annotations = "Mart_Ensembl_105.txt";
my $Control = "Controls.txt";
my $outfile = "Controls_Annotated.txt";


open IN,"<$Annotations" or die $!;
my %Annotate;
while (<IN>){
  
  chomp;
  
  if (/ENSMUSG(\d+)./){
    
    $Annotate{$1}=$_;
  }
}

close IN;

open OUT, ">$outfile" or die $!;
open IN2,"<$Control" or die $!;

while (<IN2>){
  
  if (/ENSMUSG(\d+)./){
    
    my $ID=$1;
    
    my $match=$Annotate{$ID};
    if (defined $match){
      print "$ID\n";
      
      print OUT "$match$_";
      
    }
    
  }
}
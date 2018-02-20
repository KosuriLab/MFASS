#3 bases in exon, 6 bases in intron
use strict;
$|++;

my $inputfile = $ARGV[0];
my $upstr = $ARGV[1] || 3;
my $downstr = $ARGV[2] || 6;


my $usemaxent = 1;

my %me2x5 = &makescorematrix('me2x5');
my %seq = &makesequencematrix('splicemodels/splice5sequences');

my %bgd;
$bgd{'A'} = 0.27;
$bgd{'C'} = 0.23;
$bgd{'G'} = 0.23;
$bgd{'T'} = 0.27; 

open (FILE,"<$inputfile") || die "can't open!\n";

while(<FILE>) {
    chomp;
    if (/^\s*$/) { #discard blank lines;
        print "\n";
	    next;
    } 
    elsif (/^>/) {
        print "\n";
	    next;
    }
    elsif (/[NQWERYUIOPLKJHFDSZXVBM]/) { #ignore if N or protein
        print "\n";
	    next;
    }
    else {
        $_ =~ s/\cM//g; #gets rid of carriage return
	    my $str = $_;
	    $str = uc($str);
	    
	    my $discard_upstr = $upstr - 3;
	    my $discard_downstr = $downstr - 6;
	
	    $str =~ /\w{$discard_upstr}(\w{9})\w{$discard_downstr}/;
	    my $motifstr = $1;
	
	    if (length($motifstr) != 9) {
	        die("9 mer was not found inside of string '$str' with $upstr and $downstr upstream and downstream!");
	    } 
	
	    if ($usemaxent) { 
	        print sprintf("%.2f",&log2(&scoreconsensus($motifstr)*$me2x5{$seq{&getrest($motifstr)}}))."\n"; 	
	    }
    }
}

  
sub makesequencematrix{
    my $file = shift;
    my %matrix;my $n=0;
    open(SCOREF, $file) || die "Can't open $file!\n";
    while(<SCOREF>) { 
	chomp;
	$_=~ s/\s//;
	$matrix{$_} = $n;
	$n++;
    }
    close(SCOREF);
    return %matrix;
}
sub makescorematrix{
    my $file = shift;
    my %matrix;my $n=0;
    open(SCOREF, $file) || die "Can't open $file!\n";
    while(<SCOREF>) { 
	chomp;
	$_=~ s/\s//;
	$matrix{$n} = $_;
	$n++;
    }
    close(SCOREF);
    return %matrix;
}

sub getrest{
  my $seq = shift;
  my @seqa = split(//,uc($seq));
  return $seqa[0].$seqa[1].$seqa[2].$seqa[5].$seqa[6].$seqa[7].$seqa[8];
}
sub scoreconsensus{
  my $seq = shift;
  my @seqa = split(//,uc($seq));
  my %bgd; 
  $bgd{'A'} = 0.27; 
  $bgd{'C'} = 0.23; 
  $bgd{'G'} = 0.23; 
  $bgd{'T'} = 0.27;  
  my %cons1;
  $cons1{'A'} = 0.004;
  $cons1{'C'} = 0.0032;
  $cons1{'G'} = 0.9896;
  $cons1{'T'} = 0.0032;
  my %cons2;
  $cons2{'A'} = 0.0034; 
  $cons2{'C'} = 0.0039; 
  $cons2{'G'} = 0.0042; 
  $cons2{'T'} = 0.9884;
  if (($bgd{$seqa[3]}*$bgd{$seqa[4]}) == 0) {
      die("Seq was $seq. \nPositions 3 and 4 were $seqa[4] and $seqa[5]. \nDivide by zero.");
  } 
  my $addscore = $cons1{$seqa[3]}*$cons2{$seqa[4]}/($bgd{$seqa[3]}*$bgd{$seqa[4]}); 
  return $addscore;
}

sub log2{
      my ($val) = @_;
      
      return $val > 0 ? log($val)/log(2) : -inf;    
}

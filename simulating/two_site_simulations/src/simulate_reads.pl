use strict ; 
use warnings ; 
use Math::Random ; 
use List::Util ;

### cmd line is [ploidy] [depth] [selam] [genotypes]

### set sample ploidy
my $ploidy = $ARGV[0] ;
### set mean per site depth 
my $depth = $ARGV[1] ; 
### files are selam output (uncompressed)
### macs2hmm output is genotypes

### length 
my $length = 1e7 ; 
my $num_panel = 2 ; 

### data objects to store SELAM output 
my %ancestry ; 
my %position ; 

### read SELAM file and record ancestry at relevant positions 
print STDERR "READING ANCESTRY DATA\n" ;
open IN, "<$ARGV[2]" ;
while (<IN>) { 
	if ( $_ =~ m/^#/ ) { 
		next ;
	}
	chomp $_ ; 
	my @split = split (/\t/, $_) ; 
	push @{ $position{$split[3]*2+$split[5]-2} }, int($split[7]*$length/0.427) ;
	push @{ $ancestry{$split[3]*2+$split[5]-2} }, $split[6] ; 
}
close IN ; 

### so we don't have any look ahead errors
foreach my $ind ( keys %position ) { 
	push @{ $position{$ind} }, 1000000000 ; 
}

my %panel ; 
my %rec ; 
my %depth ; 
my %geno_a ; 
my %geno_b ; 
print STDERR "ACQUIRING GENOTYPE DATA\n" ;
while (<STDIN>) { 
	chomp $_ ; 
	my @split = split ( /\t/, $_ ) ; 

	$panel{$split[0]} = $split[1] ; 
	foreach ( 2..4 ) { 
		$panel{$split[0]} = $panel{$split[0]}."\t".$split[$_] ;
	}
	$rec{$split[0]} = $split[$#split] ; 

	print 1, "\t", $split[0], "\t", $panel{$split[0]}, "\t", $rec{$split[0]} ; 

	for ( my $i = 0 ; $i < length( $split[$num_panel*2+1] ) ; $i += $ploidy ) { 
		my $geno = 0 ; 
		for ( my $ind = $i ; $ind < $i + $ploidy ; $ind ++ ) { 
			while ( $split[0] > ${ $position{$ind} }[1] ) { 
				shift ( @{ $position{$ind} } ) ; 
				shift ( @{ $ancestry{$ind} } ) ; 
			}
			$geno += substr( $split[$num_panel*2+1+$ancestry{$ind}[0]], $ind, 1 ) ;   
		}

		my $count = Math::Random::random_poisson( 1, $depth ) ; 
		my $rate = ( $geno/$ploidy ) * 0.99 + ( 1 - $geno/$ploidy ) * 0.01 ; 
		my $a = Math::Random::random_binomial( 1, $count, $rate ) ;  	
		my $b = $count - $a ; 

		print "\t", $a, "\t", $b ; 

	}
	print "\n" ;
}


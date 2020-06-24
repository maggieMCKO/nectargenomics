#!/usr/bin/env perl
use Modern::Perl ;
use File::Slurp ;
use Number::Range ;
use Data::Dumper ;
use Carp ;

my $help="USAGE: $0 [BED FILE] [DB NAME] [CHUNK SIZE] [OVERLAP SIZE] [MAX SIZE] 
	Mandatory arguments:
	BED FILE => Genes to not split inside
	DB NAME  => Species, db code, e.g. HLmyoMyo6

	Optional arguments:
	CHUNK SIZE   => Preferred size for splits in bp (default 2.5Mb)
	OVERLAP SIZE => Amount of overlap in bp (default is 10% of CHUNK SIZE)
	MAX SIZE     => Maximum allowed size of a split (default 1.75 times CHUNK SIZE)
" ;

die "Not enough arguments!\n\n$help" if (scalar(@ARGV)<2) ;

my $genes=$ARGV[0] ; # input bed file with genes not to split inside
my $db=$ARGV[1] ; # db containing chrom.sizes
my $chromFile="$ENV{genomePath}"."/gbdb-HL/$db/chrom.sizes" ;
my @chroms=read_file($chromFile,chomp=>1) ;
my $chunk=$ARGV[2] ; # maximum chunk size allowed, this will be flexible
$chunk=2500000 unless ($chunk) ;
#my $chunk=$max*0.8 ; # starting chunk size to begin search for splits of
my $inc=$chunk*0.005 ; # increment to increase by if no splits available
my $overlap=$ARGV[3] ;
$overlap=($chunk*0.1) unless (defined $overlap) ;
my $max=$ARGV[4] ;
$max=$chunk*1.75 unless (defined $max);# hard maximum for chunk size 

my $btCall="bedSort $genes stdout | bedtools merge" ;
my $splitLst=`$btCall` ;
die "ERROR: bedSort $genes stdout | bedtools merge\n\n$help" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0) ;
my @splits=split(/\n/,$splitLst) ;

foreach my $l (@chroms){
	my($chr,$size)=(split(/\t/,$l))[0,1] ;
	do {print "$chr\t0\t$size\tsize_$size\n";next} if ($size<$chunk) ;
	my @spots=grep /$chr\b/, @splits ;
	my $range=Number::Range->new() ;
	foreach my $spot (@spots){
		my ($left,$right)=(split(/\t/,$spot))[1,2] ;
		$range->addrange("$left..$right") ;			
	}
	my $w=0 ;
	my $plus=0 ;
	for (my $i=0;$i<$size;$i+=$plus){
		my $left=0 ;
		$left=($i-$overlap)if($i>0) ;
		my $right=$left+$chunk ;
		if ($right > $size){
			$right=$size ;
			$w = $right - $left ;
			print "$chr\t$left\t$right\tsize_$w\n" ;
			next
		}
		if ($range->inrange("$right")){
			for ($right=$left+$chunk;$right<$size;$right+=$inc){
				last if ($right > ($left+$max)) ;
				next if ($range->inrange("$right")) ;
				last if (! $range->inrange("$right")) ;
			}
		}
		$right=($left + $chunk) if ($right > ($left+$max)) ;
		$w = $right - $left ;
		$plus=$w-$overlap ;
		if (($range->inrange("$left"))&&($range->inrange("$right"))){
			carp "$chr\t$left\t$right\tsize_$w - interval splits in genic regions"
		} ;
		print "$chr\t$left\t$right\tsize_$w\n" ;
		
	}
}
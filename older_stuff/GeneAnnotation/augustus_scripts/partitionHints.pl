#!/usr/bin/env perl
use Modern::Perl ;
use File::Slurp ;

my $help="USAGE: $0 [BED FILE OF SPLIT REGIONS] [HINT FILE]\n\n" ;
die "Not enough arguments!\n\n$help" if (scalar(@ARGV)<2) ;

my $tmpBed=`mktemp` ;
die 'mktemp failed\n' if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0) ;
chomp $tmpBed ;
my $tmpHints=`mktemp` ;
die 'mktemp failed\n' if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0) ;
chomp $tmpHints ;

`cat $ARGV[0] | sort -n -k2,2 |sort -s -k1,1 > $tmpBed` ;
die 'ERROR: failed "cat $ARGV[0] | sort -n -k2,2 |sort -s -k1,1 > $tmpBed"' if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0) ;
`cat $ARGV[1] | sort -n -k4,4 | sort -s -k1,1 > $tmpHints` ;
die 'ERROR: failed "cat $ARGV[1] | sort -n -k4,4 | sort -s -k1,1 > $tmpHints"' if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0) ;

my $splitDir=$ARGV[2] ;
$splitDir="hints_split" unless ($splitDir) ;
if (!-d $splitDir){
	mkdir $splitDir
}

open my $hintfh, "<$tmpHints" or die "Cannot open $tmpHints\n\n" ;
my @splits=read_file("$tmpBed",chomp=>1) ;

HINT: while (my $l=<$hintfh>){
	my ($scaffold,$s,$e)=(split(/\h+/,$l))[0,3,4] ;
	SPLIT: foreach my $split (@splits){
		my ($chr,$i,$j)=(split(/\h+/,$split))[0,1,2] ;
		if (($s>=$i)&&($e<=$j)&&($scaffold eq $chr)){
			append_file("$splitDir/$chr.$i.$j.split.hints",$l) ;
			next HINT
		}
		elsif($e<$i){
			#append_file("hints_split/$chr.$i.$j.split.hints",$l) ;
			next SPLIT
		}
		elsif(($s<$j)&&($e>$j)){
			#append_file("hints_split/$chr.$i.$j.split.hints",$l) ;
			next SPLIT
		}
		elsif (($s<$i)&&($e>$i)&&($e<$j)){
			next HINT
		}
		elsif($s>$j){
			shift @splits ;
			next HINT
		}
		elsif ($scaffold ne $chr){
			next SPLIT
		}
		else{
			die "ERROR: File should be sorted??\n\nCURRENT SPLIT:$split\nHINT: $l" 
		}
	}
}

`rm $tmpHints $tmpBed` ;

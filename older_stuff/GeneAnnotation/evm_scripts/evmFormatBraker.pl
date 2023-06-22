#!/usr/bin/env perl
use warnings;
use strict ;

my $count ;
my $fh ;
if (scalar(@ARGV)>0){
	open $fh, "<$ARGV[0]";
}
else{
	$fh = *STDIN ;
}

while (<$fh>){
	my $l = $_ ;
	if ($l =~ /^##( |\w)/){
		print $l
	}
	elsif ($l =~ /^###/){
		$count = 0 ;
		print $l 
	}
	elsif($l =~/\tmRNA\t/){
		my ($start,$end) = (split(/\t/,$l))[3,4] ;
		print $l ;
	}
	elsif ($l=~/\tCDS\t/){
		$count ++ ;
		chomp $l ;
		my $e = $l ;
		$e =~ s/CDS/exon/ ;
		my ($id) = $e =~ /transcript_id=([\w.]+)/ ;
		$id=";ID=$id.$count" ;
		print $e.$id."\n";
		print $l.$id."cds\n";
	}
	else{
		print $l ;
	}

}

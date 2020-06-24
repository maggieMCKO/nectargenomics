#!/usr/bin/env perl
use Modern::Perl ;
use List::MoreUtils qw(any) ;

my $fh ;

if (scalar(@ARGV)>0){
	open $fh, "<$ARGV[0]" ;
}
else{
	$fh=*STDIN
}
while (my $line=<$fh>){
	chomp $line ;
	my @info=split(/\t/,$line) ;
	#my $strand=$info[2] ;
	my @exonStarts=split(/,/,$info[8]) ;
	my @exonEnds=split(/,/,$info[9]) ;
	if ($info[2] eq "+"){
		if (any {$_ == $info[6]} @exonEnds){
			print $line."\n"
		}
		elsif ($info[6]==0){
			print $line."\n"
		}		
		else{
			$info[6]+=3 ;
			print join("\t",@info)."\n"
		}
	}
	elsif ($info[2] eq "-"){
		if (any {$_ == $info[5]} @exonStarts){
			print $line."\n"
		}
		elsif ($info[5]==0){
			print $line."\n"
		}
		else{
			$info[5]-=3 ;
			print join("\t",@info)."\n"
		}
	}
}
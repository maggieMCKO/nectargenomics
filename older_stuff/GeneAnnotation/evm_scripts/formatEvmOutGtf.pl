#!/usr/bin/env perl
use Modern::Perl ;

while (my $l=<STDIN>){
	chomp $l ;
	next unless ($l) ;
	print $l if ($l=~ /^#/) ;
	my @info = split(/\t/,$l) ;
	my ($chr,$left)=$info[0]=~/(\w+)\_(\d+)\-\d+/ ;
	if (! $chr){
		print join("\t",@info)."\n" ;
		next
	}
	my $offset=$left-1 ;
	$info[3]+=$offset ;
	$info[4]+=$offset ;
	$info[0]=$chr ;
	print join("\t",@info)."\n"
}



 ##~/bin/EVidenceModeler-1.1.1/EvmUtils/EVM_to_GFF3.pl outs/scaffold_m19_p_92.out  scaffold_m19_p_92 | gt gff3 -tidy -sort | gt gff3_to_gtf | formatEvmOutGtf.pl 

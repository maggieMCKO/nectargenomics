#!/usr/bin/env perl
use warnings;
use strict ;

my @headers ;
my @glines ;
my $fh ;
if (scalar(@ARGV)>0){
	open $fh, "<$ARGV[0]";
}
else{
	$fh = *STDIN ;
}
while (<$fh>){
	my $line = $_ ;
	### get all the seq info and gff version lines
	if ($line =~ /^##( |\w+)/){
		push @headers, $line 
	}
	### ignore comment lines
	elsif ($line =~ /^#[^#]/){
		next
	}
	### get all the annotation lines
	else{
		push @glines, $line
	}
}

my @genes = split(/###\n/,join('',@glines)) ;
if ($glines[scalar(@glines)-1]eq "###"){
	pop @genes ;	
}

my $aln_type;
if ($ARGV[1]){
	$aln_type = "cDNA_match"
}
else{
	$aln_type = "nucleotide_to_protein_match"
}

chomp @headers ;
print join("\n",@headers)."\n" ;

my $count = 0 ;
foreach my $g (@genes){
	$count++ ;
	my @aln = split(/\n/,$g) ;
#	shift @aln ;
	my ($mrna) = grep /\tmRNA\t/, @aln ;
	my $target ;
	if ($mrna =~/Target/){
		($target)=$mrna=~/Target=([\w.\/]+)/ ;
	}
	else{
		my $dodge ;
		($dodge,$target)=$mrna=~/(transcript_id|[Nn]ame)=([\w.\/]+)/
	}
	my $start = 1;
	my $end ;
	foreach my $line(@aln){
		if ($line=~/\texon\t/){
			chomp $line ;
			my @info = split(/\t/,$line) ;
			my ($l,$r)= (@info)[3,4] ;
			$info[2] = "$aln_type" ;
			pop @info ;
			my $diff = $r -$l ;
			$end = $start +$diff ;
			my $t = "ID=CDS$count;Target=$target $start $end\n" ;
			print join("\t",@info)."\t".$t ;
			$start = $end +1 ;
		}
	}
	print "###\n"
}	
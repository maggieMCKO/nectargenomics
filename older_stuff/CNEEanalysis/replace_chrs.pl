#!/bin/env perl

use strict;
use warnings;

my $accfile = shift;
my $bed = shift; 

open my $acc, "$accfile" or die;
my %key;
while (<$acc>) {
	chomp;
	next if /^#/;
	my @fields = split;
        #acc file needs to be in from -> to format
	my $chr = $fields[0];
	$key{$chr} = $fields[1];
#	print "$fields[9]\t\t$fields[6]\n"
}

open my $infile, $bed or die;

while (<$infile>) {
	chomp;
	my @fields = split;
	my $chr = $fields[0];
	my $acc = exists($key{$chr}) ? $key{$chr} : "NA";
	next if $acc eq "NA";
	$fields[0] = $acc;
	print join("\t", @fields);
	print "\n";
} 


#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# usage: perl csv2constraint.pl -i ../data/sequencing_rpoB/raxml_prep/asv_metadata_filtered.csv -c > ../data/sequencing_rpoB/raxml_prep/contree.tre 

# process command line arguments
my $infile;
my $comprehensive;
my $otus = 'tip_label'; # which column has the OTU labels
my $taxon = 'family'; # which column has the higher taxon
GetOptions(
  'infile=s'      => \$infile,
  'otus=s'        => \$otus,
  'taxon=s'       => \$taxon,
  'comprehensive' => \$comprehensive,
);

# initialize data structure
my %constraints = ( 'NA' => [] );

# start reading file
my @header;
open my $fh, '<', $infile or die $!;
LINE: while(<$fh>) {
  chomp;
  
  # read header on first iteration
  if ( not @header ) {
    my @cols = split /,/, $_;
    for my $name ( @cols ) {
      $name =~ s/"//g;
      push @header, $name;
    }
    next LINE;
  }
  
  # read record
  my %record;
  my @record = split /,/, $_;
  for my $i ( 0 .. $#record ) {
    $record[$i] =~ s/"//g;
    $record{$header[$i]} = $record[$i];
  }
  
  # postprocess
  my $name = $record{$taxon};
  $constraints{$name} = [] if not $constraints{$name};
  push @{ $constraints{$name} }, $record{$otus};
}

# print output
print '(';
my @taxa;
TAXON: for my $taxon ( keys %constraints ) {
  my $t;
  if ( $taxon eq 'NA' ) {
    if ( $comprehensive ) {
      $t = join ',', @{ $constraints{$taxon} };
    }
    else {
      next TAXON;
    }
  }
  else {
    if ( scalar(@{ $constraints{$taxon} }) > 1 ) {
      $t = '(' . join(',', @{ $constraints{$taxon} }) . ')';
    }
    else {
      $t = $constraints{$taxon}->[0];
    }
  }
  push @taxa, $t;
}
print join ',', @taxa;
print ');';
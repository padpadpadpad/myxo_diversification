#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Text::CSV 'csv';
use Bio::Phylo::Factory;
use Bio::Phylo::IO 'parse_tree';

# usage: perl csv2constraint.pl -i asv_metadata_filtered.csv -b  backbone.tre > contree.tre 

# process command line arguments
my $infile;
my $backbone;
my $otus = 'tip_label'; # which column has the OTU labels
my $taxon = 'family'; # which column has the higher taxon
GetOptions(
  'infile=s'   => \$infile,
  'otus=s'     => \$otus,
  'taxon=s'    => \$taxon,
  'backbone=s' => \$backbone,
);

# instantiate Factory
my $fac = Bio::Phylo::Factory->new;

# read $backbone
my $tree = parse_tree(
  '-format' => 'newick',
  '-file'   => $backbone,
);

# read $infile
my $csv = csv(
  'in'      => $infile,
  'headers' => 'auto',
);

# expand terminals
for my $tip ( @{ $tree->get_terminals } ) {
  my $name = $tip->get_name;
  for my $otu ( grep { $_->{$taxon} eq $name } @$csv ) {
    $tip->set_child( $fac->create_node( '-name' => $otu->{$otus} ) );
  }
}

# raxml doesn't like nested 1-degree nodes
$tree->remove_unbranched_internals;

# print output
print $tree->to_newick;
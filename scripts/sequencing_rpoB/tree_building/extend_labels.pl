#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Text::CSV 'csv';
use Bio::Phylo::IO 'parse_tree';

# process command line arguments
my $tree_file;
my $otus_file;
my $col_otu   = 'tip_label'; # which column has the OTU labels
my $col_taxon = 'family';    # which column has the higher taxon
GetOptions(
  'tree=s'  => \$tree_file,
  'otus=s'  => \$otus_file,
);

# read tree
my $tree = parse_tree(
  '-format' => 'newick',
  '-file'   => $tree_file,
);

# read table
my $table = csv(
  'in'      => $otus_file,
  'headers' => 'auto',
);

# iterate over OTUs in table
for my $otu ( @$table ) {
  my $family = $otu->{$col_taxon};
  next if $family eq 'NA';
  my $name = $otu->{$col_otu};
  my $node = $tree->get_by_name($name);
  $node->set_name( $name . '_' . $family );
}

print $tree->to_newick;
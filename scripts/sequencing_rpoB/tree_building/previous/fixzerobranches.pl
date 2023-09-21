#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO 'parse_tree';

# process command line arguments
my $intree;
my $add = 0.000001;
GetOptions(
  'intree=s' => \$intree,
  'add=f'    => \$add,
);

# read tree
my $tree = parse_tree( '-format' => 'newick', '-file' => $intree );

# update branches
for my $node ( @{ $tree->get_entities } ) {
  my $bl = $node->get_branch_length; 
  if ( $bl == 0 ) {
    $node->set_branch_length($add);
  }
  else {
    $node->set_branch_length( $bl * 10 );
  }
}

# make ultrametric again
$tree->ultrametricize;

# write output
print $tree->to_newick;
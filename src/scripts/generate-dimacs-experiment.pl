#!/usr/bin/perl

use strict;

if (@ARGV != 2) {
    print STDERR "Usage: generate-scenario.pl [#nodes in graph] [num instances] > result_file\n";
    die;
}

print STDOUT "p aux sp p2p $ARGV[1]\n";
for (my $k = 0; $k < $ARGV[1]; $k++) {
    my $random_idx_1 = int(rand($ARGV[0]));
    my $random_idx_2 = int(rand($ARGV[0]));
    print STDOUT "q $random_idx_1 $random_idx_2\n";
}


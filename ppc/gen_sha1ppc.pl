#!/usr/bin/perl -w

use strict;

my %platform = (
    # Special extra substitutions that have to be done on this platform
    darwin => sub {
        s{sha1_core}{_sha1_core};
        s{;}{\n}g;
    },
);

my $extra = exists $platform{$^O} ? $platform{$^O} : sub {};

while (<>) {
    $extra->();
    s{R\[([^]]+)\]}{'r'.eval"$1"}ge;
    print;
}

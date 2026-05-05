#!/usr/bin/perl -w
use strict;

if (scalar(@ARGV) < 2) {
    print "Usage: <contigs.fasta> <length_threshold> [output_file]\n";
    exit;
}

my $FileName = $ARGV[0];
my $limit    = $ARGV[1] || 1;
my $output   = $ARGV[2];

open(my $IN, '<', $FileName) or die "Couldn't open file $FileName: $!\n";

my $Counter = 0;
my $Length = 0;
my $name = '';
my $str = '';
my $res = '';

while (<$IN>) {
    if (/^>(\S+)/) {
        if ($Length >= $limit) {
            $Counter++;
            print "$name\t$Length\n";
            $res .= $str;
        }
        $str = $_;
        $name = $1;
        $Length = 0;
    } else {
        $str .= $_;
        $Length += length($_);
    }
}

# Handle the last sequence
if ($Length >= $limit) {
    $Counter++;
    print "$name\t$Length\n";
    $res .= $str;
}

close($IN);

if (defined $output) {
    open(my $OUT, '>', $output) or die "Couldn't write to $output: $!\n";
    print $OUT $res;
    close($OUT);
}

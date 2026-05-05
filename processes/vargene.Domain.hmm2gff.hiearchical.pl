#!/usr/bin/perl
use strict;
use warnings;

# Usage: perl hmm2gff.hierarchical.pl LENGTH.Domains.txt < input.domtblout > output.gff

my %domain_lengths;

# Read domain lengths file
open(my $len_fh, '<', $ARGV[0]) or die "Cannot open $ARGV[0]: $!";
while (<$len_fh>) {
    next if /^#/;
    chomp;
    my ($domain, $length) = split(/\s+/);
    $domain_lengths{$domain} = $length;
}
close $len_fh;

my %domain_count;
my $source = "hmmscan";

while (<STDIN>) {
    next if /^#/;
    chomp;
    my @fields = split(/\s+/);

    # Parse domtblout columns
    my $target = $fields[0];        # target name (domain family)
    my $query = $fields[3];         # query name (sequence ID)
    my $ali_from = $fields[15];     # alignment start (1-based)
    my $ali_to = $fields[16];       # alignment end
    my $dom_i_evalue = $fields[11]; # domain i-Evalue
    my $dom_score = $fields[13];    # domain score

    # Standardize domain name
    my $domain = $target;
    $domain = "DBLn" if $domain eq "Duffy_binding_like";

    # Skip if no reference length available
    next unless exists $domain_lengths{$domain};

    my $ref_length = $domain_lengths{$domain};
    my $hit_length = $ali_to - $ali_from + 1;

    # Calculate acceptable range (±15%)
    my $min_len = int($ref_length * 0.80);
    my $max_len = int($ref_length * 1.10);

    # Filter based on domain length
    next if $hit_length < $min_len || $hit_length > $max_len;

    # GFF3 attributes - must include 'label=' for the association pipeline
    my $domain_id = "$domain.$query";
    my $label = $domain;
    my $attributes = "ID=$domain_id;Name=$domain;label=$label;coords=$ali_from-$ali_to;";
    $attributes .= "Target=$domain 1 $ref_length;Evalue=$dom_i_evalue;Score=$dom_score";

    # Print GFF3 line
    print join("\t", 
        $query,          # seqid (query name)
        $source,         # source (hmmscan)
        'domain',        # type
        $ali_from,       # start
        $ali_to,         # end
        $dom_score,      # score
        '.',             # strand
        '.',             # phase
        $attributes
    ), "\n";
}

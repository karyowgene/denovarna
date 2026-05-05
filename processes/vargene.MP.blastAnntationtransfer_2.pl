#!/usr/bin/perl

my $Annotation_f = shift;
my $blast_f = shift;

# Load annotations
open F, $Annotation_f or die "Cannot open annotation file: $!\n";
my %Annotation;
while(<F>){
    chomp;
    my @ar = split(/\t/);
    $Annotation{$ar[0]} = $ar[1];
}
close(F);

# Load BLAST results
open F, $blast_f or die "Cannot open BLAST file: $!\n";
my %blast;
while(<F>){
    chomp;
    my @ar = split(/\t/);
    my $query = $ar[0];
    my $subject = $ar[1];
    
    # Store only the first hit for each query
    if (!defined($blast{$query})) {
        $blast{$query} = $Annotation{$subject};
    }
}
close(F);

# Process FASTA input
while(<STDIN>){
    if (/^>(\S+)/) {
        my $header = $1;
        # Extract just the ID part (before first space)
        my ($id) = split(/ /, $header);
        
        # Remove any trailing description that might have been added previously
        $id =~ s/^(ERR\d+\.g\d+\.t\d+).*/$1/;
        
        if (defined($blast{$id})) {
            print ">$header $blast{$id}\n";
        } else {
            print ">$header\n";
        }
    } else {
        print;
    }
}
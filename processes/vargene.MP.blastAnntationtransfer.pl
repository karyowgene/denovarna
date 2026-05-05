

my $Annotation_f=shift;
my $blast_f=shift;

open F, $Annotation_f or die "Gimme file please!\n";

my %Annotation;
while(<F>){
    chomp;
    my @ar=split(/\t/);
    $Annotation{$ar[0]}=$ar[1];
}
close(F);

open F, $blast_f or die "Gimme file please!\n";

my %blast;
while(<F>){
    chomp;
    my @ar=split(/\t/);
    if (! defined($blast{$ar[0]})){
	$blast{$ar[0]}=$Annotation{$ar[1]};
    }
}
close(F);

while(<STDIN>){
    if (/>(\S+).(cds\d+)/){
	print ">$1.$2 $blast{$1}\n"
    }
    elsif (/>(\S+)/){
	print ">$1 $blast{$1}\n"
    }
    else {
	print
    }

}


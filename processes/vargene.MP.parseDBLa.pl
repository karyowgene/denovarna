
my $debug=0;


while(<>){
    my $oldline=$_;
    my @ar=split(/\s+/);
    if ($ar[0] eq "DBLa") {
	$n=<STDIN>;
	if ($debug){
	    print ">>$oldline";
	    print ">>$n";
	}
	@ar2=split(/\s+/,$n) ;
	if ($ar2[0] ne "DBLa") {
	    print $oldline;
	    print $n
	} elsif ($ar[16]>$ar2[15]) {
	    print $oldline
	}
	else {
	    if ($debug){
		print ">> Do the merge";
	    }
	    $ar[13]+=$ar2[13];
	    $ar[16]=$ar2[16];
	    $ar[18]=$ar2[18];
	    $ar[20]=$ar2[20];
	    print (join("\t",@ar)."\n")}
    } else
    {
	print
    }
} 

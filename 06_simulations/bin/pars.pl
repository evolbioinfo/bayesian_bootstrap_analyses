#!/usr/bin/perl -w

use strict;

my @acgtcols;

my $total=0;
# We select only columns counting A C G or T (not N nor -, etc.)
$_=<>;
chomp();
my @names=split(/\t/);
for(my $i=1;$i<=$#names;$i++){
    if($names[$i] =~/[ACGTacgt]/){
	push @acgtcols, $i;
    }
}
#print join(",",@acgtcols)."\n";
while(<>){
    chomp;

    my @cols = split(/\t/);

    my $site=0;
    for my $i (@acgtcols){
	if($cols[$i]>0){
	    $site++;
	}
    }
    $total+=$site-1;
    #print(($site-1)."\n");
}

print $total."\n";

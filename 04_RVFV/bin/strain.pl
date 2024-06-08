#!/usr/bin/perl
use strict;
use warnings;

sub segment{
    my ($sampname) = @_;
    my $seg = "";
    my $name;
    
    if($sampname =~/segment ([LMS])/){
	$seg = $1;
	$sampname =~s/segment [LMS].*//g;
    }elsif($sampname =~/\s([LMS]) segment/){
	$seg = $1;
	$sampname =~s/\s[LMS] segment.*//g;
    }elsif($sampname =~/nucleocapsid/){
	$seg = "S";
	$sampname =~s/\snucleocapsid.*//g;
    }elsif($sampname =~/non[-]{0,1}structural (protein|gene)/){
	$seg = "S";
	$sampname =~s/\snon[-]{0,1}structural (protein|gene).*//g;
    }elsif($sampname =~/glycoprotein/){
	$seg = "M";
	$sampname =~s/\sglycoprotein.*//g;
    }elsif($sampname =~/polymerase/){
	$seg = "L";
	$sampname =~s/\spolymerase.*//g;
    }elsif($sampname =~/M protein/){
	$seg = "M";
	$sampname =~s/\sM protein.*//g;
    }elsif($sampname =~/NSs (protein|gene)/){
	$seg = "S";
	$sampname =~s/\sNSs protein.*//g;
    }elsif($sampname =~/G(n|2) protein/){
	$seg = "M";
	$sampname =~s/\sG(n|2) protein.*//g;
    }elsif($sampname =~/\(([LMS])\)/){
	$seg = $1;
	$sampname =~s/\s\([LMS]\).*//g;
    }

    $name = $sampname;
    $name =~ s/isolate/strain/g;
    return ($seg,$name);
}

my $name="";
my $id="";
my $segment="";
my $sequence = "";
my $nbsamples = 0;
my %samples;
my $localid = "";


my $ncbifile = $ARGV[0];
my $gbfile = $ARGV[1];


my %meta;
my $acc;
my $sampname = "" ;
my $startsampname = 0;
open(IN,$gbfile);
while(<IN>){
    chomp;
    s/^\s+//g;
    s/\s+$//g;
    
#DEFINITION  Rift Valley fever virus strain Lunyo polymerase (L) gene, complete
#            cds.
#ACCESSION   OQ440142
# 	if(/^>([^\s]*)(.*)segment (.*?),(.*)/){

    # Extract sample name
    if(/^DEFINITION\s+(.*)$/){
	$sampname = $1;
	$startsampname = 1;
    }elsif(/^ACCESSION/){
	$startsampname = 0;
    }elsif($startsampname == 1){
	$sampname .= " ".$_;
    }

    if(/^VERSION\s+(.*)$/){
	# Extract segment from sample name
	print $sampname."\n";
	my ($seg,$name) = segment($sampname);
	print "seg: -".$seg."-\n";
	$acc = $1;
	$meta{$acc}{"strain"} = "";
	$meta{$acc}{"host"} = "";
	$meta{$acc}{"country"} = "";
	$meta{$acc}{"date"} = "";
	$meta{$acc}{"segment"} = $seg;
	# clean name
	$meta{$acc}{"name"} = $name;
    }

    #/mol_type="genomic RNA"
    #/strain="SA01-1322"
    #/host="Aedes vexans arabiensis"
    #/db_xref="taxon:11588"
    #/segment="S"
    #/country="Saudi Arabia"
    #/collection_date="2001"


    
    if(/\/strain=\"(.*)\"/){
	$meta{$acc}{"strain"} = $1;
    }
    if(/\/host=\"(.*)\"/){
	$meta{$acc}{"host"} = $1;
    }
    if(/\/country=\"(.*)\"/){
	$meta{$acc}{"country"} = $1;
    }
    if(/\/collection_date=\"(.*)\"/){
	$meta{$acc}{"date"} = $1;
    }
}
close(IN);

open(IN,$ncbifile);
while(<IN>){
    chomp;
    #if(/^>/){
    #	if(/^>.*(strain|isolate)\s(.*?)\ssegment (.*?),/){
    #	    print "$2\n";
    #	}else{
    #	    #print "not found : $_\n";
    #	}
    #}

    if(/^>/){
	if($sequence ne ""){
	    if($meta{$id}{"segment"} ne ""){
		open(OUT,">>segment_".$meta{$id}{"segment"}.".fasta");
		print OUT ">$id\n";
		print OUT "$sequence\n";
		close(OUT);
		
		open(ID,">>metadata.txt");
		print ID "$id\t$localid\t".$meta{$id}{"name"}."\t".$meta{$id}{"segment"}."\t".$meta{$id}{"strain"}."\t".$meta{$id}{"host"}."\t".$meta{$id}{"country"}."\t".$meta{$id}{"date"}."\n";
		close(ID);
	    }
	}
	if(/^>([^\s]*)/){
	    $id = $1;
	    $name = $meta{$id}{"name"};
	    $sequence = "";
	    if(defined $samples{$name."-".$meta{$id}{"country"}."-".$meta{$id}{"date"}}){
		$localid = $samples{$name."-".$meta{$id}{"country"}."-".$meta{$id}{"date"}};
	    }else{
		$localid = "Sample_".$nbsamples;
		$nbsamples += 1;
		$samples{$name."-".$meta{$id}{"country"}."-".$meta{$id}{"date"}} = $localid;
	    }
	}else{
	    $sequence = "";
	    $name = "";
	    $id = "";
	    $localid = "";
	}
    }else{
	if($name ne ""){
	    $sequence.=$_;
	}
    }
}
if($name ne ""){
    if($meta{$id}{"segment"} ne ""){
	open(OUT,">>segment_".$meta{$id}{"segment"}.".fasta");
	print OUT ">$localid\n";
	print OUT "$sequence\n";
	close(OUT);
	
	open(ID,">>metadata.txt");
	print ID "$id\t$localid\t$name\t".$meta{$id}{"segment"}."\t".$meta{$id}{"strain"}."\t".$meta{$id}{"host"}."\t".$meta{$id}{"country"}."\t".$meta{$id}{"date"}."\n";
	close(ID);
    }
}

close(IN);

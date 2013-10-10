package Packages::NOEs;

require Exporter;
use strict;
use Packages::General qw(GetBondLength);

our (@ISA, @EXPORT, $VERSION, @EXPORT_OK);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(GetNOEs SaveNOEs);
$VERSION = "1.00";

sub GetNOEs {
    my ($nFile, $atoms, $resData, $atomList) = @_;
    my ($atom1, $atom2, %NOE, $noeName, @DATA, $counter, $rec, @tmp, $i, $avg);

    $counter = 0;
    open NOES, $nFile or die "ERROR: Cannot open NOE file $nFile: $!\n";
    while (<NOES>) {
	chomp;
	if ($_ =~/(A|C|G|T)(\d+)\s+(\S+)\s\-\s(A|C|G|T)(\d+)\s+(\S+)(\s*.*)$/) {
	    $atom1 = "${1}${2}_${3}";
	    $atom2 = "${4}${5}_${6}";
	    ($atom1, $atom2) = ($atom2, $atom1) if ($atom1 gt $atom2);
	    if (! isStored(\%NOE, $atom1)) {
		$NOE{$atom1} = getAtomNum($resData, $atom1, $atoms);
	    }
	    if (! isStored(\%NOE, $atom2)) {
		$NOE{$atom2} = getAtomNum($resData, $atom2, $atoms);
	    }

	    $noeName = "${atom1} ${atom2}";
	    $rec = (
		    {
			"NAME"  => $noeName, 
			"ATOM1" => $NOE{$atom1},
			"ATOM2" => $NOE{$atom2},
		    }
		    );

	    if ($7) {
		$avg = $7;
		$i = 1;
		while ($avg =~ /\s+(\d+\.\d+)/g) {
		    $rec->{AVG}{$i} = $1;
		    $i++;
		}
	    }

	    push @DATA, $rec;
	    $counter++;
			   
	}
    }
    close NOES;
    die "ERROR: No valid NOE's found in $nFile\n" if (! @DATA);

    print "found $counter NOEs...";
    return (\@DATA);
}

sub getAtomNum {
    my ($resData, $inStr, $atoms) = @_;
    my ($resNum, $i, $atmNum, $myRes, @tmp, $resName, $atmName);

    if ($inStr =~ /(.+)\_(.+)/) {
	$resName = $1;
	$atmName = $2;
    } else {
	die "ERROR: Invalid string $inStr\n";
    }
    if ($resName =~ /^(\w)(\d+)/) {
	$resNum = $2;
    } else {
	die "ERROR: Invalid NOE specified by RES $resNum\n";
    }
    
    die "ERROR: Residue $resNum not found!\n" if (! exists($resData->{$resNum}));

    @tmp = keys %{ $resData->{$resNum}{ATOMS} };

    for $i (@tmp) {
	if (exists($atoms->{$i}) &&
	    ($atoms->{$i}{ATMNAME} eq $atmName)) {
	    $atmNum = $i;
	    last;
	}
    }

    if (! $atmNum) {
	$atmName = fixName($atmName);
	for $i (@tmp) {
	    if (exists($atoms->{$i}) &&
		$atoms->{$i}{ATMNAME} eq $atmName) {
		$atmNum = $i;
		last;
	    }
	}
    }

    die "ERROR: No atom found named $atmName in residue $resNum\n" if (! $atmNum);

    return $atmNum;
}

sub isStored {
    my ($noes, $searchStr) = @_;
    my ($returnStr) = 0;

    if (exists($noes->{$searchStr})) {
	$returnStr = 1;
    } else {
	$searchStr = fixName($searchStr);
	$returnStr = 1 if (exists($noes->{$searchStr}));
    }

    return $returnStr;
}

sub fixName {
    my ($searchStr) = $_[0];
    $searchStr =~ s/\*/\'/;
    $searchStr =~ s/5M/7/;
    if ($searchStr eq "H2'") {
	$searchStr = "H2'1";
    }
    if ($searchStr =~ /(\d+)H(\S+)/) {
	$searchStr = "H${2}${1}";
    }
    return $searchStr;
}

sub SaveNOEs {
    my ($ATOMS, $BOX, $frameNum, $noes) = @_;
    my ($i, $dim, $dist, $atom1, $atom2);

    for $i (@{ $noes }) {
	$atom1 = $i->{ATOM1};
	$atom2 = $i->{ATOM2};
	$dist = GetBondLength($ATOMS->{$atom1}, $ATOMS->{$atom2});
	push @{ $i->{VALS} }, $dist;
    }

}

1;

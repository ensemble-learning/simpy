package Packages::MESO;

require Exporter;
use strict;
use Cwd;
our (@ISA, @EXPORT, $VERSION, @EXPORT_OK);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(GetMesoParms CreateMesoModel MakeMesoBGF);
$VERSION = "1.00";

sub GetMesoParms {
    my ($parm) = $_[0];
    my ($curr_id, $rec, $mykey, %PARMS_HOLDER, $i, $j);

    $curr_id = -1;
    $rec = ();

    open PARMFILE, $parm or die "Cannot open parameter file $parm: $!\n";
    while(<PARMFILE>) {
	chomp;
	if ($_ =~ /^BEAD_ID: (\d+)/) {
	    if ($curr_id > -1) {
		%{ $PARMS_HOLDER{"BEADS"}{$curr_id} } = %{ $rec };
	    }
	    $curr_id = $1;
	    $rec = ();
        } elsif ($curr_id > -1 and $_ =~ /^BEAD_HBOND_(\w+)_(\d+):\s(\w+)/) {
	    $rec->{"HBONDS"}{$1}{$2} = $3;
	} elsif ($curr_id > -1 and $_ =~ /^BEAD_(\w+): (.+)/) {
	    $rec->{$1} = $2;
	} elsif ($_ =~ /^SAME_(\w+)_\d+: (\d+)\-(\d+)\s(\d+)\-(\d+)/) {
	    $i = sprintf("%03d%03d", $2, $3);
	    $j = sprintf("%03d%03d", $4, $5);
	    $PARMS_HOLDER{"SAME"}{$1}{$i} = $j;
	} elsif ($_ =~ /^MASS_(\w+)\s+(\d+\.*\d*)/) {
	    $PARMS_HOLDER{"MASSES"}{$1} = $2;
	} elsif ($_ =~ /^EQUIV_\d+: (\d+)\s+(\d+)/) {
	    $PARMS_HOLDER{"EQUIVALENCE"}{$1} = $2;
	}
	
    }
    
    close PARMFILE;

    %{ $PARMS_HOLDER{"BEADS"}{$curr_id} } = %{ $rec };

    die "ERROR: Invalid parameter file $parm\n"
	if ($curr_id == -1);
    return (\%PARMS_HOLDER);
}

sub getEquiv {
    my ($parms) = $_[0];
    my ($i, $j, %EQUIV);
    
    for $i (keys %{ $parms->{"EQUIVALENCE"} } ) {
	$j =  $parms->{"EQUIVALENCE"}{$i};
	$EQUIV{ $parms->{"BEADS"}{$i}{"NAME"} }= $parms->{"BEADS"}{$j}{"NAME"};
    }

    return \%EQUIV;
}


sub CreateMesoModel {
    my ($atoms, $parms) = @_;
    my ($EQUIV, $resId, $currId, $strandEnd, @tmp, $oldResId, $atmCounter, $atmLabel);
    my ($atmRes, $beadCounter, $baseRes, %MESO_MODEL, $nonMember, $beadMember, $beadName);
    my ($hbondA, $hbondD, $ff, $counter, $index, $total, $hbType);

    $EQUIV = getEquiv($parms);
    $resId = $currId = -1;
    $strandEnd = 1;
    @tmp = sort numerically keys %{ $atoms };
    $resId = $atoms->{$tmp[0]}{"RESNUM"};
    $oldResId = $resId;

    for $atmCounter (@tmp) {
	$resId = $atoms->{$atmCounter}{"RESNUM"};
	$atmLabel = Trim($atoms->{$atmCounter}{"ATMNAME"});
	$atmRes = Trim($atoms->{$atmCounter}{"RESNAME"});
	$strandEnd = getStrand($atmRes, $strandEnd, $resId, $oldResId);
	$oldResId = $resId;

	$atmRes =~ s/\d$//;
	
	for $beadCounter (keys %{ $parms->{"BEADS"} }) {
            $baseRes = $parms->{"BEADS"}{$beadCounter}{"RES"};
	    next if (! isMember($baseRes,$atmRes));
	    $beadName = $parms->{"BEADS"}{$beadCounter}{"NAME"};
	    if (exists($EQUIV->{$beadName}) and $strandEnd == 1) {
		$beadMember = $parms->{"BEADS"}{$beadCounter}{"MEMBERS"};
		$nonMember = $parms->{"BEADS"}{$beadCounter}{"NONMEMBERS"};
		$hbondA = $parms->{"BEADS"}{$beadCounter}{"HBONDS"}{"ACCEPTOR"};
		$hbondD = $parms->{"BEADS"}{$beadCounter}{"HBONDS"}{"DONOR"};
		$beadName = $EQUIV->{$beadName};
		$beadCounter = $parms->{"EQUIVALENCE"}{$beadCounter};
	    } else {
		$beadMember = $parms->{"BEADS"}{$beadCounter}{"MEMBERS"};
		$nonMember = $parms->{"BEADS"}{$beadCounter}{"NONMEMBERS"};
                $hbondA = $parms->{"BEADS"}{$beadCounter}{"HBONDS"}{"ACCEPTOR"};
                $hbondD = $parms->{"BEADS"}{$beadCounter}{"HBONDS"}{"DONOR"};		
	    }
	    $ff = $parms->{"BEADS"}{$beadCounter}{"ELEMENT"};
	    if (isMember($beadMember,$atmLabel)) {
		$MESO_MODEL{$resId}{$beadCounter}{"ATOMLIST"}{$atmLabel} = 1;
		$MESO_MODEL{$resId}{$beadCounter}{"ATOMS"}{$atmCounter} = 1;
		$MESO_MODEL{$resId}{$beadCounter}{"INDEX"} = $beadCounter;
		updateCOM(\%{ $MESO_MODEL{$resId}{$beadCounter} },\%{ $atoms->{$atmCounter} }, $parms);
		$MESO_MODEL{$resId}{$beadCounter}{"FFTYPE"} = $ff;
		$MESO_MODEL{$resId}{$beadCounter}{"SOLUTE"} = $atoms->{$atmCounter}{"SOLUTE"};
		$MESO_MODEL{$resId}{$beadCounter}{"ATMNAME"} = $beadName;
		$MESO_MODEL{$resId}{$beadCounter}{"RESNAME"} = Trim($atoms->{$atmCounter}{"RESNAME"});
		
		for $counter ("CHARGE", "RADII", "NUMBONDS", "LONEPAIRS", "ELEMENT") {
		    $MESO_MODEL{$resId}{$beadCounter}{$counter} = $parms->{"BEADS"}{$beadCounter}{$counter};
		}
		$atoms->{$atmCounter}{"MESO"}{"RES"} = $resId;
		$atoms->{$atmCounter}{"MESO"}{"BEAD"} = $beadCounter;
	    } elsif (isMember($nonMember,$atmLabel)) {
		$MESO_MODEL{$resId}{$beadCounter}{"ATOMLIST"}{$atmLabel} = 1;
		$atoms->{$atmCounter}{"MESO"}{"RES"} = $resId;
		$atoms->{$atmCounter}{"MESO"}{"BEAD"} = $beadCounter;
	    }

#	    Hydrogen Bonds Acceptors
	    $hbType = checkHB($hbondA,$atmLabel);
	    if ($hbType) {
		$MESO_MODEL{$resId}{$beadCounter}{"HBONDS"}{"ACCEPTOR"}{$hbType}{"ATOMS"} = $atmLabel;
		for $counter ("XCOORD", "YCOORD", "ZCOORD") {
		    $MESO_MODEL{$resId}{$beadCounter}{"HBONDS"}{"ACCEPTOR"}{$hbType}{$counter} = 
			$atoms->{$atmCounter}{$counter};
		}
	    }
	    
	    #Hydrogen Bond Donors
	    $hbType = checkHB($hbondD,$atmLabel);
	    if ($hbType) {
	        $MESO_MODEL{$resId}{$beadCounter}{"HBONDS"}{"DONOR"}{$hbType}{"ATOMS"} = $atmLabel;
		for $counter ("XCOORD", "YCOORD", "ZCOORD") {
		    $MESO_MODEL{$resId}{$beadCounter}{"HBONDS"}{"DONOR"}{$hbType}{$counter} = 
			$atoms->{$atmCounter}{$counter};
		}
	    }
	}
    }
    
    for $resId (keys %MESO_MODEL) {
	for $currId (keys %{ $MESO_MODEL{$resId} }) {
	    $beadCounter = $MESO_MODEL{$resId}{$currId};
            if (! $beadCounter->{"COM"}{"TotalMass"} ) {
                delete $MESO_MODEL{$resId}{$currId};
                next;
            }
	    $beadCounter->{"XCOORD"} = $beadCounter->{"COM"}{"X"}/$beadCounter->{"COM"}{"TotalMass"};
	    $beadCounter->{"YCOORD"} = $beadCounter->{"COM"}{"Y"}/$beadCounter->{"COM"}{"TotalMass"};
	    $beadCounter->{"ZCOORD"} = $beadCounter->{"COM"}{"Z"}/$beadCounter->{"COM"}{"TotalMass"};
	    next
		if (! exists($beadCounter->{"HBONDS"}));

	    for $hbType (keys %{ $beadCounter->{"HBONDS"} }) {
		$index = 50;
		$index = 100 if ($hbType eq "ACCEPTOR");
		for $counter (keys %{ $beadCounter->{"HBONDS"}{$hbType} }) {
		    createHBondBead(\%{ $MESO_MODEL{$resId} }, $index, $hbType, $counter, $currId);
		}
	    }
	}
    }
	
    return \%MESO_MODEL;	
}

sub checkHB {
    my ($hbList, $atmLabel) = @_;
    my ($counter, $returnVal);

    $returnVal = 0;
    for $counter (keys %{ $hbList }) {
	if ($hbList->{$counter} eq $atmLabel) {
	    $returnVal = $counter;
	    last;
	}
    }

    return $returnVal;
}

sub getStrand {
    my ($resName, $strandEnd, $resId, $oldResId) = @_;
    
    $resName =~ /(\d+)/;
    if (defined($1)) {
	if ($1 == 5) {
	    $strandEnd = 1;
	} elsif ($resId > $oldResId) {
	    $strandEnd *= -1;
	}
    } else {
	if ($resId > $oldResId) {
	    $strandEnd *= -1;
	}
    }

    return $strandEnd;
}

sub getBeadInfo {
    my ($beadName, $beadCounter, $strandEnd, $parms) = @_;

    if ($beadName =~ /PHP|PHO/) {
	if ($strandEnd == 1) {
	    $beadName = "PHO";
	    $beadCounter = 1;
	} else {
	    $beadName = "PHP";
	    $beadCounter = 2;
	}
    } elsif ($beadName =~ /SUG|SUS/) {
	if ($strandEnd == 1) {
	    $beadName = "SUG";
	    $beadCounter = 3;
	} else {
	    $beadName = "SUS";
	    $beadCounter = 4;
	}
    }
    
    return ($beadName, $beadCounter);
    
}

sub isMember {
    my ($searchStr, $searchItem) = @_;
    my (@tmp, $returnVal, $counter);
    
    $returnVal = 0;
    if (! defined($searchStr)) {
	$returnVal = 0;
    } else {
	
	if ($searchStr !~ /,/) {
	    push @tmp, $searchStr;
	} else {
	    @tmp = split /,/, $searchStr;
	}
    
	for $counter (@tmp) {
	    if (lc($counter) eq lc($searchItem)) {
		$returnVal =  1;
		last;
	    }
	}
    }

    return $returnVal;
}

sub createHBondBead {
    my ($res, $index, $type, $hbIndex, $parent) = @_;
    my (%Curr_Bead, $counter, $hbond, $atmName, $element, $companion);
    
    $hbond = $res->{$parent}{"HBONDS"}{$type}{$hbIndex};
    $atmName = $hbond->{"ATOMS"};
    $Curr_Bead{"ATOMLIST"}{$atmName} = 1;
    $Curr_Bead{"INDEX"} = $index + $hbIndex;
    
    for $counter ("XCOORD", "YCOORD", "ZCOORD", "ATOMS") {
	$Curr_Bead{$counter} = $hbond->{$counter};
    }
    $Curr_Bead{"CHARGE"} = 0.0;
    $Curr_Bead{"RADII"} = 1.058;
    $Curr_Bead{"NUMBONDS"} = 1;
    $Curr_Bead{"LONEPAIRS"} = 0;
    if ($type eq "ACCEPTOR") {
	$element = "H";
	$companion = "D";
    } else {
	$element = "D";
	$companion = "H";
    }
    $Curr_Bead{"ELEMENT"} = $element;
    $Curr_Bead{"ATMNAME"} = $element . $hbIndex;
    $Curr_Bead{"FFTYPE"} = $element . "_" . $hbIndex;
    $Curr_Bead{"PARENT"} = $parent;
    $Curr_Bead{"SOLUTE"} = 1;
    $Curr_Bead{"COMPANION"} = $companion;

    %{ $res->{($index + $hbIndex)} } = %Curr_Bead;
}

sub updateCOM {
    my ($curr_bead, $atmCounter, $parms) = @_;
    my ($atm_mass);

    $atm_mass = getAtmMass(Trim($atmCounter->{"ATMNAME"}), $parms);

    $curr_bead->{"COM"}{"X"} += $atmCounter->{"XCOORD"} * $atm_mass;
    $curr_bead->{"COM"}{"Y"} += $atmCounter->{"YCOORD"} * $atm_mass;
    $curr_bead->{"COM"}{"Z"} += $atmCounter->{"ZCOORD"} * $atm_mass;
    $curr_bead->{"COM"}{"TotalMass"} += $atm_mass;

#    return $curr_bead;
}

sub MakeMesoBGF {
    my ($Meso, $parms, $atoms, $atomBonds) = @_;
    my ($HEADER, $resC, $beadC, $bead, $index, $BGF, $BONDS, $resName);

    $index = 1;
    for $resC (sort numerically keys %{ $Meso } ) {
	for $beadC (sort numerically keys %{ $Meso->{$resC} } ) {
	    $bead = $Meso->{$resC}{$beadC};
	    $bead->{"INDEX"} = $index;
	    $bead->{"LABEL"} = "ATOM";
	    $bead->{"RESNUM"} = $resC;
	    $bead->{"ID"} = $beadC;
	    $bead->{"ATMNAME"} = Trim($bead->{"ATMNAME"});
	    $bead->{"FFTYPE"} = Trim($bead->{"FFTYPE"});
	    if (exists($bead->{"PARENT"})) {
		$bead->{"PARENTID"} = $Meso->{$resC}{ $bead->{"PARENT"} }{"INDEX"};
	    }
	    if (! exists($bead->{"RESNAME"})) {
		$bead->{"RESNAME"} = $resName;
	    } else {
		$resName = $bead->{"RESNAME"};
	    }
	    %{ $BGF->{$index} } = %{ $bead };
	    $index++;
	}
    }

    $BONDS = buildBondList($BGF, $Meso, $parms, $atoms, $atomBonds);
    return ($BGF, $BONDS);
}

sub Trim {
    my ($inString) = $_[0];

    for ($inString) {
        s/^\s+//;
        s/\s+$//;
    }

    return $inString;
}

sub buildBondList {
    my ($BGF, $Meso, $parms, $atoms, $atomsBonds) = @_;
    my ($mType, $beadC, $resC, $counter, $bType, $bond);
    my (%BONDS, $parent, %CONS);

    for $counter (keys %{ $atoms }) {
	$resC = $atoms->{$counter}{"MESO"}{"RES"};
	next if (! $resC);
	$beadC = $atoms->{$counter}{"MESO"}{"BEAD"};
	$mType = $Meso->{ $resC }{ $beadC }{"INDEX"};
	next if (! defined($mType));
	for $bond (@{ $atomsBonds->{$counter} }) {
	    $resC = $atoms->{$bond}{"MESO"}{"RES"};
	    next if (! $resC);
	    $beadC = $atoms->{$bond}{"MESO"}{"BEAD"};
	    $bType = $Meso->{ $resC }{ $beadC }{"INDEX"};
	    next if (! defined($bType));
	    if ($mType != $bType) {
		$BONDS{$mType}{$bType} = 1;
		$BONDS{$bType}{$mType} = 1;
	    }
	}
    }
    
    #h-bonds
    for $beadC (keys %{ $BGF } ) {
        next if (! $BGF->{$beadC}{"SOLUTE"} or ! exists($BGF->{$beadC}{"PARENT"}));
	$parent = $BGF->{$beadC}{"PARENTID"};
	$BONDS{$beadC}{$parent} = 1;
	$BONDS{$parent}{$beadC} = 1;
    }

    for $bond (keys %{ $BGF }) {
	if (! exists($BONDS{$bond})) {
	    $CONS{$bond} = ();
	} else {
	    @{ $CONS{$bond} } = sort numerically keys %{ $BONDS{$bond} };
	}
    }
    
    return \%CONS;
}

sub getSearchList {
    my ($aType, $bType, $parms) = @_;
    my (@searchList, %keyList, $key1, $key2);

    @{ $keyList{$bType} } = ($aType, "X");
    @{ $keyList{$aType} } = ($bType, "X");
    
    if (exists($parms->{"EQUIVALENCE"}{$bType})) {
	for $key1 (keys %{ $parms->{"EQUIVALENCE"}{$bType} }) {
	    push @{ $keyList{$aType} }, $key1;
	}
    }

    if (exists($parms->{"EQUIVALENCE"}{$aType})) {
	for $key1 (keys %{ $parms->{"EQUIVALENCE"}{$aType} }) {
	    push @{ $keyList{$bType} }, $key1;
	}
    }
    
    for $key1 (@{ $keyList{$bType} }) {
        for $key2 (@{ $keyList{$aType} }) {
            push @searchList, $key1 . "-" . $key2;
        }
    }
    
    return \@searchList;
}

sub getAtmMass {
    my ($atmLabel, $parms) = @_;
    my ($returnval, $in_label);

    $returnval = 0;
    
    if ($atmLabel =~ /^\d*([a-z]+)/i) {
	$in_label = $1;
   
	if ($parms->{"MASSES"}{$in_label}) {
	    $returnval = $parms->{"MASSES"}{$in_label};
	}
    }
    return $returnval;
}

sub numerically {
    ($a<=>$b);
}

1;

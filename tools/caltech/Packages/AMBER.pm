package Packages::AMBER;

require Exporter;
use Packages::General qw(PrintProgress);
use constant q2e => 18.2223;
use Math::Trig qw(pi);

use strict;

our (@EXPORT_OK, @ISA, @EXPORT, $VERSION);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(parseAmberLib getTopInfo updateAtmLabels parseCoordFile getConn CreateAmberTrj ConvertAmberBox
		getOpts loadCnvFile parseAmberFF centerAtoms GetVels ParseAmberTrj AmberLib GetAmberByteOffset);
$VERSION = "1.00";

sub numerically { ($a<=>$b); }

sub parseAmberLib {
    my ($parmFile, $PARMS) = @_;    
    my ($parmName, $parmType, @HEADERS, $parmCounter, @VALS, $i, $inStr);
    my (@ATOMS, $currParm, $offset, $tmp, $tmp1, $rec, $isValid);
    
    $parmCounter = 0;
    $parmName = "";
    open PARM, $parmFile or die "ERROR: Cannot open file $parmFile: $!\n";
    while (<PARM>) {
	chomp;
	$inStr = '';
	$inStr = $_;
	$inStr =~ s/\!|\"|\'//g;
	if ($inStr =~ /entry\.(\w+)\.\w+\.(\w+)\s+(\w+)\s+(.+)/) {
            if ($parmName ne "") {
                #print "$parmName $parmType: $parmCounter\n";
            }
	    $parmName = $1;	
	    $parmType = $2;
	    @HEADERS = split /\s+/, $4;
	    if ($3 =~ /table|array|str|dbl/) {
		$parmCounter = 1;
	    } else {
		$parmCounter = 0;
	    }
        } elsif ($parmType =~ /Orders/ and $inStr =~ /^\s\"(\S+)\"/ and $parmCounter) {
            $PARMS->{$parmType}{$parmCounter} = $1;
            $parmCounter++;
	} elsif ($inStr =~ /^\s\"(\S+)\"/ and $parmCounter) {
            @ATOMS = @VALS = ();
            $tmp = $1;
            if ($tmp eq "?") {
                $tmp = "X";
            }
            $currParm = \%{ $PARMS->{$parmType}{$tmp} };
            $tmp = $tmp1 = $';
            while ($tmp =~ /\s+\"(\S+)\"/g) {
                if ($1 eq "?") {
                    push @ATOMS, "X";
                } else {
                    push @ATOMS, $1;
                }
                $tmp1 = $';
            }
            while ($tmp1 =~ /\s+(\S+)/g) {
                push @VALS, $1;
            }
            for (0 .. ($#ATOMS - 1)) {
                $currParm = \%{ $currParm->{$ATOMS[$_]} };
            }
            
            $isValid = 0;
            if ($#ATOMS == -1) {
                $isValid = 1;
            } elsif (! exists($currParm->{ $ATOMS[$#ATOMS] })) {
                $isValid = 1;
            } elsif ($parmName eq $currParm->{ $ATOMS[$#ATOMS] }{"parmName"}) {
                $isValid  = 1;
            } elsif ($parmName eq "frcmod03") {
                $isValid = 1;
            }
            if ($isValid) {
                if ($#ATOMS > -1) {
                    $currParm = \%{ $currParm->{ $ATOMS[$#ATOMS] } };
                }
                $currParm->{"counter"} = $parmCounter;
                $currParm->{"parmName"} = $parmName;
                $offset = ($#ATOMS + 2) * 2 + 1;
                $rec = ();
 		for $i (0 .. $#VALS) {
		    last
		        if ($#HEADERS < (($i * 2) + $offset));
		    $VALS[$i] =~ s/\"//g;
		    $rec->{ $HEADERS[($i * 2) + $offset] } = $VALS[$i];
		}
               
                if ($rec) {
                    push @{ $currParm->{"VALS"} }, $rec;
                }
	    }
 	    $parmCounter++;
	}
    }
    #print "$parmName $parmType: $parmCounter\n";
    close PARM;
    
    die "ERROR: $parmFile do not contain any valid information\n"
        if (! $PARMS);
    return $PARMS;
}

sub getTopInfo {
    my ($tFile, $OPTS) = @_;
    my ($DATA, @tmp);

    $DATA = parseTopFile($tFile, $OPTS); # Get the topology information
    addResInfo(\%{ $DATA->{"ATOMS"} }, \%{ $DATA->{"RESIDUES"} }); # add the residue information
    # Now update the valence information
    $DATA->{"BONDLIST"} = getValenceList(\%{ $DATA->{"BONDLIST" } }, 3); #Bonds
    $DATA->{"ANGLELIST"} = getValenceList(\%{ $DATA->{"ANGLELIST" } }, 4); # Angles
    $DATA->{"DIHEDRALLIST"}  = getValenceList(\%{ $DATA->{"DIHEDRALLIST"} }, 5); #Torsions and inversions

    return ($DATA, $DATA->{"COUNT"}{"ATOMS"});
}

sub updateAtmLabels {
    my ($ATOMS, $cnv_file, $RES) = @_;
    my ($atom, $resName, $atmName);

    for $atom (keys %{ $ATOMS }) {
	$resName = $ATOMS->{$atom}{"RESNAME"};
	if ($resName =~ /^D(\w)3?5?/) {
	    $resName = $1;
	}
	$atmName = $ATOMS->{$atom}{"ATMNAME"};
        if ($atmName =~ /^(\w+)\'(\d+)/) {
	    $atmName = $2 . $1 . "*";
        }
	$atmName =~ s/\'/\*/g;
	if ($RES) {
	    if (exists($RES->{$resName}{$atmName})) {
		$ATOMS->{$atom}{"FFTYPE"} = $RES->{$resName}{$atmName}{"BGF_LABEL"};
		$ATOMS->{$atom}{"LONEPAIRS"} = $RES->{$resName}{$atmName}{"LONEPAIRS"};
	    } elsif (exists($RES->{"ALL"}{$atmName})) {
		$ATOMS->{$atom}{"FFTYPE"} = $RES->{"ALL"}{$atmName}{"BGF_LABEL"};
		$ATOMS->{$atom}{"LONEPAIRS"} = $RES->{"ALL"}{$atmName}{"LONEPAIRS"};
	    } else {
		if (length($resName) > 3) {
		    $ATOMS->{$atom}{"FFTYPE"} = substr($resName, 0, 3);
		} else {
		    $ATOMS->{$atom}{"FFTYPE"} = $resName;
		}
		$ATOMS->{$atom}{"LONEPAIRS"} = 0;
	    }
	}
    }
}

sub getConn {
    my ($BONDLIST, $ATOMS) = @_;
    my ($atom1, $atom2, $bondC, %CONN);
    
    for $bondC (keys %{ $BONDLIST } ) {
	($atom1, $atom2) = @{ $BONDLIST->{$bondC}{"ATOMS"} };
	push @{ $CONN{$atom1} }, $atom2;
	push @{ $CONN{$atom2} }, $atom1;
	$ATOMS->{$atom1}{"NUMBONDS"} = $#{ $CONN{$atom1} } + 1;
	$ATOMS->{$atom2}{"NUMBONDS"} = $#{ $CONN{$atom2} } + 1;
    }

    for $atom1 (keys %{ $ATOMS }) {
	if (! exists($CONN{$atom1})) {
	    $ATOMS->{$atom1}{"NUMBONDS"} = 0;
	    @{ $CONN{$atom1} } = ();
	}
    }
	
    return \%CONN;
}

sub loadCnvFile {
    my ($cnv_file) = $_[0];
    my ($in_data, $is_valid, $rec, $label, %DATA);
    $is_valid = 0;

    open INFILE, $cnv_file || die "Cannot load conversion file $cnv_file: $!\n";
    while (<INFILE>) {
	chomp;
	if ($_ =~ /^\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\-?\d+\.\d+)/) {
	    $is_valid = 1;
	    $label = $2;

	    $rec = (
		    {
			"BGF_LABEL" => $3,
			"BONDS"     => $4,
			"LONEPAIRS" => $5,
			"CHARGE"    => $6,
		    }
		    );
	    if ($1 =~ /\*\*\*/) {
		$DATA{"ALL"}{$label} = $rec;
	    } else {
		$DATA{$1}{$label} = $rec;
	    }
	    
	}
    }
    close INFILE;

    die "Unable too parse conversion file $cnv_file\n"
	if (! $is_valid);

    return \%DATA;
}

sub parseTopFile {
    my ($tFile, $OPTS) = @_;
    my ($start, %DATA, $parmC, $increment, $fieldSize);
    my ($parmType, $typeLabel);

    $start = $increment = $parmC = $fieldSize = 0;
 
    open TOPFILE, $tFile or die "ERROR: Cannot open topology file $tFile: $!\n";
    while (<TOPFILE>) {
	chomp;
	if ($_ =~ /^\%FLAG (\S+)/) {
	    if (exists($OPTS->{$1})) {
		$start = $parmC = 1;
		$increment = $fieldSize = 0;
		$parmType = $OPTS->{$1}{"TYPE"};
		if ($OPTS->{$1}{"HEADER"}) {
		    $typeLabel = $OPTS->{$1}{"HEADER"};
		} else {
		    $typeLabel = "DATA";
		}
		#print "FOUND $parmType: $typeLabel\n";
	    } else {
		$start = $increment = $parmC = $fieldSize = 0;
	    }
	} elsif ($_ =~ /^\%FORMAT\((\d+)\w(\d+)\.*\d*\)/ and $start) {
	    $increment = $1;
	    $fieldSize = $2;
	} elsif ($start and $increment) {
	    while ($_ =~ /(.{$fieldSize})/g) {
		$DATA{$parmType}{$parmC}{$typeLabel} = $1;
		$DATA{$parmType}{$parmC}{$typeLabel} =~ s/^\s+//;
		$DATA{$parmType}{$parmC}{$typeLabel} =~ s/\s+$//;
		$DATA{"COUNT"}{$parmType} = $parmC;
		$parmC++;
	    }
	}
    }
    close TOPFILE;

    die "ERROR: No valid data found while parsing $tFile\n" if (! %DATA);

    return \%DATA;
}

sub parseCoordFile {
    my ($data, $cFile, $totAtms) = @_;
    my ($patern, $atomC, $timeStamp, $datType); 
    my ($hasVel, $atoms, $box, $lastline);

    $atoms = \%{ $data->{"ATOMS"} };
    $box = \%{ $data->{"BOX"} };
    $atomC = $timeStamp = $hasVel = 0;
    $datType = "COORD";
    $patern = '\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)';
    open COORDFILE, $cFile or die "ERROR: Cannot open coordinate file $cFile: $!\n";
    while (<COORDFILE>) {
	chomp;
	if ($_ =~ /^\s*(\d+)\s*(\S*)$/) {
	    if ($2) {
		$timeStamp = $1;
	    }
	    $totAtms = $1;
	    $atomC = 1;
	} elsif ($atomC > 0 and ! $hasVel) {
	    while ($_ =~ /$patern/g) {
		$atoms->{$atomC}{"X" . $datType} = $1;
		$atoms->{$atomC}{"Y" . $datType} = $2;
		$atoms->{$atomC}{"Z" . $datType} = $3;
		if ($datType eq "VEL") {
		    $atoms->{$atomC}{XVEL} /= (1000/20.455); #amber converstion factor
                    $atoms->{$atomC}{YVEL} /= (1000/20.455);
                    $atoms->{$atomC}{ZVEL} /= (1000/20.455);
		}
		$atomC++;
		if ($atomC > $totAtms and $datType eq "COORD") {
		    $datType = "VEL";
		    $atomC = 1;
		} elsif ($atomC > $totAtms and $datType eq "VEL") {
		    $hasVel = 1;
		    last;
		}
	    }
	} 
	$lastline = $_;
    }
    
    close COORDFILE;
    
    if ($atomC > $totAtms) {
	if ($lastline =~ /$patern/) {
	    $box->{2}{"DATA"} = $1;
	    $box->{3}{"DATA"} = $2;
	    $box->{4}{"DATA"} = $3;
	}
    }

    die "ERROR: No valid data found while parsing $cFile\n"
	if (! $totAtms);

    return ($timeStamp, $hasVel);
}

sub addResInfo {
    my ($ATOMS, $RES) = @_;
    my ($resC, $startAtm, $endAtm, $atomC);
    my (@tmp, $currRes, $resName);

    $startAtm = $endAtm = 0;
    @tmp = sort numerically keys %{ $RES };
    for $resC (@tmp) {
	$currRes = $resC - 1;
	if ($startAtm == 0) {
	    $startAtm = $RES->{$resC}{"STARTATM"};
	    $resName = $RES->{$resC}{"NAME"};
	} else {
	    $endAtm =  $RES->{$resC}{"STARTATM"} - 1;
	    for $atomC ($startAtm .. $endAtm) {
		$ATOMS->{$atomC}{"RESNAME"} = $resName;
		$ATOMS->{$atomC}{"RESNUM"} = $currRes;
		$ATOMS->{$atomC}{"CHARGE"} = $ATOMS->{$atomC}{"CHARGE"}/q2e;
		$ATOMS->{$atomC}{"LONEPAIRS"} = 0;
	    }
	    $startAtm = $RES->{$resC}{"STARTATM"};
	    $resName = $RES->{$resC}{"NAME"};
	}
    }

    # get put the last residue
    $currRes++;
    @tmp = sort numerically keys %{ $ATOMS };
    $endAtm = $tmp[$#tmp];
    for $atomC ($startAtm .. $endAtm) {
	$ATOMS->{$atomC}{"RESNAME"} = $resName;
	$ATOMS->{$atomC}{"RESNUM"} = $currRes;
	$ATOMS->{$atomC}{"CHARGE"} = $ATOMS->{$atomC}{"CHARGE"}/q2e;
	$ATOMS->{$atomC}{"LONEPAIRS"} = 0;
    }
}

sub getValenceList {
    my ($VALENCE, $totTerms) = @_;
    my ($valStr, $i, $typeName, @PARMS, %tmp, $j, $k, $typeID, $atmIndex);

    $j = 0;
    for $typeName ("wHYDROGEN", "noHYDROGEN") {
	$i = 1;
	while (exists($VALENCE->{$i}{$typeName}) and exists($VALENCE->{($i + $totTerms -1)}{$typeName})) {
	    $typeID = $VALENCE->{$i + $totTerms - 1}{$typeName};
	    if ($totTerms != 5 or ($totTerms == 5 and $typeID > 0)) { #not multidihedral
		$j++;
		for $k (0 .. ($totTerms - 2)) {
		    $atmIndex = ((abs($VALENCE->{$k + $i}{$typeName})/3) + 1);
		    push @{ $tmp{$j}{"ATOMS"} }, $atmIndex;
		}
		if ($totTerms == 5) {
		    if ($atmIndex < 0) { # check for inversion
			$tmp{$j}{"INVERSION"} = 1;
		    } else {
			$tmp{$j}{"INVERSION"} = 0;
		    }
		}
	    }
	    push @{ $tmp{$j}{"TYPE"} }, $typeID;
	    $i += $totTerms;
	}
    }    
	    
    return \%tmp;   
}

sub getOpts {
    my ($inFile) = "/ul/tpascal/scripts/dat/amberTopFile.perldata";
    my (%OPTS);
    scalar eval `cat $inFile` or die "Cannot recreate data in file $inFile: $! $@\n";
    return \%OPTS;
}	      

sub parseAmberFF {
    my ($amberFile, $PARMS) = @_;
    my ($rec, $torsionCounter, $inversionCounter, $section, $i, @tmp, $atom1);

    $PARMS->{"atoms"}{"X"}{"mass"} = 0.00;
    $PARMS->{"atoms"}{"X"}{"polarizability"} = 0.00;
    $torsionCounter = $inversionCounter = 1;
    $section = 0;

    open INFILE, $amberFile or die "Cannot open $amberFile: $!\n";
    while (<INFILE>) {
	chomp;
        if ($_ =~ /^(\w+\*?\-?\+?)\s+\w+\*?\-?\+?/ && $_ !~ /\d+\./ && $section > 4) { #same vdw
	    @tmp = split /\s+/, $_;
	    $atom1 = shift @tmp;
	    next if (! exists($PARMS->{atoms}{$atom1}));
	    for $i (@tmp) {
		$PARMS->{atoms}{$i} = \%{ $PARMS->{atoms}{$atom1} };
	    }
	} elsif ($_ =~ /^\s*(\w+\*?\-?\+?)\s+(\d+\.\d*)\s*(\d*\.*\d*)/) { # atomtype
	    $section = 1;
	    if (! exists($PARMS->{atoms}{$1})) {
		$PARMS->{"atoms"}{$1}{VALS}[0]{mass} = $2;
	    } else {
		$PARMS->{atoms}{$1}{VALS}[0]{r} = $2;
		$PARMS->{atoms}{$1}{VALS}[0]{e} = $3;
	    }
	} elsif ($_ =~ /^(\w+\*?\-?\+?)\s*\-\s*(\w+\*?\-?\+?)\s+(\d+\.\d*)\s+(\d+\.\d*)/) { #bond
	    #if (exists($PARMS->{"atoms"}{$1}) && exists($PARMS->{"atoms"}{$2})) {
		$section = 2;
		$PARMS->{"bonds"}{$1}{$2}{VALS}[0]{kb} = $3;
		$PARMS->{"bonds"}{$1}{$2}{VALS}[0]{r0} = $4;
	    #}
	} elsif ($_ =~ /^(\w+\*?\-?\+?)\s*\-\s*(\w+\*?\-?\+?)\s*\-\s*(\w+\*?\-?\+?)\s+(\d+\.\d*)\s+(\d+\.\d*)/) { #angle
	    #if (exists($PARMS->{"atoms"}{$1}) && exists($PARMS->{"atoms"}{$2}) && exists($PARMS->{"atoms"}{$3})) {
		$section = 3;
		$PARMS->{"angles"}{$1}{$2}{$3}{VALS}[0]{kt} = $4;
		$PARMS->{"angles"}{$1}{$2}{$3}{VALS}[0]{t0} = $5 * pi/180;
	    #}
	} elsif ($_ =~ /^(\w+\*?\-?\+?)\s*\-\s*(\w+\*?\-?\+?)\s*\-\s*(\w+\*?\-?\+?)\s*\-\s*(\w+\*?\-?\+?)\s+(\d+)\s+(\d+\.\d*)\s+(\d+\.\d*)\s+(\-?\d+\.\d*)/) { # Torsion
	    #if (exists($PARMS->{"atoms"}{$1}) and exists($PARMS->{"atoms"}{$2}) and
		#exists($PARMS->{"atoms"}{$3}) and exists($PARMS->{"atoms"}{$4})) {
		$section = 4;    
		$rec = (
			{
			    "scale"       => $5,
			    "kp"          => $6,
			    "p0"          => $7 * pi/180,
			    "n"           => abs($8),
			}
			);
		push @{ $PARMS->{"torsions"}{$1}{$2}{$3}{$4}{"VALS"} }, $rec;
		if (! exists($PARMS->{"torsions"}{$1}{$2}{$3}{$4}{counter})) {
		    $PARMS->{"torsions"}{$1}{$2}{$3}{$4}{counter} = $torsionCounter;
		    $PARMS->{"torsionOrders"}{$torsionCounter} = "0 1 2 3";
		    $torsionCounter++;
		}
	    #}
	} elsif ($_ =~ /^(\w+\*?\-?\+?)\s*\-\s*(\w+\*?\-?\+?)\s*\-\s*(\w+\*?\-?\+?)\s*\-\s*(\w+\*?\-?\+?)\s+(\d+\.\d*)\s+(\d+\.\d*)\s+(\d+\.\d*)/) { # Inversion
	    #if (exists($PARMS->{"atoms"}{$1}) and exists($PARMS->{"atoms"}{$2}) and
		#exists($PARMS->{"atoms"}{$3}) and exists($PARMS->{"atoms"}{$4})) {
		$section = 5;
		$rec = (
			{
			    "kp"          => $5,
			    "p0"          => $6,
			    "n"           => $7,
			    "type"        => "IT_JIKL",
			}
			);
		$PARMS->{"inversions"}{$1}{$2}{$3}{$4}{VALS}[0] = $rec;
		$PARMS->{"inversions"}{$1}{$2}{$3}{$4}{counter} = $inversionCounter;
		$PARMS->{"inversionOrders"}{$inversionCounter} = "2 1 0 3";
	    #}
	}
    }
    close INFILE;

    die "Error while reading amber ff file $amberFile. No valid data found\n"
	if (! keys %{ $PARMS });
    
    return $PARMS;
}

sub centerAtoms {
    # image atoms to the origin based on the center of mass
    my ($ATOMS, $BBOX, $resStart, $resEnd) = @_;
    my ($atom, $currRes, $totalMass,%COM, $dim, %BOXTRANS);
    my (%BOX) = (
		 "X" => {
		     "hi" => $BBOX->{2}{"DATA"},
		     "lo" => 0,
		 },
		 "Y" => {
		     "hi" => $BBOX->{3}{"DATA"},
		     "lo" => 0,
		 },
		 "Z" => {
		     "hi" => $BBOX->{4}{"DATA"},
		     "lo" => 0,
		 },
		 );
    # calculate com of all atoms in residues
    %COM = ();
    $totalMass = 0;
    for $atom (sort numerically keys %{ $ATOMS }) {
	$currRes =  $ATOMS->{$atom}{"RESNUM"};
	if ($currRes <= $resEnd and $currRes >= $resStart) {
	    for $dim ("X", "Y", "Z") {
		$COM{$dim} += $ATOMS->{$atom}{$dim . "COORD"} * $ATOMS->{$atom}{"MASS"};
	    }
	    $totalMass += $ATOMS->{$atom}{"MASS"};
	}
    }
    if (! $totalMass) {
	return ();
    }

    for $dim ("X", "Y", "Z") {
	$COM{$dim} /= $totalMass;
    }

    # determine how far the coordinates are out of the box
    
}

sub GetVels {
    my ($ATOMS) = $_[0];
    my ($atomC, $returnStr, $vel);

  MAINLOOP: for $atomC (sort numerically keys %{ $ATOMS }) {
      $returnStr .= sprintf("%7d", $atomC);
      for $vel ("XVEL", "YVEL", "ZVEL") {
          if (! exists($ATOMS->{$atomC}{$vel})) {
              $returnStr = "";
              last MAINLOOP;
          } else {
              $returnStr .= sprintf("%12.6f", $ATOMS->{$atomC}{$vel});
          }
      }
      $returnStr .= "\n";
  }

    if ($returnStr) {
        $returnStr = "Velocities\n\n$returnStr";
    }
    return $returnStr;
}

sub GetAmberByteOffset {
    my ($trjSelection, $amberTrj, $totAtms) = @_;
    my ($isNull, $readStart, $atomC, $readEnd, $i, $valid, $tot, $snapByteLen, $count);

    print "Getting byte offset for AMBER trajectory $amberTrj...";
    $readStart = $readEnd = -1; # byte when starting reading 
    $atomC = 0;
    $valid = $i = 0;
    $isNull = 0;
    $isNull = 1 if (! keys %{ $trjSelection });
    open AMBERTRJ, $amberTrj or die "ERROR: Cannopt open AMBER trajectory file $amberTrj: $!\n";
    while (<AMBERTRJ>) {
	chomp;
	if ($_ =~ /^\s*(\-?\d+\.\d{3})/ && $readStart == -1) { # first time through
	    $readStart = tell(AMBERTRJ) - 81; #we've started reading so record the current byte
	}
	if ($readStart > -1 && $atomC < $totAtms) {
            while ($_ =~ /(\-?\d+\.\d{3})/g) {
                $i++;
                if ($i == 3) {
		    $i = 0;
                    $atomC++;
                }
	    }
	} elsif ($atomC == $totAtms && $_ =~ /^\s*(\d+\.\d{3})\s+(\d+\.\d{3})\s+(\d+\.\d{3})\s*$/) { #box info
	    $atomC++;
	    next;
	} elsif ($atomC	>= $totAtms) {
	    $readEnd = tell(AMBERTRJ) - 81; #this is the enda
	    last;
	}
    }
    close AMBERTRJ;
    
    $snapByteLen = $readEnd - $readStart;
    $tot = (-s $amberTrj);
    $tot -= $readStart;
    $tot /= $snapByteLen;
    
    for $i (keys %{ $trjSelection }) {
	$trjSelection->{$i} = -1;
    }
    
    for $i (1 .. $tot) {
	if ($isNull  or exists($trjSelection->{$i})) { # record byte count
	    $trjSelection->{$i} = $snapByteLen * ($i - 1) + $readStart;
	    $valid++;
	}
    }
    
    $count = 0;
    for $i (keys %{ $trjSelection }) {
	if ($trjSelection->{$i} == -1) {
	    delete $trjSelection->{$i};
	} else {
	    $count++;
	}
    }
    
    die "ERROR: No valid trajectory selection or file is not valid AMBER trajectory\n" if (! $valid);
    
    print "using $count of $tot snapshots...Done\n";
}

sub ParseAmberTrj {
    my ($ATOMS, $amberTrj, $trjSelect, $totAtms, $analFunc, $printStr, $fileOut) = @_;
    my (@dim, $ATOMSINFO, $atomC, $i, @tmp, $j);
    my ($BOX, $start, $currPos, $strLen, $tot, $count, $inStr);

    print "${printStr}Calculating time remaining\r" if (defined($printStr));
    $tot = scalar keys %{ $trjSelect };

    $start = time();
    $currPos = 0;
    @dim = ("XCOORD","YCOORD","ZCOORD");
    @tmp = keys %{ $ATOMS->{1} };
    @tmp = grep !/COORD/, @tmp;

    open AMBERTRJ, $amberTrj || die "ERROR: Cannopt open AMBER trajectory file $amberTrj: $!\n";
    for $count (sort numerically keys %{ $trjSelect }) {
	next if (! seek(AMBERTRJ, $trjSelect->{$count},0) );
	$currPos++;
	$i = 0;
	$atomC = 1;
	$BOX = ();
	while (<AMBERTRJ>) {
	    chomp;
	    $inStr = $_;
	    if ($atomC <= $totAtms) {
		while ($inStr =~ /(\-?\d+\.\d{3})/g) {
		    $ATOMSINFO->{$atomC}{ $dim[$i] } = $1;
		    $i++;
		    if ($i == 3) {
			$i = 0;
			for $j (@tmp) {
			    $ATOMSINFO->{$atomC}{$j} = $ATOMS->{$atomC}{$j};
			}
			$atomC++;
		    }
		}
	    } elsif ($atomC == ($totAtms + 1) && $inStr =~ /^\s*(\d+\.\d{3})\s+(\d+\.\d{3})\s+(\d+\.\d{3})\s*$/) { #box info
		$BOX->{2}{DATA} = $1;
		$BOX->{3}{DATA} = $2;
		$BOX->{4}{DATA} = $3;
	        $BOX->{1}{DATA} = 90.00;
		$atomC++;
		next;
	    } elsif ($atomC > $totAtms) {
		$analFunc->($ATOMSINFO, $BOX, $count, $fileOut);
		$strLen = PrintProgress($currPos, $tot, $start, $printStr);
		last;
	    }
	}
    }
    close AMBERTRJ;
    if ($BOX) {
	$count = scalar(keys %{ $trjSelect });
	$analFunc->($ATOMSINFO, $BOX, $count, $fileOut);
	$strLen = PrintProgress($currPos, $tot, $start, $printStr);
    }
    printf "$printStr%-${strLen}s\n", "Done" if (defined($printStr));
}

sub ParseAmberTrjOld {
    my ($ATOMS, $amberTrj, $trjSelect, $totAtms, $analFunc, $printStr, $fileOut) = @_;
    my ($filesize, @dim, $tot, $endFrame, $ATOMSINFO, $atomC, $i, $offSet, @tmp);
    my ($currFrame, $validFrame, $BOX, $start, $currPos, $index, $strLen, $j);

    print "${printStr}Calculating time remaining\r";
    $filesize = -s $amberTrj;
    $tot = $endFrame = $index = $i = 0;
    @tmp = keys %{ $ATOMS->{1} };

    my ($patern) = '\s+(-?\d+\.\d{3})\s+(-?\d+\.\d{3})\s+(-?\d+\.\d{3})';

    if (keys %{ $trjSelect }) {
        @dim = sort numerically keys %{ $trjSelect };
        $tot = $#dim + 1;
        $endFrame = pop @dim;
    }
    @dim = ("XCOORD","YCOORD","ZCOORD");
    $atomC = $currFrame = 1;
    $start = time();

    open AMBERTRJ, $amberTrj || die "ERROR: Cannopt open AMBER trajectory file $amberTrj: $!\n";
  MAINLOOP: while (<AMBERTRJ>) {
      chomp;
      while ($_ =~ /(\-?\d+\.\d{3})/g) {
	  if ($atomC <= $totAtms) {
	      $ATOMSINFO->{$atomC}{ $dim[$i] } = $1;
	      $i++;
	      if ($i == 3) {
		  $i = 0;
		  for $j (@tmp) {
		      $ATOMSINFO->{$atomC}{$j} = $ATOMS->{$atomC}{$j};
		  }
		  $atomC++;
	      }
	  } else {
	      $offSet = $atomC - $totAtms - 1;
	      if ($offSet < 2) { # store box information
		  $BOX->{ $dim[$offSet] } = $1;
		  $atomC++;
	      } else { # finished reading frame
		  $BOX->{ $dim[2] } = $1;
		  if (! keys %{ $trjSelect } or exists($trjSelect->{$currFrame})) { # this was a valid frame
		      $index++;
		      $analFunc->($ATOMSINFO, $BOX, $currFrame, $fileOut);
		      if (! $tot) {
			  $currPos = tell(AMBERTRJ);
			  $strLen = PrintProgress($currPos, $filesize, $start, $printStr);
		      } else {
			  $strLen = PrintProgress($index, $tot, $start, $printStr);
		      }
		  }
		  $atomC = 1;
		  $i = 0;
		  $currFrame++;
		  last MAINLOOP if ($endFrame and $currFrame > $endFrame);
	      }
	  }
      }
  }
    close AMBERTRJ;
    printf "$printStr%-${strLen}s\n", "Done";
}

sub CreateAmberTrj {
    my ($ATOMS, $BBOX, $OUTFILE) = @_;
    my ($atomC, $counter, @tmp, $dim, $BOX);

    $counter = 0;
    @tmp = sort numerically keys %{ $ATOMS };
    if (defined($BBOX)) {
	for $dim ("X", "Y", "Z") {
	    $BOX->{$dim} = $BBOX->{$dim}{"hi"} - $BBOX->{$dim}{"lo"} if (exists($BBOX->{$dim}));
	    $BOX->{$dim} = $BBOX->{"${dim}COORD"}{"hi"} - $BBOX->{"${dim}COORD"}{"lo"} if (exists($BBOX->{"${dim}COORD"}));
	}
    }

    for $atomC (@tmp) {
	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	    $counter++;
	    printf $OUTFILE "%8.3f", $ATOMS->{$atomC}{$dim};
	    if ($counter == 10) {
		print $OUTFILE "\n";
		$counter = 0;
	    }
	}
    }
    if (defined($BBOX)) {
	for $dim ("X", "Y", "Z") {
	    $counter++;
	    printf $OUTFILE "%8.3f", $BOX->{$dim};
            if ($counter == 10) {
                print $OUTFILE "\n";
                $counter = 0;
            }
	}
	#printf $OUTFILE "%8.3f%8.3f%8.3f", 90,90,90;
	print $OUTFILE "\n";
    } else {
	print $OUTFILE "\n" if ($counter > 0);
    }
}

sub AmberLib {
    my (%LIBS) = (
		      "ACE"=> 1,
		      "ALA"=> 1,
		      "ARG"=> 1,
		      "ASH"=> 1,
		      "ASN"=> 1,
		      "ASP"=> 1,
		      "CALA"=> 1,
		      "CARG"=> 1,
		      "CASN"=> 1,
		      "CASP"=> 1,
		      "CCYS"=> 1,
		      "CCYX"=> 1,
		      "CGLN"=> 1,
		      "CGLU"=> 1,
		      "CGLY"=> 1,
		      "CHCL3BOX"=> 1,
		      "CHID"=> 1,
		      "CHIE"=> 1,
		      "CHIP"=> 1,
		      "CHIS"=> 1,
		      "CILE"=> 1,
		      "CIO"=> 1,
		      "CLEU"=> 1,
		      "CLYS"=> 1,
		      "CMET"=> 1,
		      "CPHE"=> 1,
		      "CPRO"=> 1,
		      "CSER"=> 1,
		      "CTHR"=> 1,
		      "CTRP"=> 1,
		      "CTYR"=> 1,
		      "CVAL"=> 1,
		      "CYM"=> 1,
		      "CYS"=> 1,
		      "CYX"=> 1,
		      "Cl-"=> 1,
		      "Cs+"=> 1,
		      "DA"=> 1,
		      "DA3"=> 1,
		      "DA5"=> 1,
		      "DAN"=> 1,
		      "DC"=> 1,
		      "DC3"=> 1,
		      "DC4"=> 1,
		      "DC5"=> 1,
		      "DCN"=> 1,
		      "DG"=> 1,
		      "DG3"=> 1,
		      "DG5"=> 1,
		      "DGN"=> 1,
		      "DT"=> 1,
		      "DT3"=> 1,
		      "DT5"=> 1,
		      "DTN"=> 1,
		      "GLH"=> 1,
		      "GLN"=> 1,
		      "GLU"=> 1,
		      "GLY"=> 1,
		      "HID"=> 1,
		      "HIE"=> 1,
		      "HIP"=> 1,
		      "HIS"=> 1,
		      "HOH"=> 1,
		      "IB"=> 1,
		      "ILE"=> 1,
		      "K+"=> 1,
		      "LEU"=> 1,
		      "LYN"=> 1,
		      "LYS"=> 1,
		      "Li+"=> 1,
		      "MEOHBOX"=> 1,
		      "MET"=> 1,
		      "MG2"=> 1,
		      "NALA"=> 1,
		      "NARG"=> 1,
		      "NASN"=> 1,
		      "NASP"=> 1,
		      "NCYS"=> 1,
		      "NCYX"=> 1,
		      "NGLN"=> 1,
		      "NGLU"=> 1,
		      "NGLY"=> 1,
		      "NHE"=> 1,
		      "NHID"=> 1,
		      "NHIE"=> 1,
		      "NHIP"=> 1,
		      "NHIS"=> 1,
		      "NILE"=> 1,
		      "NLEU"=> 1,
		      "NLYS"=> 1,
		      "NMABOX"=> 1,
		      "NME"=> 1,
		      "NMET"=> 1,
		      "NPHE"=> 1,
		      "NPRO"=> 1,
		      "NSER"=> 1,
		      "NTHR"=> 1,
		      "NTRP"=> 1,
		      "NTYR"=> 1,
		      "NVAL"=> 1,
		      "Na+"=> 1,
		      "PHE"=> 1,
		      "PL3"=> 1,
		      "POL3BOX"=> 1,
		      "PRO"=> 1,
		      "RA"=> 1,
		      "RA3"=> 1,
		      "RA5"=> 1,
		      "RAN"=> 1,
		      "RC"=> 1,
		      "RC3"=> 1,
		      "RC5"=> 1,
		      "RCN"=> 1,
		      "RG"=> 1,
		      "RG3"=> 1,
		      "RG5"=> 1,
		      "RGN"=> 1,
		      "RU"=> 1,
		      "RU3"=> 1,
		      "RU5"=> 1,
		      "RUN"=> 1,
		      "Rb+"=> 1,
		      "SER"=> 1,
		      "SPC"=> 1,
		      "SPCBOX"=> 1,
		      "T4E"=> 1,
		      "THR"=> 1,
		      "TIP3PBOX"=> 1,
		      "TIP4PBOX"=> 1,
		      "TP4"=> 1,
		      "TP5"=> 1,
		      "WAT"=> 1,
		      "TRP"=> 1,
		      "TYR"=> 1,
		      "VAL"=> 1,
		  );
    
    return \%LIBS;
}

sub ConvertAmberBox {
    my ($box) = $_[0];
    my (%BOX, @dim, $i);

    @dim = ("XCOORD", "YCOORD", "ZCOORD");
    for $i (0 .. 2) {
        $BOX{$dim[$i]}{lo} = 0;
        $BOX{$dim[$i]}{hi} = $box->{$i + 2}{DATA};
        $BOX{$dim[$i]}{len} = $box->{$i + 2}{DATA};
	$BOX{$dim[$i]}{CENTER} = $box->{$i + 2}{DATA}/2;
    }
                                                                                                                
    return \%BOX;


}
1;	      

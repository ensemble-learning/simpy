package Packages::LAMMPS;

require Exporter;
use Packages::General qw(PrintProgress);
use File::Basename qw(basename);
our (@ISA, @EXPORT, $VERSION);
use strict;

@ISA = qw(Exporter);
@EXPORT = qw(ParseLammpsTrajectoryFile ReadDataFile ParseLAMMPSLogFile CreateLAMMPSTrj 
ParseLAMMPSTrj CreateInputFile GetLammpsByteOffset GetLammpsTrjType ConvertLammpsBox);
$VERSION = "1.00";

sub determineIfScaled {
    my ($DATA, $OPTS, $tmp2) = @_;
    my ($i, $j);
    my (@dim) = ("XCOORD", "YCOORD", "ZCOORD");

    $OPTS->{scaled} = $OPTS->{imaged} = 1;
    for $i (keys %{ $DATA->{ATOMS} }) {
        for $j (0 .. $#dim) {
            $OPTS->{scaled} = 0 if ($DATA->{ATOMS}{$i}{$dim[$j]} > 2);
            if (! $OPTS->{scaled}) {
                $OPTS->{imaged} = 0 if ($DATA->{ATOMS}{$i}{$dim[$j]} > $DATA->{"BOX BOUNDS"}[$j]{hi} ||
                                        $DATA->{ATOMS}{$i}{$dim[$j]} < $DATA->{"BOX BOUNDS"}[$j]{lo});
                last;
            }
        }
    }
    print "reading ";
    if (! $OPTS->{scaled}) {
        print "unscaled ";
    } else {
        print "scaled ";
    }
    if ($OPTS->{imaged}) {
        print "imaged ";
    } else {
        print "unimaged ";
    }
    print "coordinates...";
}

sub GetLammpsTrjType {
    my ($SELECT, $trjFile, $field, $OPTS) = @_;
    my ($i, $tmp1);

    print "Determining LAMMPS trajectory type...";
    for $i (keys %{ $SELECT }) {
	$tmp1->{$i} = $SELECT->{$i};
	last;
    }
    ParseLAMMPSTrj($OPTS, $trjFile, $tmp1, $field, \&determineIfScaled, undef, undef);
    print "Done\n";
}

sub ReadDataFile {
    my ($inFile) = $_[0];
    my (%DATA, $var, @val, $i);

    $var = "";
    open DATAFILE, $inFile or die "ERROR: Cannot open LAMMPS data file $inFile: $!\n";
    while (<DATAFILE>) {
	chomp;
	if ($_ =~ /^\s*(\d+) (atoms|bonds|angles|dihedrals|impropers)\s*$/) {
	    $DATA{TOTAL}{uc($2)} = $1;
	    $DATA{TOTAL}{count}++;
	} elsif ($_ =~ /^\s*(\d+) (atom|bond|angle|dihedral|improper) types\s*$/) {
	    $DATA{TYPE}{uc($2)} = $1;
	    $DATA{TYPE}{count}++;
	} elsif ($_ =~ /^\s*(\d+\.\d+)\s+(\d+\.\d+)\s+(x|y|z)lo (x|y|z)hi\s*$/) {
	    $DATA{BOX}{uc($2)}{lo} = $1;
	    $DATA{BOX}{uc($2)}{hi} = $2;
	    $DATA{BOX}{count}++;
	} elsif ($_ =~ /^Masses\s*$/) {
	    $var = "MASS";
	} elsif ($_ =~ /^(\w+) Coeffs/) {
	    $var = uc($1);
	} elsif ($_ =~ /^(Atoms|Velocities)\s*$/) {
	    $var = uc($1);
	} elsif (defined($var) and $_ =~ /^\s*(\d+)\s+(.+)$/) {
	    @val = split /\s+/, $2;
	    for $i (0 .. $#val) {
		$DATA{$var}{$1}{$i} = $val[$i];
	    }
	} elsif ($_ =~ /^\S+/) {
	    undef($var);
	}
    }

    for $var ("TYPE", "TOTAL", "BOX", "MASS", "ATOMS") {
	die "ERROR: Relevant header (" . uc($var) . ") not found when reading file!\n"
	    if (! exists($DATA{$var})); 
    }

    die "ERROR: Invalid format in data file\n"
	if ($DATA{TYPE}{count} < 5 || $DATA{TOTAL}{count} < 5 || $DATA{BOX}{count} < 3);

    return \%DATA;
	    
}

sub ParseLammpsTrajectoryFile {
    my ($fileName, $dType, $selection) = @_;
    my (%LINE, $timeStep, $totAtoms, %DATA, $atomC, $scaled);
    my (@dataVals, $patern, $BYTEINFO, $i);

    $timeStep = $atomC = 0;
    $scaled = 1;
    $patern = '^(\d+)\s+(\d+\s+\-?\d+\.?\d*\s+\-?\d+\.?\d*\s+\-?\d+\.?\d*)\s+';
    #$patern .= '(\-?\d+\s+\-?\d+\s+\-?\d+)?\s+';
    $patern .= '(\-?\d+\.?\d*e?\-?\d*)?\s$';

    #getByteOffset($selection, $fileName);
    open TRAJFILE, $fileName or die "ERROR: Cannot open trajectory file $fileName: $!\n";
    while (<TRAJFILE>) {
	chomp;
	if ($_ =~ /^\s*ITEM: TIMESTEP/) {
	    setZero(\%LINE);
	    $LINE{"tstep"} = 1;
	    if ($timeStep > 0) {
		delete $DATA{$timeStep} if (! exists($DATA{$timeStep}{"ATOMS"}) and lc($dType) ne "density");
	    }
	    $timeStep = $atomC = 0;
	} elsif ($_ =~ /^ITEM: NUMBER OF ATOMS/ and $timeStep > 0) {
	    setZero(\%LINE);
	    $LINE{"num_atoms"} = 1;
	} elsif ($_ =~ /^ITEM: BOX BOUNDS/ and $timeStep > 0) {
	    setZero(\%LINE);
	    $LINE{"box"} = 1;
	} elsif ($_ =~ /^ITEM: ATOMS/ and $timeStep > 0) {
	    setZero(\%LINE);
	    $LINE{"atoms"} = 1;
	} elsif ($_ =~ /^(\d+)$/ and $LINE{"tstep"}) {
	    setZero(\%LINE);
	    $timeStep = $1;
	    #print "GOT TIMESTEP $timeStep\n";
	} elsif ($_ =~ /^(\d+)$/ and $LINE{"num_atoms"} and $timeStep > 0) {
	    setZero(\%LINE);
	    $totAtoms = $1;
	    $DATA{$timeStep}{"TOTAL_ATOMS"} = $totAtoms;
	} elsif ($_ =~ /^(\-?\d+\.?\d*)\s+(\-?\d+\.?\d*)/ and $LINE{"box"} and $timeStep > 0) {
	    setZero(\%LINE);
	    $LINE{"box"} = 1;
	    push @{ $DATA{$timeStep}{"BOX"} }, (
						{
						    "lo" => $1,
						    "hi" => $2,
						}
						);
	} elsif ($LINE{"atoms"} and $timeStep > 0 and $_ =~ /^(\d+)\s(.+)/) {
	    $atomC = $1;
	    next if (defined($selection) and $selection ne "*" and ! exists($selection->{$atomC}));
	    @dataVals = split /\s/, $2;
	    if ($dType =~ /bgf|com|avg|rms/) {
		die "ERROR: Atom data not found for atom $atomC\n" if ($#dataVals < 3);
		$DATA{$timeStep}{"ATOMS"}{$atomC}{"TYPEID"} = shift @dataVals;
		$DATA{$timeStep}{"ATOMS"}{$atomC}{"XCOORD"} = $dataVals[0];
                $DATA{$timeStep}{"ATOMS"}{$atomC}{"YCOORD"} = $dataVals[1];
                $DATA{$timeStep}{"ATOMS"}{$atomC}{"ZCOORD"} = $dataVals[2];
		if ($scaled) {
		    for $i (0 .. 2) {
			if ($dataVals[$i] > 2) {
			    $scaled = 0;
			    last;
			}
		    }
		}
		if ($#dataVals > 4) {
		    $DATA{$timeStep}{"ATOMS"}{$atomC}{"XINDEX"} = $dataVals[3];
                    $DATA{$timeStep}{"ATOMS"}{$atomC}{"YINDEX"} = $dataVals[4];
                    $DATA{$timeStep}{"ATOMS"}{$atomC}{"ZINDEX"} = $dataVals[5];
		}
	    }

	    if ($dType =~ /vel/) {
		die "ERROR: Veleocity data not found for atom $atomC\n" if ($#dataVals < 2);
		$DATA{$timeStep}{"ATOMS"}{$atomC}{"XCOORD"} = $dataVals[0];
                $DATA{$timeStep}{"ATOMS"}{$atomC}{"YCOORD"} = $dataVals[1];
                $DATA{$timeStep}{"ATOMS"}{$atomC}{"ZCOORD"} = $dataVals[2];
	    }

	    if ($dType =~ /eng/) {
		$DATA{$timeStep}{"ATOMS"}{$atomC}{"ENERGY"} = pop @dataVals;
		$scaled = 0;
	    }
	}
    }
    close TRAJFILE or die "ERROR: Cannot close file $fileName: $!\n";

    die "ERROR: $fileName does not contain any valid informaton\n"
	if (!%DATA);
    
    return (\%DATA, $scaled);
}

sub setZero {
    my ($LINE, $exception) = @_;
    
    for (keys %{ $LINE }) {
	next
	    if ($exception and $exception =~ /\s+$_/);
	$LINE->{$_} = 0;
    }
}

sub numerically {
    ($a<=>$b);
}

sub ParseLAMMPSLogFile {
    my ($inFile, $selection, $saveFunc, $OUTDATA) = @_;
    my ($tStep, %DATA, $rem, $count);

    $tStep = -1;
    $count = 0;
    open LOGFILE, $inFile or die "ERROR: Cannot read from LAMMPS log file $inFile: $!\n";
    while (<LOGFILE>) {
	chomp;
	if ($_ =~ /^---------------- Step\s+(\d+)\s----- CPU =\s+(\d+\.\d+)/) {
	    if (keys %DATA) {
		$saveFunc->(\%DATA, $tStep, $OUTDATA);
		%DATA = ();
		$tStep = -1;
	    }
	    next if (defined($selection) and ! exists($selection->{$1}));
	    $tStep = $1;
	    $DATA{CPU} = $2;
	    $count++;
	} elsif ($tStep > -1 && $_ =~ /^(\S+)\s+\=\s+(\-?\d+\.\d+)(.+)/) {
	    $DATA{lc($1)} = $2;
	    $rem = $3;
	    while ($rem =~ /(\S+)\s+\=\s+(\-?\d+\.\d+)/g) {
		$DATA{lc($1)} = $2;
	    }
	} else {
            if (keys %DATA) {
                $saveFunc->(\%DATA, $tStep, $OUTDATA);
                %DATA = ();
	    }
	    $tStep = -1;
	}
    }

    $saveFunc->(\%DATA, $tStep, $OUTDATA) if (keys %DATA);
    close LOGFILE;
    
    die "ERROR: LAMMPS log file $inFile does not contain any valid information!\n"
	if (! $count);
}

sub CreateLAMMPSTrj {
    my ($ATOMS, $HEADERS, $OUTFILE) = @_;
    my ($i, $atomC, $atom, @dim, $index, $j, $bLen);

    @dim = ("TYPE", 
	    "XCOORD", "YCOORD", "ZCOORD", 
	    "XINDEX", "YINDEX", "ZINDEX",
	    "XVEL", "YVEL", "ZVEL");

    for $i ("TIMESTEP", "NUMBER OF ATOMS") {
	print $OUTFILE "ITEM: $i\n";
	for $j (@{ $HEADERS->{$i} }) {
	    print $OUTFILE "$j\n";
	}
    }
    # BOX
    print $OUTFILE "ITEM: BOX BOUNDS\n";
    for $i (@{ $HEADERS->{"BOX BOUNDS"} }) {
	print $OUTFILE "$i->{lo} $i->{hi}\n";
    }
    print $OUTFILE "ITEM: ATOMS\n";
    for $atomC (sort numerically keys %{ $ATOMS }) {
	$atom = \%{ $ATOMS->{$atomC} };
	print $OUTFILE "$atomC ";
        for $i (@dim) {
	    next if (! exists($atom->{$i}));
	    print $OUTFILE "$atom->{$i} ";
	}
	print $OUTFILE "\n";
    }
}    

sub getFormat {
    my ($type) = $_[0];
    my (@fields);

    $type = lc($type);
    if ($type eq "vel") {
	@fields = ("XVEL", "YVEL", "ZVEL");
    } elsif ($type eq "atom") {
	@fields = ("TYPE", "XCOORD", "YCOORD", "ZCOORD", "XINDEX", "YINDEX", "ZINDEX");
    } elsif ($type eq "energy") {
	@fields = ("TYPES", "ENERGY");
    }

    return \@fields;
}

sub getTotBytes {
    my ($inFile) = $_[0];
    my ($count) = 0;
    if (open(INFILE, $inFile)) {
	while (<INFILE>) {
	    chomp;
	    if ($_ =~ /^totalbytes: (\d+)/) {
		$count = $1;
	    }
	    close INFILE;
	}
    }
    return $count;
}

sub GetLammpsByteOffset {
    my ($selection, $trjFile, $junk) = @_;
    my ($counter, $isValid, $hasSelections, $count, $tot, $saveName, $offset, $last);
    my ($grepCmd) = "grep -b 'ITEM: TIMESTEP' $trjFile";

    $saveName = basename($trjFile);
    $saveName = "_byte_offset_${saveName}";
    print "Computing byte offset for LAMMPS trajectory $trjFile...";
    for $counter (keys %{ $selection }) {
	$selection->{$counter} = -1;
    }
    $counter = $isValid = $hasSelections = $count = $tot = 0;
    $tot = (-s $trjFile);
    $hasSelections = 1 if (keys %{ $selection });
    if (! -e $saveName) {
	open GREPCMD, "$grepCmd |" || die "ERROR: Cannot execute $grepCmd: $!\n";
    } else {
	$count = getTotBytes($saveName);
	if (($count != $tot) or ! open(GREPCMD, $saveName)) {
	    open GREPCMD, "$grepCmd |" || die "ERROR: Cannot execute $grepCmd: $!\n";
	} else {
	    print "reading file...";
	}
    }
    while (<GREPCMD>) {
        chomp;
	if ($_ =~ /^(\d+)/) {
	    $counter++;
	    $offset->{$counter} = $1;
            if (! $hasSelections || exists($selection->{$counter})) {
                $selection->{$counter} = $1;
                $isValid = 1;
            }
	    $last = $1;
        }
    }
    close GREPCMD;

    if (exists($selection->{"-1"})) {
	delete $selection->{-1};
	$isValid = 1;
    }
    if (! -e $saveName and open(BYTEOFFSET, "> $saveName")) {
	print BYTEOFFSET "totalbytes: $tot\n";
	for $count (sort { $a<=>$b} keys %{ $offset }) {
	    print BYTEOFFSET "$offset->{$count}\n";
	}
	close BYTEOFFSET;
    }
    $tot = $counter;
    $count = 0;
    for $counter (keys %{ $selection }) {
	$selection->{$counter} = ();
	$selection->{$counter}{OFFSET} = $offset->{$counter};
	if ($counter < $tot) {
	    $selection->{$counter}{LENGTH} = $offset->{($counter + 1)} - $offset->{$counter};
	} else {
	    $selection->{$counter}{LENGTH} = $last - $offset->{$counter};
	}
	$count++;
    }

    die "ERROR: No valid frame found in trajectory!\n" if (! $isValid);

    print "using $count of $tot snapshots...Done\n";
}

sub ParseLAMMPSTrj {
    my ($LOGDATA, $lammpsFile, $SELECT, $trjType, $doAnal, $printStr, $fileHandle) = @_;
    my ($inStr, $tStepData, $fields, $counter, $tot, $vals, $filesize, $frame);
    my ($atomC, $currPos, $start, $strLen, $i, $patern2);
    my ($patern1) = 'ITEM: TIMESTEP\n(\d+)\n' . 
	'ITEM: NUMBER OF ATOMS\n(\d+)\nITEM: BOX BOUNDS\n' .
	'(\-?\d+\.?\d*e?\-?\d*)\s+(\-?\d+\.?\d*e?\-?\d*)\n' .
	'(\-?\d+\.?\d*e?\-?\d*)\s+(\-?\d+\.?\d*e?\-?\d*)\n' .
	'(\-?\d+\.?\d*e?\-?\d*)\s+(\-?\d+\.?\d*e?\-?\d*)\n';
    my ($patern2) = '(\-?\d+\.?\d*e?\-?\d*)';

    print "${printStr}Calculating time remaining\r" if (defined($printStr));
    $strLen = length("Calculating time remaining");
    $fields =  getFormat($trjType);
    $patern2 = '(\d+)';
    for $i (0 .. $#{ $fields }) {
	$patern2 .= '\s+(\-?\d+\.?\d*e?\-?\d*)';
    }
    $patern2 .= '.*\n';

    open TRAJFILE, $lammpsFile or die "ERROR: Cannot open trajectory file $lammpsFile: $!\n";
    $tot = scalar keys %{ $SELECT };
    $start = time();
    $currPos = 0;

    for $i (sort numerically keys %{ $SELECT }) {
	if (read(TRAJFILE, $inStr, $SELECT->{$i}{LENGTH}, $SELECT->{$i}{OFFSET})) {
	    $currPos++;
	    if ($inStr =~ qr/$patern1/) {
		$tStepData->{FRAME} = $currPos;
		$tStepData->{TIMESTEP}[0] = $1;
		$tStepData->{"NUMBER OF ATOMS"}[0] = $2;
		$tStepData->{"BOX BOUNDS"}[0]{lo} = $3;
                $tStepData->{"BOX BOUNDS"}[0]{hi} = $4;
                $tStepData->{"BOX BOUNDS"}[1]{lo} = $5;
                $tStepData->{"BOX BOUNDS"}[1]{hi} = $6;
                $tStepData->{"BOX BOUNDS"}[2]{lo} = $7;
                $tStepData->{"BOX BOUNDS"}[2]{hi} = $8;
	    } 
	    if ($inStr =~ /ITEM: ATOMS\n/g and exists($tStepData->{TIMESTEP})) {
		while ($inStr =~ /$patern2/g) {
		    $atomC = $1;
		    last if ($atomC > $tStepData->{"NUMBER OF ATOMS"}[0]);
		    for $counter (1 .. $#{ $fields }) {
			$tStepData->{ATOMS}{$atomC}{ $fields->[$counter] } = eval('$' . $counter);
		    }
		}
		$doAnal->($tStepData, $LOGDATA, $fileHandle);
		$strLen = PrintProgress($currPos, $tot, $start, $printStr);
	    }
	    undef($tStepData);
	}
    }
    close TRAJFILE or die "ERROR: Cannot close $lammpsFile: $!\n";

    printf "$printStr%-${strLen}s\n", "Done" if (defined($printStr));
}
 
sub CreateInputFile {
    my ($parms) = $_[0];
    my ($header, $middle, $tail, $fileName);

    $fileName = "in." . $parms->{PARMS}{SUFFIX};
    $header = createInputHeader($parms, $parms->{PARMS}{FFTYPE}, $parms->{PARMS}{NUM_FILES});
    $middle = createMiddle($parms, $parms->{PARMS}{SOLUTE}, $parms->{PARMS}{SUFFIX});
    $tail = readTail($parms->{PARMS}{SOLUTE}, $parms);
    open LAMMPSINPUT, "> $fileName" || die "ERROR: Cannot create $fileName:$!\n";
    print LAMMPSINPUT "$header\n$middle\n$tail";
    close LAMMPSINPUT;
    open SINGLEPOINT, "> ${fileName}_singlepoint" || die "ERROR: Cannot create ${fileName}_singlepoint:$!\n";
    print SINGLEPOINT "$header\n$middle\nrun 0";
    close SINGLEPOINT;
}

sub createMiddle {
    my ($parms, $soluteAtms, $suffix) = @_;
    my ($middle, $numAtms, $i, $j);
    
    $middle = "read_data       data.${suffix}\n\n$parms->{PARMS}{OFF_DIAG}";
    $middle .= "pair_modify     mix $parms->{PARMS}{mix_rule}\n";
    $middle .= "neighbor        2.0 multi\n";
    $middle .= "neigh_modify    every 2 delay 4 check yes\n";
    $middle .= "thermo_style    multi\n\n";
    $middle .= "variable        input index in.${suffix}\n";
    $middle .= "variable        sname index $suffix\n";
    $numAtms = scalar keys %{ $soluteAtms };
    if ($numAtms > 0) {
	$middle .= "variable        sAtoms index $numAtms\n";
	$middle .= "group           solute id <> 1 \${sAtoms}\n";
	$middle .= "group           solvent subtract all solute\n";
    }

    return $middle;
}

sub readTail {
    my ($soluteAtms, $parms) = @_;
    my ($hasSolute, $outStr, $inStr, $inputFile, $shakeOpts);
    
    $hasSolute = 1;
    $hasSolute = 0 if (! keys %{ $soluteAtms });
    if (exists($parms->{PARMS}{INPUTTYPE}))  {
	$inputFile = "/ul/tpascal/scripts/dat/LAMMPS/in.lammps." . lc($parms->{PARMS}{INPUTTYPE});
    } else {
        $inputFile = "/ul/tpascal/scripts/dat/LAMMPS/in.lammps.full";
    }
    $shakeOpts = "";
    $shakeOpts = $parms->{PARMS}{SHAKE_MASS} . $parms->{PARMS}{SHAKE_ANGLE};
    open TAIL, $inputFile || die "ERROR: Cannot open $inputFile: $!\n";
    while (<TAIL>) {
	chomp;
	$inStr = $_;
	if (! $hasSolute) {
	    next if $inStr =~ /solute/;
	    $inStr =~ s/solvent/all/g;
	}
	if ($inStr =~ /shakeOpts/) {
	     if ($shakeOpts ne "") {
		$inStr =~ s/shakeOpts/$shakeOpts/;
	     } else {
		$inStr = "";
	     }
	}
	$inStr = "" if ($inStr =~ /unfix\s+shakeH/ && $shakeOpts eq "");
	$outStr .= "$inStr\n";
    }

    close TAIL;
    return $outStr;
}

sub createInputHeader {
    my ($parms, $ffType, $numFiles) = @_;
    my ($header, $i, $parm, $j);

    $header = "units           real\natom_style      full\nboundary        p p p\n";
    if ($parms->{PARMS}{isTriclinic}) {
	$header .= "special_bonds   0 0 0.5\n\n";
	$header .= "pair_style      lj/coul long long 10.0";
    } elsif ($numFiles == 1 && $ffType == 1) { #amber
	$header .= "special_bonds   amber\n\n";
	$header .= "pair_style      lj/charmm/coul/long/opt 9.0 10.0";
    } elsif ($numFiles == 1 && $ffType == 2) { #charmm
	$header .= "special_bonds   charmm\n\n";
	$header .= "pair_style      lj/charmm/coul/long/opt 9.0 10.0";
#    } elsif ($numFiles == 1 && $ffType == 4) { #dreiding
#	$header .= "special_bonds   0.0 0.0 0.5\n\n";
#	$header .= "pair_style      hybrid/overlay dreiding/hb 5.0 90 lj/cut/coul/long " . 
#		    "$parms->{PARMS}{cut_coul} $parms->{PARMS}{cut_vdw}";
    } elsif ($numFiles == 1 && $ffType == 3) { #MESODNA
	$header .= "special_bonds   amber\n\n";
	$header .= "pair_style      hybrid yukawa 3.0 6.0  dreiding/hb 6.0 120.0 morse/opt 12.0";
    } else {
	$header .= "special_bonds   ";
	for $i (2 .. 4) {
	    $header .= sprintf("%2.1f ", $parms->{PARMS}{"scale_cou_1${i}"});
	}
	if (! $parms->{PARMS}{same_scale}) {
	    for $i (2 .. 4) {
		$header .= sprintf("%2.1f ", $parms->{PARMS}{"scale_vdw_1${i}"});
	    }
	}
	$parms->{PARMS}{dielectric} += 0;
	$header .= "\ndielectric      $parms->{PARMS}{dielectric}";
	$header .= "\n\npair_style      ";
	if (scalar keys %{ $parms->{VDW}{TYPE} } > 1) {
	    if (exists($parms->{PARMS}{HYBRID_OVERLAP})) {
		$header .= "hybrid/overlay ";
	    } else {
		$header .= "hybrid ";
	    }
	}elsif (scalar keys %{ $parms->{VDW}{TYPE } } == 0) {
	    $header .= "none";
	}
	for $j (keys %{ $parms->{VDW}{TYPE} }) {
	    $header .= "$j ";
	}
    }
    
    for $i ("BONDS", "ANGLES", "TORSIONS", "INVERSIONS") {
	$parm = lc(substr($i, 0, -1));
	$parm = "dihedral" if ($parm eq "torsion");
	$parm = "improper" if ($parm eq "inversion");
	$header .= sprintf("\n%-16s", $parm . "_style");
	if (scalar keys %{ $parms->{$i}{TYPE} } > 1) {
	    $header .= "hybrid ";
	} elsif (scalar keys %{ $parms->{$i}{TYPE } } == 0) {
	    $header .= "none";
	}
	for $j (keys %{ $parms->{$i}{TYPE} }) {
	    if ($i eq "TORSIONS" and $j eq "charmm " and $header !~ /charmm/) {
		$header  .= "shifted ";
	    } else {
		$header .= "$j ";
	    }
	}
    }
    if ($parms->{PARMS}{isTriclinic}) {
	$header .= sprintf("\nkspace_style    ewald/n %-8.4g\n",$parms->{PARMS}{coul_accuracy});
    } elsif ($header =~ /coul.*long/) {
	$header .= sprintf("\nkspace_style    pppm %-8.4g\n",$parms->{PARMS}{coul_accuracy});
    } else {
	$header .= "\nkspace_style    none\n";
    }

    return $header;
}

sub ConvertLammpsBox {
    my ($box) = $_[0];
    my ($generalBox, $i, @dim, %BOX);

    @dim = ("XCOORD", "YCOORD", "ZCOORD");

    for $i (0 .. 2) {
        $BOX{$dim[$i]}{lo} = $box->[$i]{lo};
        $BOX{$dim[$i]}{hi} = $box->[$i]{hi};
        $BOX{$dim[$i]}{len} = $box->[$i]{hi} - $box->[$i]{lo};
	$BOX{$dim[$i]}{CENTER} = $BOX{$dim[$i]}{len}/2 - $box->[$i]{lo};
    }
                                                                                                                
    return \%BOX;
}   
1;

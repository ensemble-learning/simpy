package Packages::LAMMPS;

require Exporter;
use Packages::General qw(PrintProgress GetTime);
use File::Basename qw(basename);
our (@ISA, @EXPORT, $VERSION, @EXPORT_OK);
use strict;

my $scripts_dir = "/ul/tpascal/scripts"; #change here as necessary
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(ParseLammpsTrajectoryFile ReadDataFile ParseLAMMPSLogFile CreateLAMMPSTrj 
                ParseLAMMPSTrj CreateInputFile GetLammpsByteOffset GetLammpsTrjType 
                ConvertLammpsBox ParseLAMMPSTrjSave);
$VERSION = "1.00";

sub determineIfScaled {
    my ($DATA, $OPTS, $tmp2) = @_;
    my ($i, $j);
    my (@dim) = ("XCOORD", "YCOORD", "ZCOORD");

    $OPTS->{ISIFT} = 0;
    $OPTS->{scaled} = $OPTS->{imaged} = 1;
    for $i (keys %{ $DATA->{ATOMS} }) {
	$OPTS->{ISIFT} = 1 if ($DATA->{ATOMS}{$i}{SYZ});
        for $j (0 .. $#dim) {
            $OPTS->{scaled} = 0 if ($DATA->{ATOMS}{$i}{$dim[$j]} > 2 or $DATA->{ATOMS}{$i}{$dim[$j]} < -2);
            if (! $OPTS->{scaled}) {
                $OPTS->{imaged} = 0 if ($DATA->{ATOMS}{$i}{$dim[$j]} > $DATA->{"BOX BOUNDS"}[$j]{hi} or
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
	} elsif ($tStep > -1 && $_ =~ /^(\S+)\s+\=\s+(\-?\d+\.?\d*E?\-?\d*)(.+)/) {
	    $DATA{lc($1)} = $2;
	    $rem = $3;
	    while ($rem =~ /(\S+)\s+\=\s+(\-?\d+\.?\d*E?\-?\d*)/g) {
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
    print $OUTFILE "ITEM: ATOMS id xu yu zu\n";
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
	@fields = ("TYPE", "XCOORD", "YCOORD", "ZCOORD", "XINDEX", "YINDEX", "ZINDEX", "CHARGE");
    } elsif ($type eq "ift") {
        #@fields = ("TYPE", "XCOORD", "YCOORD", "ZCOORD", "XINDEX", "YINDEX", 
			#"ZINDEX", "SXX", "SYY", "SZZ", "SXY", "SXZ", "SYZ");
	@fields = ("TYPE", "XCOORD", "YCOORD", "ZCOORD", "SXX", "SYY", "SZZ", "SXY", "SXZ", "SYZ");
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
		last;
	    }
	}
	close INFILE;
    }
    return $count;
}

sub GetLammpsByteOffset {
    my ($selection, $trjFile, $junk) = @_;
    my ($counter, $isValid, $hasSelections, $count, $tot, $saveName, $offset, $last, $bSelect, $shouldWrite);
    my ($grepCmd) = "grep -b 'ITEM: TIMESTEP' $trjFile";

    $shouldWrite = 1;
    $saveName = basename($trjFile);
    $saveName = ".byte.offset.${saveName}";
    print "Computing byte offset for LAMMPS trajectory $trjFile...";
    for $counter (keys %{ $selection }) {
	$bSelect->{$counter}{start} = -1;
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
	    $shouldWrite = 0;
	}
    }
    while (<GREPCMD>) {
        chomp;
	if ($_ =~ /^(\d+)/) {
	    $counter++;
	    $offset->{$counter} = $1;
            if (! $hasSelections || exists($selection->{$counter})) {
                $bSelect->{$counter}{start} = $1;
                $isValid = 1;
            }
	    $last = $1;
        }
    }
    close GREPCMD;

    if (exists($selection->{"-1"})) {
	$bSelect->{-1}{start} = $last;
	$isValid = 1;
    }
    if ($shouldWrite and open(BYTEOFFSET, "> $saveName")) {
	print BYTEOFFSET "totalbytes: $tot\n";
	for $count (sort { $a<=>$b} keys %{ $offset }) {
	    print BYTEOFFSET "$offset->{$count}\n";
	}
	close BYTEOFFSET;
    }
    $tot = $counter;
    $count = 0;
    for $counter (keys %{ $selection }) {
        if ($bSelect->{$counter}{start} == -1) {
	    delete $bSelect->{$counter};
	    delete $selection->{$counter};
	} else {
	    if (exists($offset->{ ($counter + 1) })) {
		$bSelect->{$counter}{end} = $offset->{ ($counter + 1) } -1;
	    } else {
		$bSelect->{$counter}{end} = $last;
	    }
	    $count++;
	}
    }

    die "ERROR: No valid frame found in trajectory!\n" if (! $isValid);

    $count = $tot if (! $hasSelections);
    print "using $count of $tot snapshots...Done\n";
    %{ $selection } = %{ $bSelect };
}

sub getNewFormat {
    my ($fmtStr, $oldFmt) = @_;
    my (@fields, $curr);
    while ($fmtStr =~ /(\w+)/g) {
	$curr = uc $1;
	if ($curr =~ /type/i) {
	    push @fields, "TYPE";
	} elsif ($curr =~ /^(x|y|z)(u|s)?/i) {
	    push @fields, "${1}COORD";
	} elsif ($curr =~ /^i(x|y|z)/i) {
	    push @fields, "${1}INDEX";
	} elsif ($curr =~ /^v(x|y|z)/i) {
            push @fields, "${1}VEL";
	} elsif ($curr =~ /^(s)(x|y|z)(x|y|z)/i) {
            push @fields, "${1}${2}${3}";
	} elsif ($curr =~/q/i) {
	    push @fields, "CHARGE";
	} elsif ($curr !~ /^id/i) {
	    push @fields, $curr;
	}
    }
    if (! @fields) {
	@fields = @{ $oldFmt };
    }
    return \@fields;
}

sub ParseLAMMPSTrj {
    my ($LOGDATA, $lammpsFile, $SELECT, $trjType, $doAnal, $printStr, $fileHandle) = @_;
    my ($inStr, $tStepData, $fields, $counter, $tot, $field, $filesize, $frame, $totF);
    my ($atomC, $currPos, $start, $strLen, $tStep, %header, $rec, $coords, $i, @vals);

    print "${printStr}Calculating time remaining\r" if (defined($printStr));
    $strLen = length("Calculating time remaining");
    $fields =  getFormat($trjType);
    $totF = scalar(@{ $fields });
    %header = ("TIMESTEP"=>1,"NUMBER OF ATOMS"=>1,"BOX BOUNDS"=>1,"ATOMS"=>1);
    #getByteOffset($SELECT, $lammpsFile);

    open TRAJFILE, $lammpsFile or die "ERROR: Cannot open trajectory file $lammpsFile: $!\n";
    $tot = scalar keys %{ $SELECT };
    $start = time();
    $currPos = 0;

    for $i (sort numerically keys %{ $SELECT }) {
        next if ( ! seek(TRAJFILE, $SELECT->{$i}{start}, 0));
	$currPos++;
	while (<TRAJFILE>) {
	    chomp;
	    $inStr = $_;
	    study;
	    if ($inStr =~ /^ITEM:\s+(TIMESTEP|NUMBER OF ATOMS|BOX BOUNDS|ATOMS)(.*)$/io) {
		if (exists($header{$1})) {
		    $field = $1;
		    if ($field eq "TIMESTEP" && keys %{ $tStepData }) {
			$tStepData->{FRAME} = $currPos;
			$doAnal->($tStepData, $LOGDATA, $fileHandle);
			$strLen = PrintProgress($currPos, $tot, $start, $printStr);
			undef($tStepData);
			last;
		    }elsif ($field eq "ATOMS" and defined($2)) {
			$fields = getNewFormat($2, $fields);
			$totF = scalar(@{ $fields });
		    }
		} else {
		    undef($field);
		}
	    } elsif (defined($field) && $inStr =~ /^\s*(\-?\d+\.?\d*e?[\-|\+]?\d*)\s*(\-?\d*.*)$/o) {
		$atomC = $1;
		$coords = $2;
		if ($field eq "BOX BOUNDS") {
		    $rec = (
			    {
				"lo" => $1,
				"hi" => $2,
			    }
			    );
		    push @{ $tStepData->{$field} }, $rec;
		} elsif ($field =~ /^ATOMS/o) {
		    @vals = split /\s+/, $coords;
		    $counter = 0;
		    while ($counter < $totF) {
			$tStepData->{ATOMS}{$atomC}{ $fields->[$counter] } = $vals[$counter];
			$counter++;
		    }
		    #while ($coords =~ /\s(\-?\d+\.?\d*e?[\-|\+]?\d*)/g && ($counter <= $#{ $fields })) {
			#$tStepData->{ATOMS}{$atomC}{ $fields->[$counter] } = $1;
			#$counter++;
		    #}
		} else {
		    $tStepData->{$field}[0] = $1;
		}
	    }
	}
    }
    close TRAJFILE or die "ERROR: Cannot close $lammpsFile: $!\n";

    if (keys %{ $tStepData }) {
        $doAnal->($tStepData, $LOGDATA, $fileHandle);
    }
    $i = GetTime(time() - $start);
    printf "$printStr%-${strLen}s\n", "${i}s elapsed..Done" if (defined($printStr));
}

sub ParseLAMMPSTrj_new {
    my ($LOGDATA, $lammpsFile, $SELECT, $trjType, $doAnal, $printStr, $fileHandle) = @_;
    my ($strLen, $fields, $totF, $tot, $start, $currPos, $headerPtn, $aC, $field, $totL);
    my ($i, $j, $k, $tStepData, $bLen, $buf, $counter, @tmp, @vals, $atomC, $aT, $vArrayTot);

    print "${printStr}Calculating time remaining\r" if (defined($printStr));
    $strLen = length("Calculating time remaining");
    $fields =  getFormat($trjType);
    $totF = scalar(@{ $fields });
    open TRAJFILE, $lammpsFile or die "ERROR: Cannot open trajectory file $lammpsFile: $!\n";
    $tot = scalar keys %{ $SELECT };
    $start = time();
    $currPos = $vArrayTot = 0;
    $headerPtn = qr/^ITEM: TIMESTEP\n(\d+)\nITEM: NUMBER OF ATOMS\n(\d+)\nITEM: BOX BOUNDS\n([^\n]+)\n([^\n]+)\n([^\n]+)\nITEM: ATOMS(?:[^\n]*)\n/;
    $field = qr/(\S+)/;
    for $i (sort numerically keys %{ $SELECT }) {
	$tStepData = ();
	$bLen = $SELECT->{$i}{end} - $SELECT->{$i}{start};
	next if ( ! sysseek(TRAJFILE, $SELECT->{$i}{start}, 0));	
        next if ( ! sysread(TRAJFILE, $buf, $bLen,0));
	next if ($buf !~ /$headerPtn/g);
	$currPos++;
	$tStepData->{TIMESTEP}[0] = $1;
	$tStepData->{"NUMBER OF ATOMS"}[0] = $2;
	$aT = $2;
	$counter = 0;
	for $i ($3, $4, $5) {
 	    @tmp = split /\s+/, $i;
	    $tStepData->{"BOX BOUNDS"}[$counter]{lo} = $tmp[0];
	    $tStepData->{"BOX BOUNDS"}[$counter]{hi} = $tmp[1];
	    $counter++;
	}
	$counter = -1;
	$aC = 1;
	if (! $vArrayTot) {
	    $buf =~ /(.+)\n/g;
	    @tmp = split /\s+/, $1;
	    $totL = scalar(@tmp); # total number of fields per line
	    @vals = (@tmp, split /\s+/, $');
	    $vArrayTot = scalar(@vals) - $totL;
	    $k = $totF;
	    $k = $totL if ($totL < $k);
	} else {
	    @vals = split /\s+/, $';
	}

	$j = 0;
	while ($j <= $vArrayTot) {
	    $aC++;
	    #last if ($aC > $aT);
	    $atomC = $vals[$j];
	    $counter = 0;
	    while ($counter < $k) {
		$tStepData->{ATOMS}{$atomC}{ $fields->[$counter] } = $vals[($j + $counter + 1)];
		$counter++;
	    }
	    $j += $totL;
	}
	$tStepData->{FRAME} = $currPos;
	$doAnal->($tStepData, $LOGDATA, $fileHandle);
	$strLen = PrintProgress($currPos, $tot, $start, $printStr);
    }
    close TRAJFILE or die "ERROR: Cannot close $lammpsFile: $!\n";

    $i = GetTime(time() - $start);
    printf "$printStr%-${strLen}s\n", "${i}s elapsed..Done" if (defined($printStr));
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
    my ($middle, $numAtms, $i, $j, $totAtms);
    
    $parms->{HYBRID_VALENCE} = "" if (! exists($parms->{HYBRID_VALENCE}));
    $middle = "read_data       data.${suffix}\n\n$parms->{PARMS}{OFF_DIAG}";
    $middle .= "pair_modify     mix $parms->{PARMS}{mix_rule}\n" if (exists($parms->{PARMS}{mix_rule}));
    $middle .= "neighbor        2.0 multi\n";
    $middle .= "neigh_modify    every 2 delay 4 check yes\n";
    $middle .= "thermo_style    multi\nthermo_modify	line multi format float %14.6f\n";
    $middle .= "$parms->{HYBRID_VALENCE}";
    $middle .= "variable        input index in.${suffix}\n";
    $middle .= "variable        sname index $suffix\n";
    $numAtms = scalar keys %{ $soluteAtms };
    $totAtms = $parms->{PARMS}{NUM_ATOMS};
    if ($numAtms > 0 and $numAtms < $totAtms) {
	$middle .= "variable        sAtoms index $numAtms\n";
	$middle .= "group           solute id <> 1 \${sAtoms}\n";
	$middle .= "group           solvent subtract all solute\n";
    }
    if($parms->{PARMS}{USE_HBOND}) {
        $middle .= "\ncompute   hb all pair hbond/dreiding/lj\n" if ($parms->{PARMS}{USE_HBOND} == 1);
        $middle .= "\ncompute   hb all pair hbond/dreiding/morse\n" if ($parms->{PARMS}{USE_HBOND} == 2);
        $middle .= "variable    E_hbond equal c_hb[1]\n";
	$middle .= "thermo_style 	custom etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong v_E_hbond press vol\n";
	$middle .= "thermo_modify   line multi\n";
    }
    if(exists($parms->{PARMS}{REAX})) {
	$middle .= "\nfix	charge all qeq/reax 1 0.0 10.0 1.0e-6 reax/c\n";
    } elsif(exists($parms->{QEq}) and keys %{ $parms->{QEq} }) {
	$middle .= "\nfix	charge all qeq/reax 1 0.0 10.0 1.0e-06 param.qeq\n";
    }
    return $middle;
}

sub readTail {
    my ($soluteAtms, $parms) = @_;
    my ($hasSolute, $outStr, $inStr, $inputFile, $shakeOpts, $solAtms, $numAtms);
    
    $hasSolute = 1;
    if (keys %{ $soluteAtms }) {
	$solAtms = scalar(keys %{ $soluteAtms });
	$numAtms = $parms->{PARMS}{NUM_ATOMS};
	$hasSolute = 0 if ($numAtms == $solAtms);
    } else {
	$hasSolute = 0;
    }

    foreach $inputFile (@{ $parms->{PARMS}{INPUTLOC} }) {
	$shakeOpts = "";
	$shakeOpts = $parms->{PARMS}{SHAKE_MASS} . $parms->{PARMS}{SHAKE_ANGLE} if ($parms->{PARMS}{OPTIONS}{SHAKE});
	open TAIL, $inputFile || die "ERROR: Cannot open $inputFile: $!\n";
	while (<TAIL>) {
	    chomp;
	    $inStr = $_;
	    if (! $hasSolute) {
		next if $inStr =~ /solute|restraint|unfix\s+3/;
		$inStr =~ s/solvent/all/g;
	    }
	    if ($inStr =~ /shakeOpts/) {
		$inStr =~ s/ all / $parms->{PARMS}{OPTIONS}{SHAKE_WHAT} / if($hasSolute && exists($parms->{PARMS}{OPTIONS}{SHAKE_WHAT}));
		if ($shakeOpts ne "") {
		    $inStr =~ s/shakeOpts/$shakeOpts/;
	 	} else {
		    $inStr = "";
	 	}
	    }
	    $inStr = "thermo_style	custom etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong v_E_hbond press vol" 
		if ($inStr =~ /^thermo_style/ and $parms->{PARMS}{USE_HBOND});
	    $inStr =~ s/all npt.*$/all npt temp 300.0 300.0 100.0 x 1.0 1.0 2000.0 y 1.0 1.0 2000.0 couple xy/ 
		if ($parms->{PARMS}{OPTIONS}{PERIODICITY} == 2);
	    $inStr =~ s/wall_str/fix zconfine all wall\/lj93 zlo EDGE 1.0 1.0 2.5 zhi EDGE 1.0 1.0 2.5 units box/
		if ($parms->{PARMS}{OPTIONS}{PERIODICITY} == 2);
	    $inStr =~ s/wall_str// if($parms->{PARMS}{OPTIONS}{PERIODICITY} != 2);
	    $inStr = "" if ($inStr =~ /unfix\s+shakeH/ && $shakeOpts eq "");
	    $outStr .= "$inStr\n";
	}

	close TAIL;
    }
    return $outStr;
}

sub createInputHeader {
    my ($parms, $ffType, $numFiles) = @_;
    my ($header, $i, $parm, $j);

    $parms->{PARMS}{dielectric} += 0;
    $header = "units           real\natom_style      full\n";
    if($parms->{PARMS}{OPTIONS}{PERIODICITY} == 0) { 
	$header .= "boundary        f f f\n";
    }elsif ($parms->{PARMS}{OPTIONS}{PERIODICITY} == 2) {
	$header .= "boundary        p p f\n";
    }else {
	$header .= "boundary        p p p\n";
    }
    $header .= "dielectric      $parms->{PARMS}{dielectric}\n";
    #$header .= "newton          off\n" if ($parms->{PARMS}{HAS_HBONDS});

    #if ($parms->{PARMS}{isTriclinic}) {
	#$header .= "special_bonds   lj/coul 0 0 0.5\n\n";
	#$header .= "pair_style      lj/coul long long 10.0";
#    } elsif ($ffType == 1) { #amber
    if ($ffType == 1 and ! $parms->{PARMS}{is_tip4p}) {
	$header .= "special_bonds   amber\n\n";
	$header .= "pair_style      lj/charmm/coul/long/opt 9.0 10.0";
    } elsif ($ffType == 2 && ! $parms->{PARMS}{is_tip4p}) { #charmm
	$header .= "special_bonds   charmm\n\n";
	$header .= "pair_style      lj/charmm/coul/long/opt 9.0 10.0";
#    } elsif ($numFiles == 1 && $ffType == 4) { #dreiding
#	$header .= "special_bonds   0.0 0.0 0.5\n\n";
#	$header .= "pair_style      hybrid/overlay dreiding/hb 5.0 90 lj/cut/coul/long " . 
#		    "$parms->{PARMS}{cut_coul} $parms->{PARMS}{cut_vdw}";
    } elsif ($ffType == 3) { #MESODNA
	$header .= "special_bonds   amber\n\n";
	$header .= "pair_style      hybrid yukawa 3.0 6.0  dreiding/hb 6.0 120.0 morse/opt 12.0";
    } elsif ($ffType == 5) { #reax
	#$header .= "\npair_style      reax 10.0 1.0e-5";
	$header .= "\npair_style      reax/c NULL";
    } elsif ($ffType == 6  and ! $parms->{PARMS}{is_tip4p}) { #opls
        $header .= "special_bonds   lj/coul 0.0 0.0 0.5\n\n";
        $header .= "pair_style      lj/charmm/coul/long/opt 9.0 10.0";	
    } else {
	$header .= "special_bonds   lj";
	$header .= "/coul" if ($parms->{PARMS}{same_scale});
	$header .= " ";
	for $i (2 .. 4) {
	    $header .= sprintf("%2.1f ", $parms->{PARMS}{"scale_cou_1${i}"});
	}
	if (! $parms->{PARMS}{same_scale}) {
	    $header .= "coul ";
	    for $i (2 .. 4) {
		$header .= sprintf("%2.1f ", $parms->{PARMS}{"scale_vdw_1${i}"});
	    }
	}
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
	    $header .= "$j $parms->{VDW}{TYPE}{$j} ";
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
	    #if ($i eq "TORSIONS" and $j eq "charmm " and $header !~ /charmm/) {
		#$header  .= "shifted ";
	    #} else {
		$header .= "$j ";
	    #}
	}
    }
    if ($parms->{PARMS}{isTriclinic} and $parms->{PARMS}{OPTIONS}{PERIODICITY} and ($parms->{PARMS}{HAS_CHARGE} or $parms->{PARMS}{OPTIONS}{VDW_EWALD}) ) {
	$header .= sprintf("\nkspace_style    ewald/n %-8.4g\n",$parms->{PARMS}{coul_accuracy});
    } elsif ($header =~ /tip4p/ and $parms->{PARMS}{OPTIONS}{PERIODICITY}) {
        $header .= sprintf("\nkspace_style    pppm/tip4p %-8.4g\n",$parms->{PARMS}{coul_accuracy});
    } elsif ($header =~ /coul.*long/ and $parms->{PARMS}{OPTIONS}{PERIODICITY} and ! $parms->{PARMS}{OPTIONS}{VDW_EWALD}) {
	$header .= sprintf("\nkspace_style    pppm %-8.4g\n",$parms->{PARMS}{coul_accuracy});
    } elsif (($header =~ /coul.*long/ and $parms->{PARMS}{OPTIONS}{PERIODICITY}) or $parms->{PARMS}{OPTIONS}{VDW_EWALD}) {
	$header .= sprintf("\nkspace_style    ewald/n %-8.4g\n",$parms->{PARMS}{coul_accuracy});
    } else {
	$header .= "\nkspace_style    none\n";
    }
    $header .= "kspace_modify	slab 2.0\n" if($parms->{PARMS}{OPTIONS}{PERIODICITY} == 2 and $header =~ /coul\/long/);
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

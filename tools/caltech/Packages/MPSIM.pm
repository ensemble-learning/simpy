package Packages::MPSIM;

require Exporter;
use strict;

our(@ISA, @EXPORT, $VERSION, @EXPORT_OK);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(parseMPSimFF);
$VERSION = "1.00";

my $scripts_dir = "/ul/tpascal/scripts"; #change here as necessary
sub updateMass {
    my ($currMass, $modMass) = @_;

    if ($modMass =~ /(\d+\.\d+)/) {
	$currMass = $1;
    }

    return $currMass;
}

sub loadCerius2LammpsConverter {
    my ($inFile) = "${scripts_dir}/dat/ceriusParms2Lammps.perldata";
    my (%SHASH);
    scalar eval `cat $inFile` or die "Cannot recreate data in file $inFile: $! $@\n";
    return \%SHASH;
}

sub loadMpsim2CeriusConverter {
    my ($inFile) = "${scripts_dir}/dat/mpsimTypes2Cerius.perldata";
    my (%SHASH);
    scalar eval `cat $inFile` or die "Cannot recreate data in file $inFile: $! $@\n";
    return \%SHASH;
}

sub getLammpsOpts {
    my ($parmName, $parmType, $SHASH, $istip4p) = @_;
    my (%RET);
    my ($searchStr) = lc($parmType) . "_" . lc($parmName);
    
    $istip4p = 0 if (! defined($istip4p));
    $searchStr .= "_tip4p" if ($parmType eq "vdw" and $istip4p);
    if (exists($SHASH->{$searchStr})) {
	%RET = %{ $SHASH->{$searchStr} };
    } else {
	$RET{name} = lc($parmName);
	$RET{opts} = "";
    }
    return \%RET;
}

sub getCeriusTypes {
    my ($mpsimName, $parmType, $SHASH)  = @_;
    my ($ret);

    my ($searchStr) = lc($parmType) . "_" . lc($mpsimName);

    return $SHASH->{$searchStr};
}

sub getEquilVal {
    my ($parms, @vals) = @_;
    my ($equil_val, $i, $tmp, $currect_parm, $isVal);
    
    $tmp = $parms;
    $isVal = 1;
    for $i (@vals) {
	if (exists($tmp->{$i})) {
	    $tmp = $tmp->{$i};
	} else {
	    $isVal = 0;
	    last;
	}
    }
    if (! $isVal) {
	@vals = reverse @vals;
	$tmp = $parms;
	$isVal = 1;
	for $i (@vals) {
	    if (exists($tmp->{$i})) {
		$tmp = $tmp->{$i};
	    } else {
		$isVal = 0;
		last;
	    }
	}
    }

    if ($isVal) {
	for $i (sort {$a<=>$b} keys %{ $tmp }) {
	    $tmp = $tmp->{$i};
	    last;
	}
	$equil_val = $tmp->{VALS}[1];
	$equil_val = $tmp->{VALS}[2] if ($tmp->{TYPE} eq "MORSE" and ! exists($tmp->{CTYPE}));
    }
    
    return $equil_val;
}

sub findDuplicate {
    my ($parm, $currType) = @_;
    my ($i, $typeCounter);
                                                                                                                                 
    $typeCounter = 0;
    $currType = lc($currType);
    for $i (keys %{ $parm }) {
        if (lc($parm->{$i}{TYPE}) eq $currType) {
            $typeCounter = $i;
            last;
        }
    }
                                                                                                                                 
    return $typeCounter;
}

sub parseMPSimFF {
    my ($ff_file, $ELEMENTS, $alter, $oldFF) = @_;
    my (%PARMS, $which_var, $in_data, $type_counter, @vdws, $hb_counter, $pi, $M2Cnv, $parmType);
    my ($atom1, $atom2, @bonds, $atom3, @angles, $tmp1, @inversions, $currParm, $istip4p);
    my ($atom4, @torsions, $torsion_type, $inversion_type, $counter, $use_hb, $use_charge);
    my ($bond_counter, $angle_counter, $torsion_counter, @tmp, $inversion_counter, $bool, $tmp);
    my ($vdwType, $vdwDat, $i, $j, $dihdr_scale, $CNV, $vdw_counter, $crossType);
    
    $CNV = &loadCerius2LammpsConverter;
    $M2Cnv = &loadMpsim2CeriusConverter;
    $which_var = $bond_counter = $angle_counter = $torsion_counter = $counter = 0;
    $use_hb = $hb_counter = $crossType = $istip4p = 0;
    if (defined($oldFF)) {
	%PARMS = %{ $oldFF };
    }
    $istip4p = 1 if ($ff_file =~ /tip4/);
    $PARMS{"PARMS"}{"is_tip4p"} = 1 if ($istip4p);
    $PARMS{"PARMS"}{"cut_vdw"} = 10.0;
    $PARMS{"PARMS"}{"cut_coul"} = 10.0;
    $PARMS{"PARMS"}{"coul_accuracy"} = 0.0001;
    $PARMS{"PARMS"}{"hbond_distance"} = 2.5;
    $PARMS{"PARMS"}{"hbond_angle"} = 90;
    $PARMS{"PARMS"}{"scale_torsions"} = 0;
    $PARMS{"PARMS"}{"single_inversion"} = 1;
    $PARMS{"PARMS"}{"dielectric"} = 1;
    $PARMS{"PARMS"}{"mix_rule"} = "geometric";
    $PARMS{"PARMS"}{"EXO_CYCLIC"} = 1.0;
    $PARMS{"PARMS"}{"scale_cou_12"} =  $PARMS{"PARMS"}{"scale_cou_13"} = 0;
    $PARMS{"PARMS"}{"scale_vdw_12"} =  $PARMS{"PARMS"}{"scale_vdw_13"} = 0;
    $PARMS{"PARMS"}{"scale_cou_14"} = $PARMS{"PARMS"}{"scale_vdw_14"} = 1;
    $PARMS{"PARMS"}{"same_scale"} = 1;
    $PARMS{"PARMS"}{"lbond"} = $PARMS{"PARMS"}{"langle"} = 1;
    $PARMS{"PARMS"}{"ltorsion"} = $PARMS{"PARMS"}{"linversn"} = 1;
    $PARMS{"PARMS"}{"lnonbond"} = $PARMS{"PARMS"}{"lhbond"} = 1;
    $PARMS{"PARMS"}{"bndxang"} = $PARMS{"PARMS"}{"angxang"} = 1;
    $PARMS{"PARMS"}{"use_angle_bond_cosine_real"} = 1;
    $PARMS{"PARMS"}{"anganginv"} = 1;
    $PARMS{"PARMS"}{"bndbndto"} = $PARMS{"PARMS"}{"angangto"} = 1;
    $use_charge = 0;
    $pi = 3.141592654;
    $type_counter = 0;
    $type_counter = ($PARMS{"PARMS"}{"type_counter"} + 1) if (exists($PARMS{"PARMS"}{"type_counter"}));
    open FORCEFIELD, $ff_file or die "Cannot open force field file $ff_file: $!\n";
    while (<FORCEFIELD>) {
	chomp;
	$in_data = $_;
        if ($in_data =~ /^END|AUTOTYPE|PITWIST|ADDED H|LONE PAIRS|GASTEIGER/) {
	    $which_var = 0;
	} elsif ($in_data =~ /^\s*ETOR SCAL\s+(\d+\.\d+)/) {
	    $PARMS{"PARMS"}{"EXO_CYCLIC"} = $1;
	} elsif ($in_data =~ /^\s*DIELCTRIC\s+(\d+\.\d+)/) {
	    $PARMS{"PARMS"}{"dielectric"} = $1;
	} elsif ($in_data =~ /^\s*ALL INVER\s+T/) {
	    $PARMS{"PARMS"}{"single_inversion"} = 0;
 	} elsif ($in_data =~ /^\s*TORS SCAL\s+T/) {
	    $PARMS{"PARMS"}{"scale_torsions"} = 1;
	} elsif ($in_data =~ /^\s*RNB GEOMN\s+F/) {
	    $PARMS{"PARMS"}{"mix_rule"} = "arithmetic";
	} elsif ($in_data =~ /^\s*NBEXBND\s+(\w)/) {
	    if (lc($1) eq "t") {
		$bool = 0;
	    } else {
		$bool = 1;
	    }
	    $PARMS{"PARMS"}{"scale_cou_12"} = $PARMS{"PARMS"}{"scale_vdw_12"} = $bool;
        } elsif ($in_data =~ /^\s*NBEXANG\s+(\w)/) {
            if (lc($1) eq "t") {
                $bool = 0;
            } else {
                $bool = 1;
            }
            $PARMS{"PARMS"}{"scale_cou_13"} = $PARMS{"PARMS"}{"scale_vdw_13"} = $bool;
 	} elsif ($in_data =~ /^\s*SCAL NB14\s+(\d+\.\d+)/) {
            $PARMS{"PARMS"}{"scale_coulomb"} = $1;
	    $PARMS{"PARMS"}{"scale_coulomb_14"} = $PARMS{"PARMS"}{"scale_vdw_14"} = $1;
	} elsif ($in_data =~ /^\s*LHBOND\s+T/) {
	    $use_hb = 1;
	    $PARMS{"PARMS"}{"hbond_distance"} = 4.5;
	    $PARMS{"PARMS"}{"hbond_angle"} = 60;
	} elsif ($in_data =~ /\s*ANGX 2 K\s+F/) {
	    $PARMS{"PARMS"}{"use_angle_bond_cosine_real"} = 0;
	} elsif ($in_data =~ /^\s*(\w+)\s+(T|F)\s*$/) {
	    if (exists($PARMS{PARMS}{lc($1)})) {
		if (lc($2) eq "t") {
		    $bool = 1;
		} else {
		    $bool = 0;
		}
		$PARMS{PARMS}{lc($1)} = $bool;
	    }
	} elsif ($in_data =~ /^FFLABEL/) {
	    $which_var = 1;
	} elsif ($in_data =~ /^MPSIM_HB/ and $PARMS{PARMS}{lhbond}) {
	    $which_var = 2;	    
	} elsif ($in_data =~ /^VDW AT ITY/ and $PARMS{PARMS}{lnonbond}) {
	    $which_var = 3;
	} elsif ($in_data =~ /^BONDSTRTCH/ and $PARMS{PARMS}{lbond}) {
	    $which_var = 4;
	} elsif ($in_data =~ /^ANGLE/ and $PARMS{PARMS}{langle}) {
	    $which_var = 5;
	    $crossType = "";
	} elsif ($in_data =~ /^TORSION/ and $PARMS{PARMS}{ltorsion}) {
	    $which_var = 6;
	    $crossType = "";
	} elsif ($in_data =~ /^INVERSION/ and $PARMS{PARMS}{linversn}) {
	    $which_var = 7;
	    $crossType = "";
	} elsif ($in_data =~ /^NONBOND\-OFF/) {
	    $which_var = 8;
	} elsif ($in_data =~ /^ANGANGINV/) {
	    $which_var = 9;
	    $crossType = "AngleAngle";
	} elsif ($in_data =~ /^(\S+)\s+(\d+)(.{10})\s+(\-?\d+\.?\d*)\s+(\-?\d+)\s+(\-?\d+)\s+\d+\s+\d+\s+(\d+)/ and ($which_var == 1)) {
	    $type_counter += 1;
	    next if (! exists($ELEMENTS->{$2}));
	    $i = $ELEMENTS->{$2}{SYMBOL};
	    $tmp1 = $ELEMENTS->{$2}{MASS};
	    updateMass($tmp1, $3);
	    $currParm = \%{ $PARMS{"ATOMTYPES"} };
	    $currParm->{$1} = (
			       {
				   "TYPEID"     => $type_counter,
				   "ATOM"       => $i,
				   "ATOMNUM"    => $2,
				   "MASS"       => $tmp1,
				   "CHARGE"     => $4,
				   "HYBRID"     => $5,
				   "NUMBONDS"   => $6,
				   "LONEPAIRS"  => $7,
				   "LABEL"      => $1,
				   "USED"       => 0,
				   "USE_CHARGE" => $use_charge,
			       }
			       );
#	    print "ATOMTYPES: $1: $type_counter\n";
        } elsif ($in_data =~ /^\s*(\S+)\s*\-\s*(\S+)\s*\-\s*(\S+)\s+(\d+)\s+(\-?\d+.*)/i and $which_var == 2) { #hbonds
            next if (! exists($PARMS{"ATOMTYPES"}{$1}));
            next if (! exists($PARMS{"ATOMTYPES"}{$2}));
            next if (! exists($PARMS{"ATOMTYPES"}{$3}));
            $PARMS{PARMS}{HAS_HBONDS} = 1;
	    $atom1 = $1;
            $atom2 = $3;
            $vdwType = $4;
	    $tmp = $5;
	    @vdws = ($2);
	    $counter = 0;
            while ($tmp =~ /(\-?\d+\.?\d*)/g) {
		push @vdws, $1;
		$counter++;
		last if ($counter == 4);
	    }
            if ($vdwType == 2) { #morse
		($vdws[2],$vdws[3]) = ($vdws[3],$vdws[2]);
                $vdws[2] /= ($vdws[3] * 2); # changed from div to multi 07/28/2007
                #$vdws[0] *= 2;
                #$vdws[1] /= ($vdws[2]/2);
	    }
            if (! exists($PARMS{VDW}{$atom1}{$atom2})) {
                $vdw_counter = 1;
            } else {
                $vdw_counter = scalar(keys %{ $PARMS{VDW}{$atom1}{$atom2} }) + 1;
            }
	    if (($vdwType == 1 and scalar @vdws == 2) or ($vdwType == 2 and scalar @vdws == 3)) { 
		# append the power for the cosine - 4 - to the end of the data
		push @vdws, "4";
	    }
	    $vdwType = getCeriusTypes($4, "hbond", $M2Cnv);
	    $currParm = \%{ $PARMS{"VDW"}{$atom1}{$atom2} };
            $currParm->{$vdw_counter} = (
					 {
					     "TYPE"   => $vdwType,
					     "Lammps" => getLammpsOpts($vdwType,"hbond", $CNV),
					     "KEY"    => "$atom1 $atom2 ",
					     "VALS"   => [@vdws],
					     "ATOM"   => "$atom1 $atom2 ",
					     "USED"   => 0,
					     "IT"     => "hbond",
					 }
					 );
	    
	} elsif ($in_data =~ /^\s*(\S+)\s*\-?\s*(\S+)\s+(\S+)\s+(\S+)\s*(.*)/i and ($which_var == 3 or $which_var == 8)) {
	    $atom1 = $1;
	    if ($which_var == 3) {
		$atom2 = $1;
		$vdwType = $2;
		$vdwDat = "$3 $4 $5";
	    } else {
		$atom2 = $2;
		$vdwType = $3;
		$vdwDat = "$4 $5";
	    }
            if ($atom1 lt $atom2) {
                ($atom1, $atom2) = ($atom2, $atom1);
            }
   	    next if (! exists($PARMS{"ATOMTYPES"}{$atom1}));
	    next if (! exists($PARMS{"ATOMTYPES"}{$atom2}));
	    if ($vdwType != 6) {
		if ($vdwDat =~ /(.*)\# 1\-4 scaling: (.*)/) {
		    @vdws = split /\s+/, "$1 $2";
		} else {
		    @vdws = split /\s+/, $vdwDat;
		}
		@vdws = GetAbs(\@vdws);
		while ($#vdws > 3) {
		    pop @vdws;
	 	}
		#pop @vdws if ($vdwType != 33);
		if (! defined($alter) || $alter != 0) {
		    if ($vdwType == 3 || $vdwType == 33) { #morse or stretch morse
			($vdws[0], $vdws[1]) = ($vdws[1], $vdws[0]);
			($vdws[1], $vdws[2]) = ($vdws[2], $vdws[1]);
			$vdws[1] /= ($vdws[2] * 2); # changed from div to multi 07/28/2007
			$vdws[3] /= ($vdws[2] * 2) if (scalar(@vdws) > 3); #aplha2 for stretch morse
			#$vdws[0] *= 2;
			#$vdws[1] /= ($vdws[2]/2);
		    } elsif ($vdwType == 1) { #lj6_12
			($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
			$vdws[1] = $vdws[1] / (2**(1/6));
			while ($#vdws > 1) {
			    pop @vdws;
			}
#			$vdws[1] = $vdws[1]/2;
		    } elsif ($vdwType == 2) { #exponential-6
			($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
			@tmp = @vdws;
			@vdws = ();
			#$vdws[0] = 6 * $tmp[1] * exp($tmp[2]) / ($tmp[2] - 6);
			#$vdws[1] = $tmp[0]/$tmp[2];
			#$vdws[2] = $tmp[1] * $tmp[2] * $tmp[0]**6/($tmp[2] - 6);
			$vdws[0] = $tmp[0] * (6/($tmp[2]-6)) * exp($tmp[2]);
			$vdws[1] = $tmp[1]/$tmp[2];
			$vdws[2] = $tmp[1]**6 * $tmp[0] * ($tmp[2]/($tmp[2] - 6));
		    } else {
			next;
		    }
		}
		if (! exists($PARMS{VDW}{$atom1}{$atom2})) {
		    $vdw_counter = 1;
		} else {
		    $vdw_counter = scalar(keys %{ $PARMS{VDW}{$atom1}{$atom2} }) + 1;
		}
		$tmp1 = "vdw";
		$currParm = \%{ $PARMS{"VDW"}{$atom1}{$atom2} };
		$currParm->{$vdw_counter} = (
					     {
						 "TYPE"   => getCeriusTypes($vdwType, "vdw", $M2Cnv),
						 "Lammps" => getLammpsOpts(getCeriusTypes($vdwType, "vdw", $M2Cnv),"vdw", $CNV, $istip4p),
						 "KEY"    => "$atom1 $atom2 ",
						 "VALS"   => [@vdws],
						 "ATOM"   => "$atom1 $atom2 ",
						 "USED"   => 0,
						 "IT"     => "vdw",
						 "NUM"    => $vdwType,
					     }
					     );
#		print "VDW: $1: $atom\n"
	    }
	} elsif ($in_data =~ /^\s*(\S+)\s*\-\s*(\S+)\s+(\d+)\s+(.+)/ and ($which_var == 4)) { #bond
	    if ($3 >0 and $3 < 3) { #only allow harmonic and morse type bonds
		next if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
		next if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
		
		@bonds = split /\s+/, $4;
		next if (! @bonds);
#		GetAbs(\@bonds);
		if (! defined($alter) || $alter != 0) {
                    $bonds[0] = $bonds[0] / 2; #fix for the 1/2 factor in harmonic eqns between cerius and lammps
		    if ($3 == 2) {
			$bonds[0] *= 2;
			($bonds[0], $bonds[2]) = ($bonds[2], $bonds[0]); # kb is now bonds[0]
			($bonds[1], $bonds[2]) = ($bonds[2], $bonds[1]); # r0 is now bonds[2]
			$bonds[1] = sqrt($bonds[1]/(2 * $bonds[0]));
			while ($#bonds > 2) {
			    pop @bonds;
			}
		    } else {
			while ($#bonds > 1) {
			    pop @bonds;
			}
		    }
		    ($atom1, $atom2) = ($1, $2);
                    ($atom1, $atom2) = ($2, $1) if (exists($PARMS{BONDS}{$2}{$1}));
                    $bond_counter = 1;
                    if (exists($PARMS{BONDS}{$atom1}{$atom2})) {
                        $bond_counter = findDuplicate($PARMS{BONDS}{$atom1}{$atom2}, $3);
                        $bond_counter = scalar(keys %{ $PARMS{BONDS}{$atom1}{$atom2} }) + 1 if (! $bond_counter);
                    }

		    $currParm = \%{ $PARMS{"BONDS"}{$atom1}{$atom2} };
		    $parmType = getCeriusTypes($3, "bond", $M2Cnv);
		    $currParm->{$bond_counter} = (
						  {
						      "INDEX"    => $bond_counter,
						      "TYPE"     => $parmType,
						      "Lammps"   => getLammpsOpts($parmType,"bond", $CNV),
						      "VALS"     => [@bonds],
						      "USED"     => 0,
						      "KEY"      => "$1 $2 ",
						  }
						  );
		    
#		print "BOND $bond_counter: $key_code\n";
		}
	    }
	} elsif ($in_data =~ /^\s*(\S+)\s*\-\s*(\S+)\s*\-\s*(\S+)\s+(\d+)\s+(.+)/ and ($which_var == 5)) { #angle
	    if ($4 <= 31) { #no support yet for mm2/polyene/schleyer angles
		next if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
		next if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
		next if (! exists($PARMS{"ATOMTYPES"}{$3}) and lc($3) ne "x");
		$tmp1 = ();
		@angles = split /\s+/, $5;
		next if (! @angles);
#		GetAbs(\@angles);
		if ($crossType eq "" && (! defined($alter) || $alter != 0)) {
	            $angles[0] = $angles[0]/2; #same fix as bonds
		}
		if ($4 == 1) { #cosine harmonic
		    if (sin($angles[1] * $pi/180) > 0.001) {
			$angles[0] /= sin($angles[1] * $pi/180)**2;
		    }
		    while ($#angles > 1) {
			pop @angles;
		    }
		} elsif ($4 == 21) { #harmonic
                    while ($#angles > 1) {
                        pop @angles;
                    }
		} elsif ($4 == 4) { #cosine periodic
		    $tmp1 = $angles[7];
                    while ($#angles > 1) {
                        pop @angles;
                    }
		    push @angles, $tmp1;
		} elsif ($4 == 11) {
                    if (sin($angles[1] * $pi/180) > 0.001) {
                        $angles[0] /= sin($angles[1] * $pi/180)**2;
                    }
                    $tmp1 = ();
                    while ($#angles > 1) {
                        push @{ $tmp1 }, pop @angles;
                    }
                    @{ $tmp1 } = reverse @{ $tmp1 };
		} elsif ($4 == 31) {
		    $tmp1 = ();
                    while ($#angles > 1) {
                        push @{ $tmp1 }, pop @angles;
                    }
		    @{ $tmp1 } = reverse @{ $tmp1 };
		}
                ($atom1, $atom2, $atom3) = ($1, $2, $3);
                ($atom1, $atom2, $atom3) = ($3, $2, $1) if (exists($PARMS{ANGLES}{$3}{$2}{$1}));
                $angle_counter = 1;
                if (exists($PARMS{ANGLES}{$atom1}{$atom2}{$atom3})) {
                    $angle_counter = findDuplicate($PARMS{ANGLES}{$atom1}{$atom2}{$atom3}, $4);
                    $angle_counter = scalar(keys %{ $PARMS{ANGLES}{$atom1}{$atom2}{$atom3} }) + 1 if (! $angle_counter);
                }
		$currParm = \%{ $PARMS{"ANGLES"}{$atom1}{$atom2}{$atom3} };
		$parmType = getCeriusTypes($4, "angle", $M2Cnv);
		$currParm->{$angle_counter} =  (
						{
						    "INDEX"    => $angle_counter,
						    "TYPE"     => $parmType,
						    "Lammps"   => getLammpsOpts($parmType,"angle", $CNV),
						    "VALS"     => [@angles],
						    "USED"     => 0,
						    "KEY"      => "$1 $2 $3 ",
						    "CTYPE"    => $crossType,
						}
						);
		if ($tmp1) { # cross terms
		    $atom1 = getEquilVal($PARMS{BONDS}, $1, $2);
		    $atom2 = getEquilVal($PARMS{BONDS}, $2, $3);
		    $atom3 = $angles[1];
		    next if (! defined($atom1) or ! defined($atom2));
		    if ($PARMS{PARMS}{bndxang} and $tmp1->[2] != 0) { #stretch stretch
			@angles = ();
			@angles = ($atom1, $atom2, $tmp1->[2]); # stretch stretch force constant
			$angle_counter++;
			$parmType = getCeriusTypes("a$4", "angle", $M2Cnv);
                        $angle_counter = findDuplicate($currParm, $parmType);
                        $angle_counter = scalar(keys %{ $currParm }) + 1 if (! $angle_counter);
			$currParm->{$angle_counter} =  (
							{
							    "INDEX"    => $angle_counter,
							    "TYPE"     => $parmType,
							    "Lammps"   => getLammpsOpts($parmType,"angle", $CNV),
							    "VALS"     => [@angles],
							    "USED"     => 0,
							    "KEY"      => "$1 $2 $3 ",
							    "CTYPE"    => "BndBnd",
							}
							);
		    }
		    if ($PARMS{PARMS}{angxang}) {
			# angle stretch
                        @angles = ();
                        @angles = ($atom1, $atom2, $atom3, $tmp1->[0], $tmp1->[1]); 
			if ($PARMS{"PARMS"}{"use_angle_bond_cosine_real"} and $4 == 11) { #fix for bond angle cosine cross
			    $angles[3] = -$angles[3]*sin($angles[2] * $pi/180);
			    $angles[4] = -$angles[4]*sin($angles[2] * $pi/180);
			}			     
			$parmType = getCeriusTypes("b$4", "angle", $M2Cnv);
                        $angle_counter = findDuplicate($currParm, $parmType);
                        $angle_counter = scalar(keys %{ $currParm }) + 1 if (! $angle_counter);
                        $currParm->{$angle_counter} =  (
							{
							    "INDEX"    => $angle_counter,
							    "TYPE"     => $parmType,
							    "Lammps"   => getLammpsOpts($parmType,"angle", $CNV),
							    "VALS"     => [@angles],
							    "USED"     => 0,
							    "KEY"      => "$1 $2 $3 ",
							    "CTYPE"    => "BondAngle",
							}
							);
                    }
#		print "ANGLE $angle_counter: $key_code\n";
		}
	    }
	} elsif ($in_data =~ /^\s*(\S+)\s*\-\s*(\S+)\s*\-\s*(\S+)\s*\-\s*(\S+)\s+(\d+)\s*(.+)/ and ($which_var == 6)) { #torsions
	    $torsion_type = 1;
	    next if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
	    next if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
	    next if (! exists($PARMS{"ATOMTYPES"}{$3}) and lc($3) ne "x");
	    next if (! exists($PARMS{"ATOMTYPES"}{$4}) and lc($4) ne "x");
	    
	    $torsion_counter++;
            @torsions = split /\s+/, $6;
#		GetAbs(\@torsions);
	    
            if ($crossType eq "" && (! defined($alter) || $alter != 0)) {
                $torsions[0] = $torsions[0] / 2; #1/2 fix again
	    }
            if ($torsions[2] == 1) {
		$torsions[2] = 180;
            } else {
		$torsions[2] = 0;
            }
	    $tmp1 = ();
	    while ($#torsions > 2) {
		push @{ $tmp1 }, pop @torsions;
	    }
	    @{ $tmp1 } = reverse @{ $tmp1 } if ($tmp1);  
	    #$torsions[2] = 360 - $torsions[2];
	    $currParm = \%{ $PARMS{"TORSIONS"}{$1}{$2}{$3}{$4} };
	    $parmType = getCeriusTypes($torsion_type, "dihedral", $M2Cnv);
	    if (! exists($currParm->{$5})) {
		$currParm->{$5} =  (
				    {
					"TYPE"     => $parmType,
					"Lammps"   => getLammpsOpts($parmType,"dihedral", $CNV),
					"INDEX"    => $torsion_counter,
					"VALS"     => [@torsions],
					"USED"     => 0,
					"KEY"      => "$1 $2 $3 $4 ",
					"NUM"      => 1,
					"PER"      => ($#torsions + 1),
					"1_4scale" => $PARMS{"PARMS"}{"scale_vdw_14"},
					"CTYPE"    => $crossType,
					"do_scale" => $PARMS{"PARMS"}{"scale_torsions"},
				    }
				    );
	    } else {
	        push @{ $currParm->{$5}{"VALS"} }, @torsions;
		$currParm->{$5}{NUM}++;
#		print "TORSION $torsion_counter: $key_code\n";
	    }
	    # cross terms
	    if ($PARMS{PARMS}{angangto} and $tmp1 and $tmp1->[0] != 0) { # angle angle torsion (TORSION_BEND_BEND)
                $atom1 = getEquilVal($PARMS{ANGLES}, $1, $2, $3);
                $atom2 = getEquilVal($PARMS{ANGLES}, $2, $3, $4);
		next if (! defined($atom1) or ! defined($atom2));
		@torsions = ($atom1, $atom2, $tmp1->[0]);
		$parmType = getCeriusTypes("a$torsion_type", "dihedral", $M2Cnv);
                $currParm->{$torsion_counter} =  (
						  {
						      "TYPE"     => $parmType,
						      "Lammps"   => getLammpsOpts($parmType,"dihedral", $CNV),
						      "INDEX"    => $torsion_counter,
						      "VALS"     => [@torsions],
						      "USED"     => 0,
						      "KEY"      => "$1 $2 $3 $4 ",
						      "NUM"      => 1,
						      "PER"      => ($#torsions + 1),
						      "1_4scale" => $PARMS{"PARMS"}{"scale_vdw_14"},
						      "CTYPE"    => $crossType,
						      "do_scale" => $PARMS{"PARMS"}{"scale_torsions"},
						  }
						  );
		$torsion_counter++;
	    }
	} elsif ($in_data =~ /^\s*(\S+)\s*\-\s*(\S+)\s*\-\s*(\S+)\s*\-\s*(\S+)\s+(\d+)\s+(.+)/ and ($which_var == 7)) { #inversion
	    $inversion_type = $5;
	    if ($inversion_type < 4) {
		next if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
		next if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
		next if (! exists($PARMS{"ATOMTYPES"}{$3}) and lc($3) ne "x");
		next if (! exists($PARMS{"ATOMTYPES"}{$4}) and lc($4) ne "x");
		
		$inversion_counter++;
                @inversions = split /\s+/, $6;
		next if (! @inversions);
#		GetAbs(\@torsions);
 		while ($#inversions > 1) { # there is an extra parameter that i don't know the use
		    pop @inversions;
		}
		if ($inversion_type == 3) {
		    $atom1 = $2;
		    $atom2 = $3;
		    $atom3 = $1;
		    $atom4 = $4;
		} else {
		    ($atom1, $atom2, $atom3, $atom4)  = ($1, $2, $3, $4);
		}

                if ($crossType eq "" && (! defined($alter) || $alter != 0)) {
                    if ($inversion_type == 3) {
			if ($inversions[1] == 180) {
			    $inversions[1] = -1;
			} else {
			    $inversions[1] = 1;
			}
                    } 		    $inversions[0] = $inversions[0] / 2 if ($inversion_type !=  2); #1/2 fix again
		}
		$currParm = \%{ $PARMS{"INVERSIONS"}{$atom1}{$atom2}{$atom3}{$atom4} };
		$parmType = getCeriusTypes($inversion_type, "inversion", $M2Cnv);
		$currParm->{$inversion_counter} = (
						   {
						       "TYPE"     => $parmType,
						       "Lammps"   => getLammpsOpts($parmType,"inversion", $CNV),
						       "INDEX"    => $inversion_counter,
						       "VALS"     => [@inversions],
						       "USED"     => 0,
						       "KEY"      => "$1 $2 $3 $4 ",
						       "CTYPE"    => $crossType,
						   }
						   );
#		print "TORSION $torsion_counter: $key_code\n";
	    }
	} elsif ($in_data =~ /^\s*(\S+)\s*\-\s*(\S+)\s*\-\s*(\S+)\s*\-\s*(\S+)\s+(\d+)\s+(.+)/ and ($which_var == 9)) {	# angle angle inversion
            $inversion_type = 1;
            next if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
            next if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
            next if (! exists($PARMS{"ATOMTYPES"}{$3}) and lc($3) ne "x");
            next if (! exists($PARMS{"ATOMTYPES"}{$4}) and lc($4) ne "x");
	    
            @{ $j } = split /\s+/, $6;
            next if (! @{ $j});
	    
            @{ $tmp1 } = ([$1, $2, $3, $4],[$1, $3, $2, $4],[$1, $4, $3, $3]);
	    for $i (0 .. 2) {
		$inversion_counter++;
		$atom1 = getEquilVal($PARMS{ANGLES}, $tmp1->[$i][1], $tmp1->[$i][0], $tmp1->[$i][2]);
		$atom2 = getEquilVal($PARMS{ANGLES}, $tmp1->[$i][2], $tmp1->[$i][0], $tmp1->[$i][3]);
		next if (! defined($atom1) or ! defined($atom2));
		@inversions = ($atom1, $atom2, $j->[$i]);
		$currParm = \%{ $PARMS{"INVERSIONS"}{$1}{$2}{$3}{$4} };
		$parmType = getCeriusTypes("a$inversion_type", "inversion", $M2Cnv);
		$currParm->{$inversion_counter} = (
						   {
						       "TYPE"     => $parmType,
						       "Lammps"   => getLammpsOpts($parmType,"inversion_cross", $CNV),
						       "INDEX"    => $inversion_counter,
						       "VALS"     => [@inversions],
						       "USED"     => 0,
						       "KEY"      => "$1 $2 $3 $4 ",
						       "CTYPE"    => $crossType,
						   }
						   );
	    }
	}
    }
    close FORCEFIELD;
    
    die "ERROR: Force Field file $ff_file is invalid!\n"
	if (! %PARMS);
    
    $PARMS{"PARMS"}{"type_counter"} = $type_counter;
    return (\%PARMS);
}

sub saveMPSimFF {
    my ($ffData, $save_name, $ELEMENTS) = @_;
    my ($c, $shft_angle, $torsions, $counter, $element, $vdws, $i);
    my ($atom1, $atom2, $atom3, $atom4, @ATMS, $hk, $inversions, @TMP);
    my ($pi) = atan2(1,1) *4;

    open OUTFILE, "> $save_name" or die "Cannot create file $save_name: $!\n";
    print OUTFILE "ATOMTYPES\n";
    for $atom1 (keys %{ $ffData->{"atoms"} }) {
	next
	    if ($atom1 eq "?");
	$c = \%{ $ffData->{"atoms"}{$atom1}{"VALS"}[0] };
	next if (! exists($c->{r}));
        $element = $ELEMENTS->{ $c->{"element"} }{"SYMBOL"};
	printf OUTFILE " %-11s%-3s%12.5f%8.4f%4s%4s%4s",
	$atom1,$element,$c->{"mass"},0,$c->{"hybrid"},0,0;
        $vdws .= sprintf(" %-11s LJ_6_12%14.4f%8.4f",
	$atom1, ($c->{"r"}  * 2), $c->{"e"});
	$c->{r} = 0.00001 if ($c->{r} == 0);
	#if (! exists($c->{"1_4"})) {
	    #$c->{"1_4"}{"r"} = $c->{"r"};
	    #$c->{"1_4"}{"e"} = $c->{"e"};
	#}
        if (exists($c->{"1_4"})) {
	    $vdws .= sprintf(" # 1-4 scaling: %8.5f%8.5f\n",($c->{"1_4"}{"r"}*2),abs($c->{"1_4"}{"e"}));
	} else {
	    $vdws .= "\n";
        }
	if (exists($c->{"name"})) {
	    printf OUTFILE " # $c->{name}";
	}
	printf OUTFILE "\n";
    }
    print OUTFILE "END\n#\nDIAGONAL_VDW\n$vdws";
    print OUTFILE "END\n#\nATOM_TYPING_RULES\nEND\n\#\nOFF_DIAGONAL_VDW\nEND\n#\nBOND_STRETCH\n";
    for $atom1 (keys %{ $ffData->{"bonds"} }) {
        for $atom2 (keys %{ $ffData->{"bonds"}{$atom1} }) {
            $c = \%{ $ffData->{"bonds"}{$atom1}{$atom2}{"VALS"}[0] };
            printf OUTFILE " %-9s%-11s HARMONIC%13.4f%10.4f\n",
            $atom1, $atom2, ($c->{"kb"}  * 2), $c->{"r0"};
        }
    }
    print OUTFILE "END\n#\nANGLE_BEND\n";
    for $atom1 (keys %{ $ffData->{"angles"} }) {
        for $atom2 (keys %{ $ffData->{"angles"}{$atom1} }) {
            for $atom3 (keys %{ $ffData->{"angles"}{$atom1}{$atom2} }) {
                $c = \%{ $ffData->{"angles"}{$atom1}{$atom2}{$atom3}{"VALS"}[0] };
                $c->{"t0"} = int($c->{"t0"} * 180/$pi);
                printf OUTFILE " %-9s%-9s%-11s THETA_HARM%10.4f%10.4f\n",
                $atom1, $atom2, $atom3, ($c->{"kt"} * 2), $c->{"t0"};
            }
        }
    }
    print OUTFILE "END\n#\nUREY_BRADLEY\n";
    for $atom1 (keys %{ $ffData->{"urey_bradley"} }) {
        for $atom2 (keys %{ $ffData->{"urey_bradley"}{$atom1} }) {
            for $atom3 (keys %{ $ffData->{"urey_bradley"}{$atom1}{$atom2} }) {
                $c = \%{ $ffData->{"urey_bradley"}{$atom1}{$atom2}{$atom3}{"VALS"}[0] };
                printf OUTFILE " %-9s%-9s%-11s HARMONIC%13.4f%10.4f\n",
                $atom1, $atom2, $atom3, ($c->{"ku"} * 2), $c->{"su"};
            }
        }
    }
    print OUTFILE "END\n#\nTORSIONS\n";
    for $atom1 (keys %{ $ffData->{"torsions"} }) {
        for $atom2 (keys %{ $ffData->{"torsions"}{$atom1} }) {
            for $atom3 (keys %{ $ffData->{"torsions"}{$atom1}{$atom2} }) {
                for $atom4 (keys %{ $ffData->{"torsions"}{$atom1}{$atom2}{$atom3} }) {
                    @TMP = @ATMS = ();
                    $torsions = \@{ $ffData->{"torsions"}{$atom1}{$atom2}{$atom3}{$atom4}{"VALS"} };
                    @TMP = ( $atom1, $atom2, $atom3, $atom4 );
                    $i = $ffData->{"torsions"}{$atom1}{$atom2}{$atom3}{$atom4}{"counter"};
                    while ($ffData->{"torsionOrders"}{$i} =~ /(\d)/g) {
                        push @ATMS, $TMP[$1];
                    }
                    printf OUTFILE " %-9s%-9s%-9s%-11s SHFT_DIHDR", 
                    $ATMS[0], $ATMS[1], $ATMS[2], $ATMS[3];
                    for $counter (0 .. $#{ $torsions }) {
                        $c = \%{ $torsions->[$counter] };
                        $c->{"p0"} = (int($c->{"p0"} * 180/$pi));# % 360);
                        if ($counter == 0) {
                            printf OUTFILE "%11.4f%10.4f%10.4f\n", ($c->{"kp"} * 2), $c->{"n"}, $c->{"p0"};
                        } else {
                            printf OUTFILE "%50s%11.4f%10.4f%10.4f\n", "", ($c->{"kp"} * 2), $c->{"n"}, $c->{"p0"};
                        }
                    }
                }
            }
        }
    }
    
    print OUTFILE "END\n#\nINVERSIONS\n";
    for $atom1 (keys %{ $ffData->{"inversions"} }) {
        for $atom2 (keys %{ $ffData->{"inversions"}{$atom1} }) {
            for $atom3 (keys %{ $ffData->{"inversions"}{$atom1}{$atom2} }) {
                for $atom4 (keys %{ $ffData->{"inversions"}{$atom1}{$atom2}{$atom3} }) {
                    @TMP = @ATMS = ();
                    $torsions = \@{ $ffData->{"inversions"}{$atom1}{$atom2}{$atom3}{$atom4}{"VALS"} };
                    @TMP = ( $atom1, $atom2, $atom3, $atom4 );
                    $i = $ffData->{"inversions"}{$atom1}{$atom2}{$atom3}{$atom4}{"counter"};
                    while ($ffData->{"inversionOrders"}{$i} =~ /(\d)/g) {
                        push @ATMS, $TMP[$1];
                    }
                    $c = \%{ $torsions->[0] };
		    if (! defined($c->{type}) || $c->{type} eq "IT_JIKL") {
			printf OUTFILE " %-9s%-9s%-9s%-11s IT_JIKL%14.4f%10.4f%10.4f\n",
			     $ATMS[0], $ATMS[1], $ATMS[2], $ATMS[3],
			     ($c->{"kp"} * 2), $c->{"p0"},$c->{"n"};
		    } else {
			printf OUTFILE " %-9s%-9s%-9s%-11s IT_IJKL%14.4f%10.4f\n",
			     $ATMS[0], $ATMS[1], $ATMS[2], $ATMS[3],
			     ($c->{"kp"} * 2), $c->{"p0"};
		    }
                }  
            }
        }
    }
    print OUTFILE "END\n\#\nCOULOMBIC\n X        X           CONST-EPS\nEND\n";
    close OUTFILE;
}

sub GetAbs {
    my ($inVals) = @_;
    my ($counter, @ret);
                                                                                                                              
    for $counter (@{ $inVals }) {
        if ($counter =~ /(\-\d+\.?\d*)e(\-\d+\.?\d*)/i) {
            push @ret, $1 * 10 ** $2;
        } elsif ($counter =~ /\d+\.?\d*/) {
            push @ret, $counter;
        }
    }
    return (@ret);
}


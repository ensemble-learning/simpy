package Packages::CERIUS2;

require Exporter;
use strict;

our(@ISA, @EXPORT, $VERSION);

@ISA = qw(Exporter);
@EXPORT = qw(parseCerius2FF saveCeriusFF LoadFFs ReadFFs parseMPSimFF);
$VERSION = "1.00";

sub loadConverter {
    my ($inFile) = "/net/hulk/ul3/tpascal/scripts/dat/ceriusParms2Lammps.perldata";
    my (%SHASH);
    scalar eval `cat $inFile` or die "Cannot recreate data in file $inFile: $! $@\n";
    return \%SHASH;
}

sub getLammpsOpts {
    my ($parmName, $parmType, $SHASH) = @_;
    my (%RET);
    my ($searchStr) = lc($parmType) . "_" . lc($parmName);
    
    if (exists($SHASH->{$searchStr})) {
	%RET = %{ $SHASH->{$searchStr} };
    } else {
	$RET{name} = lc($parmName);
	$RET{opts} = "";
    }
    return \%RET;
}
    
sub parseCerius2FF {
    my ($ff_file, $alter, $oldFF) = @_;
    my (%PARMS, $which_var, $in_data, $type_counter, @vdws, $hb_counter, $pi);
    my ($type_id, $type_id2, @bonds, $type_id3, @angles, $tmp1, @inversions);
    my ($type_id4, @torsions, $torsion_type, $inversion_type, $counter, $use_hb);
    my ($bond_counter, $angle_counter, $torsion_counter, @tmp, $inversion_counter, $bool);
    my ($atom1, $atom2, $vdwType, $vdwDat, $i, $dihdr_scale, $CNV, $vdw_counter, $crossType);

    $CNV = &loadConverter;
    $which_var = $bond_counter = $angle_counter = $torsion_counter = $counter = 0;
    $use_hb = $hb_counter = $crossType = 0;
    if (defined($oldFF)) {
	%PARMS = %{ $oldFF };
    }
    $PARMS{"PARMS"}{"cut_vdw"} = 14.0;
    $PARMS{"PARMS"}{"cut_coul"} = 15.0;
    $PARMS{"PARMS"}{"coul_accuracy"} = 0.0001;
    $PARMS{"PARMS"}{"hbond_distance"} = 2.5;
    $PARMS{"PARMS"}{"hbond_angle"} = 90;
    $PARMS{"PARMS"}{"scale_torsions"} = 0;
    $pi = 3.141592654;

    open FORCEFIELD, $ff_file or die "Cannot open force field file $ff_file: $!\n";
    while (<FORCEFIELD>) {
	chomp;
	$in_data = $_;
        if ($in_data =~ /^END/) {
	    $which_var = 0;
 	} elsif ($in_data =~ /^\s*SCALE_TORSIONS_ABOUT_COMMON_BOND\s+T/) {
	    $PARMS{"PARMS"}{"scale_torsions"} = 1;
	} elsif ($in_data =~ /^\s+COU_DIRECT_CUT-OFF\s+(\d+\.\d+)/) {
	    $PARMS{"PARMS"}{"cut_coul"} = $1;
	} elsif ($in_data =~ /^\s+VDW_SPLINE_OFF\s+(\d+\.\d+)/) {
	    $PARMS{"PARMS"}{"cut_vdw"} = $1;
	} elsif ($in_data =~ /^\s+VDW_COMBINATION_RULE\s+(\w+)/) {
	    $PARMS{"PARMS"}{"mix_rule"} = lc($1);
	} elsif ($in_data =~ /^\s+EWALD_SUM_COU_ACCURACY\s+(\d+\.\d+)/) {
	    $PARMS{"PARMS"}{"coul_accuracy"} = $1;
	} elsif ($in_data =~ /^\s+(COU|VDW)_EXCLUDE_(\d)\-(\d)\s+(\w)/) {
	    if (lc($4) eq "t") {
		$bool = 0;
	    } else {
		$bool = 1;
	    }
	    $PARMS{"PARMS"}{"scale_" . lc($1) . "_" . $2 . $3} = $bool;
 	} elsif ($in_data =~ /^\s*COU_1-4_SCALE_FACTOR\s+(\d+\.\d+)/) {
            $PARMS{"PARMS"}{"scale_coulomb"} = $1;
	    $PARMS{"PARMS"}{"scale_coulomb_14"} = $1 if ($PARMS{"PARMS"}{"scale_coulomb_14"})
 	} elsif ($in_data =~ /^\s*VDW_1-4_SCALE_FACTOR\s+(\d+\.\d+)/) {
            $PARMS{"PARMS"}{"scale_vdw"} = $1;
	    $PARMS{"PARMS"}{"scale_vdw_14"} = $1 if ($PARMS{"PARMS"}{"scale_vdw_14"});
            if($PARMS{"PARMS"}{"scale_coulomb"} ne $PARMS{"PARMS"}{"scale_vdw"}) { 
		$PARMS{"PARMS"}{"same_scale"} = 0;
	    } else {
		$PARMS{"PARMS"}{"same_scale"} = 1;
	    }
	} elsif ($in_data =~ /^ HYDROGEN_BONDS\s+T/) {
	    $use_hb = 1;
	    $PARMS{"PARMS"}{"hbond_distance"} = 4.5;
	    $PARMS{"PARMS"}{"hbond_angle"} = 60;
	} elsif ($use_hb and $in_data =~ /^ H-BOND_LIST_DISTANCE_OFF\s+(\d+\.\d+)/) {
	    $PARMS{"PARMS"}{"hbond_distance"} = $1;
	} elsif ($use_hb and $in_data =~ /^ H-BOND_LIST_ANGLE_OFF\s+(\d+\.\d+)/) {
	    $PARMS{"PARMS"}{"hbond_angle"} = $1;
	} elsif ($in_data =~ /^ATOMTYPES/) {
	    $which_var = 1;
	} elsif ($in_data =~ /HYDROGEN_BONDS/) {
		$which_var = 2;	    
	} elsif ($in_data =~ /^DIAGONAL_VDW/) {
		$which_var = 3;
	} elsif ($in_data =~ /^BOND_STRETCH/) {
	    $which_var = 4;
	} elsif ($in_data =~ /^ANGLE_BEND/) {
	    $which_var = 5;
	    $crossType = "";
	} elsif ($in_data =~ /^TORSIONS/) {
	    $which_var = 6;
	    $crossType = "";
	} elsif ($in_data =~ /^INVERSIONS/) {
	    $which_var = 7;
	    $crossType = "";
	} elsif ($in_data =~ /^STRETCH_STRETCH/) {
	    $which_var = 5;
	    $crossType = "BondBond";
	} elsif ($in_data =~ /^STRETCH_BEND_STRETCH/) {
	    $which_var = 5;
	    $crossType = "BondAngle";
	} elsif ($in_data =~ /^OFF_DIAGONAL_VDW/) {
	    $which_var = 8;
	} elsif ($in_data =~ /^UREY_BRADLEY/) {
	    $which_var = 9;
	    $crossType = "";
	} elsif ($in_data =~ /^\s*(\S+)\s+(\w+)\s+(\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\d+)\s+(\d+)\s+(\d+)/ and ($which_var == 1)) {
	    $type_counter += 1;
	    $PARMS{"ATOMTYPES"}{$1} = (
				       {
					   "TYPEID"    => $type_counter,
					   "ATOM"      => $2,
					   "MASS"      => $3,
					   "CHARGE"    => $4,
					   "NUMBONDS"  => $5,
					   "LONEPAIRS" => $6,
					   "OTHER"     => $7,
					   "LABEL"     => $1,
					   "USED"      => 0,
				       }
				       );
#	    print "ATOMTYPES: $1: $type_counter\n";
        } elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+.*)/i and $which_var == 2) { #hbonds
            next if (! exists($PARMS{"ATOMTYPES"}{$1}));
            next if (! exists($PARMS{"ATOMTYPES"}{$2}));
            next if (! exists($PARMS{"ATOMTYPES"}{$3}));
	    $atom1 = $1;
            $atom2 = $3;
	    if ($atom1 lt $atom2) {
		($atom1, $atom2) = ($atom2, $atom1);
	    }
            $vdwType = $4;
            @vdws = split /\s+/, "$2 $5";
            if (! exists($PARMS{VDW}{$atom1}{$atom2})) {
                $vdw_counter = 1;
            } else {
                $vdw_counter = scalar(keys %{ $PARMS{VDW}{$atom1}{$atom2} }) + 1;
            }
            $PARMS{"VDW"}{$atom1}{$atom2}{$vdw_counter} = (
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

	} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*(.*)/i and ($which_var == 3 or $which_var == 8)) {
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
	    if (uc($vdwType) ne "IGNORE") {
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$1}));

		$type_id = $PARMS{"ATOMTYPES"}{$atom1}{"TYPEID"};
		if ($vdwDat =~ /(.*)\# 1\-4 scaling: (.*)/) {
		    @vdws = split /\s+/, "$1 $2";
		} else {
		    @vdws = split /\s+/, $vdwDat;
		}
		@vdws = GetAbs(\@vdws);

		if (! defined($alter) || $alter != 0) {
		    if (uc($vdwType) eq "VDW_MORSE") {
			($vdws[0], $vdws[1]) = ($vdws[1], $vdws[0]);
			($vdws[1], $vdws[2]) = ($vdws[2], $vdws[1]);
			$vdws[1] /= ($vdws[2] * 2); # changed from div to multi 07/28/2007
			#$vdws[0] *= 2;
			#$vdws[1] /= ($vdws[2]/2);
		    } elsif (uc($vdwType) eq "LJ_6_12") {
			($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
			$vdws[1] = $vdws[1] / (2**(1/6));
			if ($#vdws > 1) {
			    ($vdws[2],$vdws[3]) = ($vdws[3],$vdws[2]);
			    $vdws[3] = $vdws[3] /(2**(1/6));
			}
#			$vdws[1] = $vdws[1]/2;
		    } elsif (uc($vdwType) eq "EXPO_6") {
			($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
			@tmp = @vdws;
			@vdws = ();
			#$vdws[0] = 6 * $tmp[1] * exp($tmp[2]) / ($tmp[2] - 6);
			#$vdws[1] = $tmp[0]/$tmp[2];
			#$vdws[2] = $tmp[1] * $tmp[2] * $tmp[0]**6/($tmp[2] - 6);
			$vdws[0] = $tmp[0] * (6/($tmp[2]-6)) * exp($tmp[2]);
			$vdws[1] = $tmp[1]/$tmp[2];
			$vdws[2] = $tmp[1]**6 * $tmp[0] * ($tmp[2]/($tmp[2] - 6));
		    } 
		}
		if (! exists($PARMS{VDW}{$atom1}{$atom2})) {
		    $vdw_counter = 1;
		} else {
		    $vdw_counter = scalar(keys %{ $PARMS{VDW}{$atom1}{$atom2} }) + 1;
		}
		$tmp1 = "vdw";
		$PARMS{"VDW"}{$atom1}{$atom2}{$vdw_counter} = (
						 {
						     "TYPE"   => $vdwType,
						     "Lammps" => getLammpsOpts($vdwType,"vdw", $CNV),
						     "KEY"    => "$atom1 $atom2 ",
						     "VALS"   => [@vdws],
						     "ATOM"   => "$atom1 $atom2 ",
						     "USED"   => 0,
						     "IT"     => "vdw",
						 }
						 );
#		print "VDW: $1: $type_id\n"
	    }
	} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\w+)\s+(.+)/ and ($which_var == 4)) {
	    if (uc($3) ne "IGNORE") {
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
		
		$bond_counter++;
		@bonds = split /\s+/, $4;
		next if (! @bonds);
#		GetAbs(\@bonds);
		if (! defined($alter) || $alter != 0) {
                    $bonds[0] = $bonds[0] / 2; #fix for the 1/2 factor in harmonic eqns between cerius and lammps
		    if (uc($3) eq "MORSE") {
			$bonds[2] = sqrt($bonds[0]/(2 * $bonds[2]));
			($bonds[1], $bonds[2]) = ($bonds[2], $bonds[1]);
			pop @bonds;
		    }
		}
		$PARMS{"BONDS"}{$1}{$2} = (
					   {
					       "INDEX"    => $bond_counter,
					       "TYPE"     => $3,
					       "Lammps"   => getLammpsOpts($3,"bond"),
					       "VALS"     => [@bonds],
					       "USED"     => 0,
					       "KEY"      => "$1 $2 ",
					   }
					   );

#		print "BOND $bond_counter: $key_code\n";
	    }
	} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)/ and ($which_var == 5)) {
	    if (uc($4) ne "IGNORE") {
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$3}) and lc($3) ne "x");
		
		$angle_counter++;
		@angles = split /\s+/, $5;
		next if (! @angles);
#		GetAbs(\@angles);
		if ($crossType eq "" && (! defined($alter) || $alter != 0)) {
	            $angles[0] = $angles[0]/2; #same fix as bonds
		}
		if ($4 eq "COS_HARMON") {
		    if (sin($angles[1] * $pi/180) > 0.001) {
			$angles[0] /= sin($angles[1] * $pi/180)**2;
		    }
		}
		$PARMS{"ANGLES"}{$1}{$2}{$3} = (
						{
						    "INDEX"    => $angle_counter,
						    "TYPE"     => $4,
						    "Lammps"   => getLammpsOpts($4,"angle", $CNV),
						    "VALS"     => [@angles],
						    "USED"     => 0,
						    "KEY"      => "$1 $2 $3 ",
						    "CTYPE"    => $crossType,
						}
						);
#		print "ANGLE $angle_counter: $key_code\n";
	    }
        } elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)/ and ($which_var == 9)) {
            if (uc($4) ne "IGNORE") {
                next
                    if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
                next
                    if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
                next
                    if (! exists($PARMS{"ATOMTYPES"}{$3}) and lc($3) ne "x");

                $angle_counter++;
                @angles = split /\s+/, $5;
#               GetAbs(\@angles);
                if ($crossType eq "" && (! defined($alter) || $alter != 0)) {
                    $angles[0] = $angles[0]/2; #same fix as bonds
                }
		if (! exists($PARMS{"ANGLES"}{$1}{$2}{$3})) {
	            $PARMS{"ANGLES"}{$1}{$2}{$3} = (
            					   {
                                                    "INDEX"    => $angle_counter,
                                                    "TYPE"     => $4,
						    "Lammps"   => getLammpsOpts($4,"angle"),
                                                    "VALS"     => [@angles],
                                                    "USED"     => 0,
                                                    "KEY"      => "$1 $2 $3 ",
						    "CTYPE"    => $crossType,
                                                   }
                                                   );
		} else {
		    $PARMS{"ANGLES"}{$1}{$2}{$3}{TYPE} = "CHARMM";
		    $PARMS{"ANGLES"}{$1}{$2}{$3}{Lammps} = getLammpsOpts("CHARMM","angle");
		    for $i (@angles) {
			push @{ $PARMS{"ANGLES"}{$1}{$2}{$3}{VALS} }, $i;
		    }
		}
#               print "ANGLE $angle_counter: $key_code\n";
            }
	} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)/ and ($which_var == 6)) {
	    $torsion_type = $5;
	    if (uc($torsion_type) ne "IGNORE") {
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$3}) and lc($3) ne "x");
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$4}) and lc($4) ne "x");

		$torsion_counter++;
                @torsions = split /\s+/, $6;
#		GetAbs(\@torsions);
 
                if ($crossType eq "" && (! defined($alter) || $alter != 0)) {
                    if (lc($torsion_type) eq "dihedral") {
			if ($torsions[2] == 1) {
			    $torsions[2] = 180;
		    	} else {
			    $torsions[2] = 0;
		    	}
                    }
                    $torsions[0] = $torsions[0] / 2; #1/2 fix again
		}
		#$torsions[2] = 360 - $torsions[2];
		$PARMS{"TORSIONS"}{$1}{$2}{$3}{$4} = (
						      {
							  "TYPE"     => $torsion_type,
							  "Lammps"   => getLammpsOpts($torsion_type,"dihedral", $CNV),
							  "INDEX"    => $torsion_counter,
							  "VALS"     => [@torsions],
							  "USED"     => 0,
							  "KEY"      => "$1 $2 $3 $4 ",
							  "NUM"      => 1,
							  "PER"      => ($#torsions + 1),
							  "1_4scale" => $PARMS{"PARMS"}{"scale_vdw_14"},
							  "CTYPE"    => $crossType,
						      }
						      );
		($type_id, $type_id2, $type_id3, $type_id4)  = ($1, $2, $3, $4);
#		print "TORSION $torsion_counter: $key_code\n";
	    }
	} elsif ($in_data =~ /^\s+(\-?\d+\.\d+)\s+(.+)/ and $which_var == 6) { #multiple torsions
	    @torsions = ();
	    @torsions = split /\s+/, $2;
#	    GetAbs(\@torsions);
            if ($crossType eq "" && (! defined($alter) || $alter != 0)) {
                if (lc($torsion_type) eq "dihedral") {
		    if ($torsions[1] == 1) {
			$torsions[1] = 180;
		    } else {
			$torsions[1] = 0;
		    }
            	}

                $tmp1 = $1 / 2; #1/2 fix again
	    } else {
		$tmp1 = $1;
	    }

	    push @{ $PARMS{"TORSIONS"}{$type_id}{$type_id2}{$type_id3}{$type_id4}{"VALS"} }, $tmp1;
	    push @{ $PARMS{"TORSIONS"}{$type_id}{$type_id2}{$type_id3}{$type_id4}{"VALS"} }, @torsions;
	    $PARMS{"TORSIONS"}{$type_id}{$type_id2}{$type_id3}{$type_id4}{NUM}++;
	    $torsion_counter++;
#	    print "FOUND MULTIPLE TORSIONS FOR $key_code\n";
	} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)/ and ($which_var == 7)) {
	    $inversion_type = $5;
	    if (uc($inversion_type) ne "IGNORE") {
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$3}) and lc($3) ne "x");
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$4}) and lc($4) ne "x");

		$inversion_counter++;
                @inversions = split /\s+/, $6;
		next if (! @inversions);
#		GetAbs(\@torsions);
 
                if ($crossType eq "" && (! defined($alter) || $alter != 0)) {
                    if (lc($inversion_type) eq "it_jikl") {
			$type_id = $2;
			$type_id2 = $3;
			$type_id3 = $1;
			$type_id4 = $4;
			if ($inversions[1] == 180) {
			    $inversions[1] = -1;
			} else {
			    $inversions[1] = 1;
			}
                    } else {
			($type_id, $type_id2, $type_id3, $type_id4)  = ($1, $2, $3, $4);
	   	    }
		    $inversions[0] = $inversions[0] / 2; #1/2 fix again
		}
		$PARMS{"INVERSIONS"}{$1}{$2}{$3}{$4} = (
							{
							    "TYPE"     => $inversion_type,
							    "Lammps"   => getLammpsOpts($inversion_type,"inversion", $CNV),
							    "INDEX"    => $inversion_counter,
							    "VALS"     => [@inversions],
							    "USED"     => 0,
							    "KEY"      => "$1 $2 $3 $4 ",
							    "CTYPE"    => $crossType,
							}
							);
#		print "TORSION $torsion_counter: $key_code\n";
	    }
	}	    
	
    }
    close FORCEFIELD;

    die "ERROR: Force Field file $ff_file is invalid!\n"
	if (! %PARMS);
    
    if ($PARMS{"PARMS"}{"cut_vdw"} > $PARMS{"PARMS"}{"cut_coul"}) {
	$PARMS{"PARMS"}{"cut_vdw"} = $PARMS{"PARMS"}{"cut_coul"};
    }

    return (\%PARMS);
}

sub saveCeriusFF {
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
	if ($counter =~ /(.+)e(.+)/i) {
	    push @ret, $1 * 10 ** $2;
	} elsif ($counter =~ /\d+\.?\d*/) {
	    push @ret, $counter;
	}
    }
    return (@ret);
}

sub LoadFFs {
    my ($ffs, $step) = @_;
    my ($i, $PARMS, $counter, $len, $printStr, $ffType);

    $PARMS = ();
    $counter = 0;
    for $i (@{ $ffs }) {
	print "Step $step: " if (defined($step));
	if (ref($i) ne "HASH") {
	    if (! defined($ffType)) {
		$ffType = "CERIUS2";
	    } elsif ($ffType ne "CERIUS2") {
		$ffType = "MULTIPLE";
	    }
	    print "Loading CERIUS2 force field $i...";
            $len = 43 + length($i);
            $PARMS = parseCerius2FF($i, 1, $PARMS);
	} else {
	    print "Loading $i->{FFTYPE} force field $i->{FF}...";
	    $len = 43 + length($i->{FF});
	    if ($i->{FFTYPE} eq "CERIUS2") {
	        $PARMS = &parseCerius2FF($i->{FF});
	    } else {
		$PARMS = &parseMPSimFF($i->{FF});
	    }
            if (! defined($ffType)) {
                $ffType = $i->{FFTYPE};
            } elsif ($ffType ne "CERIUS2") {
                $ffType = "MULTIPLE";
            }

	}
        $counter++;
        print "Done\r";
    }
    die "ERROR: No valid $ffType forcefields found!\n" if (! keys %{ $PARMS });
    $printStr = "Step $step: " if (defined($step));
    $printStr .= "Loading $ffType force field..sucessfully loaded $counter force field";
    $printStr .= "s" if ($counter > 1);
    print "$printStr";
    printf("%" . ($len - length($printStr) + 5) . "s", " ") if ($len > length($printStr));
    print "\n";

    return $PARMS;
}

sub ReadFFs {
    my ($FF) = $_[0];
    my ($i, @FFILES, @tmp, $rec, $isCerius2, $ffType);
    
    if ($FF =~ /\s+/) {
        @tmp = split /\s+/, $FF;
    } else {
        @tmp = ($FF);
    }
                                                                                                                                         
    $ffType = 0;
    for $i (@tmp) {
        $rec = ();
        if (-e $i &&  -r $i && -T $i) {
            $isCerius2 = 0;
            next if (! open(TESTFF, $i));
            while (<TESTFF>) {
                if ($_ =~ /CERIUS2/i) {
                    $isCerius2 = 1;
                    last;
                }
            }
            close TESTFF;
            $rec->{FF} = $i;
            if ($isCerius2) {
                $rec->{FFTYPE} = "CERIUS2";
            } else {
                $rec->{FFTYPE} = "MPSIM";
            }
            if ($#FFILES == 0) {
                $ffType = 1 if ($i =~ /amber/i);
                $ffType = 2 if ($i =~ /charmm/i);
            }
        } elsif (uc($i) eq "AMBER03") {
            $rec = (
                        {
                            "FF"      => "/ul/tpascal/ff/AMBER03.ff",
                            "FFTYPE"  => "CERIUS2",
                        }
                    );
            $ffType = 1;
        } elsif (uc($i) eq "AMBER03_SPC") {
            $rec = (
                        {
                            "FF"      => "/ul/tpascal/ff/AMBER03_SPC.ff",
                            "FFTYPE"  => "CERIUS2",
                        }
                    );
            $ffType = 1;
        } elsif (uc($i) eq "AMBER96") {
            $rec = (
                        {
                            "FF"      => "/ul/tpascal/ff/AMBER96.ff",
                            "FFTYPE"  => "CERIUS2",
                        }
                    );
            $ffType = 1;
        } elsif (uc($i) eq "AMBER99") {
            $rec = (
                        {
                            "FF"      => "/ul/tpascal/ff/AMBER99.ff",
                            "FFTYPE"  => "CERIUS2",
                        }
                    );
            $ffType = 1;
        } elsif (uc($i) eq "AMBER91") {
            $rec = (
                        {
                            "FF"      => "/ul/tpascal/ff/AMBER91.ff",
                            "FFTYPE"  => "CERIUS2",
                        }
                    );
            $ffType = 1;
        } elsif (uc($i) eq "CHARMM_LIPID") {
            $rec = (
                        {
                            "FF"      => "/ul/tpascal/ff/charmm_par_all27_prot_lipid.ff",
                            "FFTYPE"  => "CERIUS2",
                        }
                    );
            $ffType = 2;
        } elsif (uc($i) eq "CHARMM") {
            $rec = (
                        {
                            "FF"      => "/ul/tpascal/ff/charmm_par_all27_prot_na.ff",
                            "FFTYPE"  => "CERIUS2",
                        }
                    );
            $ffType = 2;
        } elsif (uc($i)  eq "MESODNA") {
            #push @FFILES, "/ul/tpascal/MesoDNA/ff/MesoDNA_ver6_mod.par";
            #push @FFILES, "/project/dna_nano/simulations/helices/newFF/MesoDNA_ver7.par";
            $rec = (
                        {
                            "FF"      => "/project/dna_nano/simulations/helices/newFF/MesoDNA_ver7.par",
                            "FFTYPE"  => "CERIUS2",
                        }
                    );
            $ffType = 3;
        } elsif (uc($i) eq "DREIDING") {
            $rec = (
                        {
                            "FF"      => "/ul/tpascal/ff/DREIDING2.21.ff",
                            "FFTYPE"  => "CERIUS2",
                        }
                    );
            $ffType = 4;
        }
        push @FFILES, $rec if ($rec);
    }
    die "ERROR: No valid files found!\n" if (! @FFILES);

    return (\@FFILES, $ffType);
}

sub parseMPSimFF {
    my ($ff_file, $alter, $oldFF) = @_;
    my (%PARMS, $which_var, $in_data, $type_counter, @vdws, $hb_counter, $pi);
    my ($type_id, $type_id2, @bonds, $type_id3, @angles, $tmp1, @inversions, $j);
    my ($type_id4, @torsions, $torsion_type, $inversion_type, $counter, $use_hb);
    my ($bond_counter, $angle_counter, $torsion_counter, @tmp, $inversion_counter, $bool);
    my ($atom1, $atom2, $vdwType, $vdwDat, $i, $dihdr_scale, $CNV, $vdw_counter, $ELEMENTS);

    $CNV = &loadConverter;
    $ELEMENTS = &loadElements;
    $which_var = $bond_counter = $angle_counter = $torsion_counter = $counter = 0;
    $use_hb = $hb_counter = 0;
    if (defined($oldFF)) {
	%PARMS = %{ $oldFF };
    }
    $PARMS{"PARMS"}{"cut_vdw"} = 10.0;
    $PARMS{"PARMS"}{"cut_coul"} = 10.0;
    $PARMS{"PARMS"}{"coul_accuracy"} = 0.0001;
    $PARMS{"PARMS"}{"hbond_distance"} = 2.5;
    $PARMS{"PARMS"}{"hbond_angle"} = 90;
    $PARMS{"PARMS"}{"scale_torsions"} = 0;
    $PARMS{"PARMS"}{"same_scale"} = 1;
    $pi = 3.141592654;

    open FORCEFIELD, $ff_file or die "Cannot open force field file $ff_file: $!\n";
    while (<FORCEFIELD>) {
	chomp;
	$in_data = $_;
        if ($in_data =~ /^(END|AUTOTYPE|GASTEIGER|DEL)/i) {
	    $which_var = 0;
	} elsif ($in_data =~ /^RNB GEOMN\s+(\w)/) {
	    if (lc($1) eq "t") {
		$PARMS{PARMS}{mix_rule} = "geometric";
	    } else {
		$PARMS{PARMS}{mix_rule} = "arithmetic";
	    }
 	} elsif ($in_data =~ /^TORS SCAL+T/) {
	    $PARMS{"PARMS"}{"scale_torsions"} = 1;
	} elsif ($in_data =~ /^BNDBNDTOR\s*(\w)/) {
	    if (lc($1) ne "t") {
		$bool = 0;
	    } else {
		$bool = 1;
	    }
	    $PARMS{PARMS}{scale_vdw_12} = $PARMS{PARMS}{scale_cou_12} = $bool;
        } elsif ($in_data =~ /^ANGANGTOR\s*(\w)/) {
            if (lc($1) ne "t") {
                $bool = 0;
            } else {
                $bool = 1;
            }
            $PARMS{PARMS}{scale_vdw_13} = $PARMS{PARMS}{scale_cou_13} = $bool;
 	} elsif ($in_data =~ /^SCAL NB14\s+(\d+\.\d+)/) {
            $PARMS{"PARMS"}{"scale_coulomb"} = $PARMS{"PARMS"}{"scale_vdw"} = $1;
	    $PARMS{"PARMS"}{"scale_cou_14"} = $PARMS{"PARMS"}{"scale_vdw_14"} = $1;
	} elsif ($in_data =~ /^LHBOND\s+T/) {
	    $use_hb = 1;
	    $PARMS{"PARMS"}{"hbond_distance"} = 4.5;
	    $PARMS{"PARMS"}{"hbond_angle"} = 120;
	} elsif ($in_data =~ /^FFLABEL    ATNO MODIFD/) {
	    $which_var = 1;
	} elsif ($in_data =~ /MPSIM_HB \(A\-H\-D\)  TYPE/) {
		$which_var = 2;	    
	} elsif ($in_data =~ /^VDW AT ITY       RNB      DENB     SCALE/) {
		$which_var = 3;
	} elsif ($in_data =~ /^BONDSTRTCH  TYPE/) {
	    $which_var = 4;
	} elsif ($in_data =~ /^ANGLE\-\(L\-C\-R\)     TYPE/) {
	    $which_var = 5;
	} elsif ($in_data =~ /^TORSION                 CASE   BARRIER/) {
	    $which_var = 6;
	} elsif ($in_data =~ /^INVERSION \(CENT AT 1ST\)/) {
	    $which_var = 7;
	} elsif ($in_data =~ /^NONBOND-OFF TYPE      RVDW/) {
	    $which_var = 8;
	} elsif ($in_data =~ /^UREY_BRADLEY/) {
	    $which_var = 9;
	} elsif ($in_data =~ /^\s*(\S+)\s+(\d+)\s+(\d*\.?\d*)\s+(\-?\d+\.\d+)\s+(\-?\d+)\s+(\-?\d+)\s+(\-?\d+)/ and ($which_var == 1)) {
	    $type_counter += 1;
	    $PARMS{"ATOMTYPES"}{$1} = (
				       {
					   "TYPEID"    => $type_counter,
					   "ATOM"      => $ELEMENTS->{$2}{SYMBOL},
					   "MASS"      => $ELEMENTS->{$2}{MASS},
					   "CHARGE"    => $4,
					   "NUMBONDS"  => $5,
					   "LONEPAIRS" => $6,
					   "OTHER"     => $7,
					   "LABEL"     => $1,
					   "USED"      => 0,
				       }
				       );
#	    print "ATOMTYPES: $1: $type_counter\n";
        } elsif ($in_data =~ /^\s*(\S+)\s*\-\s*(\S+)\s*\-\s*(\S+)\s+(\d+)\s+(\-?\d+\.\d+e?\-?\d*)\s+(\d+\.\d+e?\-?\d*)/ and $which_var == 2) {
	    $atom1 = $2;
	    $atom2 = $3;
	    $vdwType = "LJ_12_10";
	    @vdws = ($5, $6);
	    @vdws = GetAbs(\@vdws);
	    if (exists($PARMS{ATOMTYPES}{$atom1}) and exists($PARMS{ATOMTYPES}{$atom2})) {
                if (! exists($PARMS{VDW}{$atom1}{$atom2})) {
                    $vdw_counter = 1;
                } else {
                    $vdw_counter = scalar(keys %{ $PARMS{VDW}{$atom1}{$atom2} }) + 1;
                }
                $tmp1 = "hbond";
                $PARMS{"VDW"}{$atom1}{$atom2}{$vdw_counter} = (
                                                 {
                                                     "TYPE"   => $vdwType,
                                                     "Lammps" => getLammpsOpts($vdwType,$tmp1, $CNV),
                                                     "KEY"    => "$atom1 $atom2 ",
                                                     "VALS"   => [@vdws],
                                                     "ATOM"   => "$atom1 $atom2 ",
                                                     "USED"   => 0,
                                                     "IT"     => $tmp1,
                                                 }
                                                 );
	    }
	} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*(.*)/i and ($which_var == 3 || $which_var == 8)) {
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
	    if ($vdwType == 1 and $which_var == 3) {
		$vdwType = "EXPO_6";
	    } else {
		$vdwType = "LJ_6_12";
	    }
	    if (uc($vdwType) ne "IGNORE") {
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$1}));

		$type_id = $PARMS{"ATOMTYPES"}{$atom1}{"TYPEID"};
		if ($vdwDat =~ /(.*)\# 1\-4 scaling: (.*)/) {
		    @vdws = split /\s+/, "$1 $2";
		} else {
		    @vdws = split /\s+/, $vdwDat;
		}
		@vdws = GetAbs(\@vdws);

		if (! defined($alter) || $alter != 0) {
		    if (uc($vdwType) eq "VDW_MORSE") {
			($vdws[0], $vdws[1]) = ($vdws[1], $vdws[0]);
			($vdws[1], $vdws[2]) = ($vdws[2], $vdws[1]);
			$vdws[1] /= ($vdws[2] * 2); # changed from div to multi 07/28/2007
			#$vdws[0] *= 2;
			#$vdws[1] /= ($vdws[2]/2);
		    } elsif (uc($vdwType) eq "LJ_6_12") {
			@vdws = splice(@vdws, 0, 2);
			($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
			$vdws[1] = $vdws[1] / (2**(1/6));
#			$vdws[1] = $vdws[1]/2;
		    } elsif (uc($vdwType) eq "EXPO_6") {
			($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
			@tmp = @vdws;
			@vdws = ();
			#$vdws[0] = 6 * $tmp[1] * exp($tmp[2]) / ($tmp[2] - 6);
			#$vdws[1] = $tmp[0]/$tmp[2];
			#$vdws[2] = $tmp[1] * $tmp[2] * $tmp[0]**6/($tmp[2] - 6);
			$vdws[0] = $tmp[0] * (6/($tmp[2]-6)) * exp($tmp[2]);
			$vdws[1] = $tmp[1]/$tmp[2];
			$vdws[2] = $tmp[1]**6 * $tmp[0] * ($tmp[2]/($tmp[2] - 6));
		    } 
		}
		if (! exists($PARMS{VDW}{$atom1}{$atom2})) {
		    $vdw_counter = 1;
		} else {
		    $vdw_counter = scalar(keys %{ $PARMS{VDW}{$atom1}{$atom2} }) + 1;
		}
		$tmp1 = "vdw";
		$PARMS{"VDW"}{$atom1}{$atom2}{$vdw_counter} = (
						 {
						     "TYPE"   => $vdwType,
						     "Lammps" => getLammpsOpts($vdwType,$tmp1, $CNV),
						     "KEY"    => "$atom1 $atom2 ",
						     "VALS"   => [@vdws],
						     "ATOM"   => "$atom1 $atom2 ",
						     "USED"   => 0,
						     "IT"     => $tmp1,
						 }
						 );
#		print "VDW: $1: $type_id\n"
	    }
	} elsif ($in_data =~ /^\s*(\S+)\s+\-(\S+)\s+(\d+)\s+(.+)/ and ($which_var == 4)) {
	    if ($3 and $3 < 3) { # can only handle morse and harmonic
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
		
		$bond_counter++;
		@bonds = split /\s+/, $4;
		next if (! @bonds);
#		GetAbs(\@bonds);
		if (! defined($alter) || $alter != 0) {
                    $bonds[0] = $bonds[0] / 2; #fix for the 1/2 factor in harmonic eqns between cerius and lammps
		    if ($3 == 2) {
			$bonds[2] = sqrt($bonds[0]/(2 * $bonds[2]));
			($bonds[1], $bonds[2]) = ($bonds[2], $bonds[1]);
			pop @bonds;
		    }
		}
		$i = "HARMONIC";
		$i = "MORSE" if ($3 == 2);
		if ($i eq "HARMONIC" and $#bonds > 1) {
		    @bonds = splice(@bonds, 0, 2);
		} elsif ($i eq "MORSE" and $#bonds > 2) {
		    @bonds = splice(@bonds, 0, 3);
		}
		$PARMS{"BONDS"}{$1}{$2} = (
					   {
					       "INDEX"    => $bond_counter,
					       "TYPE"     => $i,
					       "Lammps"   => getLammpsOpts($i,"bond"),
					       "VALS"     => [@bonds],
					       "USED"     => 0,
					       "KEY"      => "$1 $2 ",
					   }
					   );

#		print "BOND $bond_counter: $key_code\n";
	    }
	} elsif ($in_data =~ /^\s*(\S+)\s+\-(\S+)\s+\-(\S+)\s+(\d+)\s+(.+)/ and ($which_var == 5)) {
	    if ($4) {
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$3}) and lc($3) ne "x");
		
		$angle_counter++;
		@angles = split /\s+/, $5;
		next if (! @angles);
#		GetAbs(\@angles);
		if (! defined($alter) || $alter != 0) {
	            $angles[0] = $angles[0]/2; #same fix as bonds
		}
		if ($4 == 11) {
		    if (sin($angles[1] * $pi/180) > 0.001) {
			$angles[0] /= sin($angles[1] * $pi/180)**2;
		    }
		}
		$i = "HARMONIC";
		$i = "COS_HARMON" if ($4 == 11);
                if ($#angles > 1) {
                    @angles = splice(@angles, 0, 2);
                }

		$PARMS{"ANGLES"}{$1}{$2}{$3} = (
						{
						    "INDEX"    => $angle_counter,
						    "TYPE"     => $i,
						    "Lammps"   => getLammpsOpts($i,"angle", $CNV),
						    "VALS"     => [@angles],
						    "USED"     => 0,
						    "KEY"      => "$1 $2 $3 ",
						}
						);
#		print "ANGLE $angle_counter: $key_code\n";
	    }
	} elsif ($in_data =~ /^\s*(\S+)\s+\-(\S+)\s+\-(\S+)\s+\-(\S+)\s+(\d+)\s+(.+)/ and ($which_var == 6)) {
	    $torsion_type = "DIHEDRAL";
	    if (uc($torsion_type) ne "IGNORE") {
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$3}) and lc($3) ne "x");
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$4}) and lc($4) ne "x");

		$torsion_counter++;
                @torsions = split /\s+/, $6;
#		GetAbs(\@torsions);
 
                if (! defined($alter) || $alter != 0) {
		    if ($torsions[2] == 1) {
		        $torsions[2] = 180;
		    } else {
		        $torsions[2] = 0;
                    }
		}
		#$torsions[2] = 360 - $torsions[2];
		if ($#torsions > 2) {
		    @torsions = splice(@torsions, 0, 3);
		}
		if (! exists($PARMS{"TORSIONS"}{$1}{$2}{$3}{$4})) {
		    $PARMS{"TORSIONS"}{$1}{$2}{$3}{$4} = (
						      {
							  "TYPE"     => $torsion_type,
							  "Lammps"   => getLammpsOpts($torsion_type,"dihedral", $CNV),
							  "INDEX"    => $torsion_counter,
							  "VALS"     => [@torsions],
							  "USED"     => 0,
							  "KEY"      => "$1 $2 $3 $4 ",
							  "NUM"      => 1,
							  "PER"      => ($#torsions + 1),
							  "1_4scale" => $PARMS{"PARMS"}{"scale_vdw_14"},
						      }
						      );
		} else {
	            push @{ $PARMS{"TORSIONS"}{$1}{$2}{$3}{$4}{"VALS"} }, @torsions;
        	    $PARMS{"TORSIONS"}{$1}{$2}{$3}{$4}{NUM}++;
		}
#		print "TORSION $torsion_counter: $key_code\n";
	    }
	} elsif ($in_data =~ /^\s*(\S+)\s+\-(\S+)\s+\-(\S+)\s+\-(\S+)\s+(\d+)\s+(.+)/ and ($which_var == 7)) {
	    if ($5 == 1) {
		$inversion_type = "IT_IJKL";
	    }elsif ($5 == 3) {
		$inversion_type = "IT_JIKL";
	    } else {
		$inversion_type = "UMBRELLA";
	    }
	    if (uc($inversion_type) ne "IGNORE") {
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$3}) and lc($3) ne "x");
		next
		    if (! exists($PARMS{"ATOMTYPES"}{$4}) and lc($4) ne "x");

		$inversion_counter++;
                @inversions = split /\s+/, $6;
		next if (! @inversions);
#		GetAbs(\@torsions);
 
                if (! defined($alter) || $alter != 0) {
                    if (lc($inversion_type) eq "it_jikl") {
			$type_id = $2;
			$type_id2 = $3;
			$type_id3 = $1;
			$type_id4 = $4;
			if ($inversions[1] == 180) {
			    $inversions[1] = -1;
			} else {
			    $inversions[1] = 1;
			}
                    } else {
			($type_id, $type_id2, $type_id3, $type_id4)  = ($1, $2, $3, $4);
	   	    }
		    $inversions[0] = $inversions[0] / 2; #1/2 fix again
		}
		$PARMS{"INVERSIONS"}{$1}{$2}{$3}{$4} = (
							{
							    "TYPE"     => $inversion_type,
							    "Lammps"   => getLammpsOpts($inversion_type,"inversion", $CNV),
							    "INDEX"    => $inversion_counter,
							    "VALS"     => [@inversions],
							    "USED"     => 0,
							    "KEY"      => "$1 $2 $3 $4 ",
							}
							);
#		print "TORSION $torsion_counter: $key_code\n";
	    }
	}	    
	
    }
    close FORCEFIELD;

    die "ERROR: Force Field file $ff_file is invalid!\n"
	if (! %PARMS);
    
    if ($PARMS{"PARMS"}{"cut_vdw"} > $PARMS{"PARMS"}{"cut_coul"}) {
	$PARMS{"PARMS"}{"cut_vdw"} = $PARMS{"PARMS"}{"cut_coul"};
    }

    return (\%PARMS);
}

sub loadElements {
    my (%ELEMENTS, $indata, $eleNum, $myDir);
    $myDir = "/ul/tpascal/scripts/Packages";
    my ($datFile) = "$myDir/elementList.txt";
    die "ERROR: $!\n" if (! -e $datFile or ! -r $datFile or ! -T $datFile);
                                                                                                                                         
    open INDATA, $datFile or die "Cannot open file $datFile: $!\n";
    while (<INDATA>) {
        chomp;
        $indata = $_;
        if ($indata =~ /^(\d+)\s+(\*?)\s*\d+/) {
            $eleNum = $1;
            if ($2) {
                $ELEMENTS{$1}{"NATURAL"} = 0;
            } else {
                $ELEMENTS{$1}{"NATURAL"} = 1;
            }
            #$indata = $';
            if ($indata =~ /^\d+\s+\*?\s*(\d+\.?\d*)\s+(\w+)\s+(\w+)/) {
                $ELEMENTS{$eleNum}{"MASS"} = $1;
                $ELEMENTS{$eleNum}{"NAME"} = $2;
                $ELEMENTS{$eleNum}{"SYMBOL"} = $3;
            } else {
                delete $ELEMENTS{$eleNum};
            }
        }
    }
                                                                                                                                         
    close INDATA;
                                                                                                                                         
    die "ERROR: No valid data found in file $datFile\n"
        if (! %ELEMENTS);

    return \%ELEMENTS;

}
1;

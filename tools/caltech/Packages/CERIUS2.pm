package Packages::CERIUS2;

require Exporter;
use strict;
use Cwd 'abs_path';
use Math::Trig qw(pi);

our(@EXPORT_OK, @ISA, @EXPORT, $VERSION);

my $scripts_dir = "/home/tao/Nutstore/code/simupy/tools/caltech"; #change here as necessary
my $ff_dir = "/ul/tpascal/ff"; #change here as necessary
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(parseCerius2FF saveCeriusFF ReadFFs parseMPSimFF parseReaxFF);
$VERSION = "1.00";

sub loadConverter {
	my ($inFile) = "${scripts_dir}/dat/ceriusParms2Lammps.perldata";
	my (%SHASH);
	scalar eval `cat $inFile` or die "Cannot recreate data in file $inFile: $! $@\n";
	return \%SHASH;
}

sub getLammpsOpts {
	my ($parmName, $parmType, $SHASH, $istip4p, $parms) = @_;
	my (%RET, $same_cut);
	my ($searchStr) = lc($parmType) . "_" . lc($parmName);

	$istip4p = 0 if (! defined($istip4p));
	$searchStr .= "_tip4p" if ($parmType eq "vdw" and $istip4p and $parmName eq "LJ_6_12");
	
	$same_cut = 1;
	#$same_cut = 1 if ($parms->{cut_coul} == $parms->{cut_vdw});
	if (exists($SHASH->{$searchStr})) {
		%RET = %{ $SHASH->{$searchStr} };
	} else {
		$RET{name} = lc($parmName);
		$RET{opts} = "";
		$RET{missing} = 1;
	}

	if($parmType eq "vdw") {
		if($RET{name} =~ /charmm/) {
			if (! $same_cut) {
				$RET{opts} = ($parms->{cut_vdw}-1) . " $parms->{cut_vdw} " . ($parms->{cut_coul}-1) . " $parms->{cut_coul}";
			} else {
				$RET{opts} = ($parms->{cut_coul}-1) . " $parms->{cut_coul}";
			}
		} elsif ($RET{name} =~ /lj|gromacs/ && $RET{name} !~ /dreiding/) {
			$RET{opts} = " $parms->{cut_vdw}  $parms->{cut_coul}";
			$RET{opts} = $parms->{cut_coul} if ($same_cut);
		}
	}
	return \%RET;
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

sub parseReaxFF {
	my ($ff_file, $elements, $alter, $oldFF) = @_;
	my (%PARMS, $type_counter, $currParm);

	$type_counter = 0;
	if (defined($oldFF)) {
		%PARMS = %{ $oldFF };
	}

	open FORCEFIELD, $ff_file or die "Cannot open $ff_file: $!\n";
	while (<FORCEFIELD>) {
		chomp;
		if ($_ =~ /^\s+([A-Za-z][A-Za-z]?)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)/) { # atom line
			$type_counter++;
			$PARMS{"ATOMTYPES"}{$1} = (
									   {
										"TYPEID"		=>	$type_counter,
										"ATOM"			=>	$1,
										"MASS"			=>	$4,
										"CHARGE"		=>	undef,
										"NUMBONDS"		=>	$3,
										"LONEPAIRS"		=>	undef,
										"OTHER"			=>	undef,
										"LABEL"			=>	$1,
										"USED"			=>	0,
										"USE_CHARGE"	=>	0,
									   }
									);
#		   print "ATOMTYPES: $1: $type_counter\n";
				$currParm = \%{ $PARMS{"VDW"}{$1}{$1} };
				$currParm->{1} =(
									{
										"TYPE"		=>	"REAX",
										"Lammps"	=>	{
														"name"	=>	"reax",
														"opts"	=>	"",
														},
										"KEY"		=>	"$1 $1 ",
										"VALS"		=>	[0],
										"ATOM"		=>	"$1 $1 ",
										"USED"		=>	0,
										"IT"		=>	"vdw",
									}
								);
#				print "VDW: $1: $type_id\n"

		}
	}
	close FORCEFIELD;
	die "ERROR: Force Field file $ff_file is invalid!\n"
		if (! %PARMS);
	return (\%PARMS);
}

sub findElementNum {
	my ($symbol, $elements) = @_;
	my ($i, $elementNum);

	$symbol = lc $symbol;
	$elementNum = 0;
	for $i (keys %{ $elements }) {
		if(lc($elements->{$i}{SYMBOL}) eq $symbol) {
			$elementNum = $i;
			last;
		}
	}

	return $elementNum;
}

sub getHbondList {
	my ($parms, $donor, $acceptor, $hydrogen, $elements) = @_;
	my ($hbondList, $i, $j, $k, $iNum, $jNum, $dList); 
	my ($aList, $hList, $rec, $stored);

	$hbondList = ();

	if($donor ne "X" && $acceptor ne "X" && $hydrogen ne "X") {
		$hbondList = [{ "donor" => $donor, "acceptor" => $acceptor, "hydrogen" => $hydrogen, }];
		return $hbondList;
	}

	push @{ $dList }, $donor if($donor ne "X");
	push @{ $hList }, $hydrogen if($hydrogen ne "X");
	push @{ $aList }, $acceptor if($acceptor ne "X");

	#first search for donors and hydrogens
	for $i (keys %{ $parms->{BONDS} }) {
		next if($i =~ /^\s*H_\s*$/);
		$iNum = findElementNum($parms->{ATOMTYPES}{$i}{ATOM}, $elements);
		for $j (keys %{ $parms->{BONDS}{$i} }) {
			next if($j =~ /^\s*H_\s*$/);
			$jNum = findElementNum($parms->{ATOMTYPES}{$j}{ATOM}, $elements);
			next if((! $iNum && ! $jNum) || ($iNum != 1 && $jNum != 1));
			next if($donor ne "X" && ($i ne $donor && $j ne $donor));
			if($iNum == 1) {
				push @{ $hList }, $i if($hydrogen eq "X");
				push @{ $dList }, $j if($donor eq "X" && $jNum =~ /^(7|8|9)$/);
			} else {
				push @{ $hList }, $j if($hydrogen eq "X");
				push @{ $dList }, $i if($donor eq "X" && $iNum =~ /^(7|8|9)$/);
			}
		}
	}

	#now search for acceptors
	if($acceptor eq "X") {
		for $i (keys %{ $parms->{ATOMTYPES} }) {
			$iNum = findElementNum($parms->{ATOMTYPES}{$i}{ATOM}, $elements);
			next if ($iNum !~ /^(6|7|8|9|16|17)$/); # not C,N,O,F,S,Cl
			next if ($parms->{ATOMTYPES}{$i}{LONEPAIRS} == 0); # must have lonepairs
			push @{ $aList }, $i;
		}
	}

	for $i (@{ $dList }) {
		for $j (@{ $aList }) {
		   #next if ($i eq $j);
			for $k (@{ $hList }) {
				next if ($k !~ /H___|H_F/);
				next if (exists($stored->{$i}{$j}{$k}));
				$stored->{$i}{$j}{$k} = 1;
				$rec = ({ "donor" => $i, "acceptor" => $j, "hydrogen" => $k, });
				push @{ $hbondList }, $rec;
		   }
		}
	}

	return $hbondList;
}

sub searchForHBparm {
	my ($currParm, $hatom) = @_;
	my ($index, $i, $count);

	$index = $count = 0;
	for $i (keys %{ $currParm }) {
		$count = $i;
		next if (! exists($currParm->{$i}{HATOM}));
		if($currParm->{$i}{HATOM} eq $hatom) {
			$index = $i;
			last;
		}
	}

	$index = $count + 1 if(! $index );

	return $index;
}

sub getVals {
	my ($parmStr) = $_[0];
	my (@ret);

	$parmStr =~ s/#.*$//;
	@ret = split /\s+/,$parmStr;

	return @ret;
}
	
sub parseCerius2FF {
	my ($ff_file, $alter, $oldFF) = @_;
	my (%PARMS, $which_var, $in_data, $type_counter, @vdws, $hb_counter);
	my (@bonds, @angles, $tmp1, @inversions, $atom4, $currParm, $calcExp6, $istip4p);
	my (@torsions, $torsion_type, $inversion_type, $counter, $use_hb, $use_charge);
	my ($atom3, $bond_counter, $angle_counter, $torsion_counter, @tmp, $inversion_counter, $bool);
	my ($atom1, $atom2, $vdwType, $vdwDat, $i, $dihdr_scale, $CNV, $vdw_counter, $crossType);
	my ($donors, $acceptors, $hydrogens, $elements, $tmp, $scale_flag);

	$CNV = &loadConverter;
	$elements = loadElements();
	$which_var = $bond_counter = $angle_counter = $torsion_counter = $counter = 0;
	$use_hb = $hb_counter = $crossType = 0;
	$calcExp6 = 1;
	$istip4p = 0;
	$istip4p = 1 if ($ff_file =~ /tip4/);
	if (defined($oldFF)) {
		%PARMS = %{ $oldFF };
	} else {
		$PARMS{PARMS}{SW} = 0;
		$PARMS{"PARMS"}{"cut_vdw"} = 14.0;
		$PARMS{"PARMS"}{"cut_coul"} = 15.0;
		$PARMS{"PARMS"}{"coul_accuracy"} = 0.001;
		$PARMS{"PARMS"}{"hbond_distance"} = 2.5;
		$PARMS{"PARMS"}{"hbond_angle"} = 90;
		$PARMS{"PARMS"}{"scale_torsions"} = 0;
		$PARMS{PARMS}{scale_torsions_by_n_defined_torsions} = 0;
		$PARMS{"PARMS"}{"single_inversion"} = 0;
		$PARMS{"PARMS"}{"dielectric"} = 1;
		$PARMS{"PARMS"}{"EXO_CYCLIC"} = 1.0;
		$PARMS{"PARMS"}{QEq} = 0;
		$PARMS{PARMS}{USE_HBOND} = 0;
		$PARMS{PARMS}{HAS_HBONDS} = 0;
	}
	$PARMS{"PARMS"}{"tip4_om_dist"} = 0.125 if ($istip4p);
	$PARMS{"PARMS"}{"is_tip4p"} = 1 if ($istip4p);	
	$type_counter = 0;
	$type_counter = ($PARMS{"PARMS"}{"type_counter"} + 1) if (exists($PARMS{"PARMS"}{"type_counter"}));
	$use_charge = 0;

	open FORCEFIELD, $ff_file or die "Cannot open force field file $ff_file: $!\n";
	while (<FORCEFIELD>) {
		chomp;
		$in_data = $_;
		$in_data =~ s/#.*$// if ($which_var != 3 and $which_var != 8);
		$PARMS{PARMS}{SW} = 1 if ($in_data =~ / SW /);
		if ($in_data =~ /^END/) {
			$which_var = 0;
		} elsif ($in_data =~ /^\s+EXOCYCLIC_TORSIONS_SCALE_FACTOR\s+(\d+\.?\d*)/) {
			$PARMS{"PARMS"}{EXO_CYCLIC} = $1;
		} elsif ($in_data =~ /^\s+USE_CURRENT_EXP6_VALS\s+T/) {
			$calcExp6 = 0;
		} elsif ($in_data =~ /^\s*COU_DIELETRIC_CONSTANT\s+(\d+\.\d+)/) {
			$PARMS{"PARMS"}{"dielectric"} = $1;
		} elsif ($in_data =~ /^\s*SINGLE_INVERSION\s+T/) {
			$PARMS{"PARMS"}{"single_inversion"} = 1;
		 } elsif ($in_data =~ /^\s*SCALE_TORSIONS_ABOUT_COMMON_BOND\s+T/) {
			$PARMS{"PARMS"}{"scale_torsions"} = 1;
		 } elsif ($in_data =~ /^\s*SCALE_BY_N_DEFINED_TORSIONS\s+T/) {
			$PARMS{"PARMS"}{"scale_torsions_by_n_defined_torsions"} = 1;
		} elsif ($in_data =~ /^\s*COU_SPLINE_OFF\s+(\d+\.\d+)/) {
			$PARMS{"PARMS"}{"cut_coul"} = $1;
		} elsif ($in_data =~ /^\s*VDW_SPLINE_OFF\s+(\d+\.\d+)/) {
			$PARMS{"PARMS"}{"cut_vdw"} = $1;
		} elsif ($in_data =~ /^\s*VDW_COMBINATION_RULE\s+(\w+)/) {
			$PARMS{"PARMS"}{"mix_rule"} = lc($1);
		} elsif ($in_data =~ /^\s*EWALD_SUM_COU_ACCURACY\s+(\d+\.\d+)/) {
			$PARMS{"PARMS"}{"coul_accuracy"} = $1;
		} elsif ($in_data =~ /^\s*(COU|VDW)_EXCLUDE_(\d)\-(\d)\s+(\w)/) {
			if (lc($4) eq "t") {
				$bool = 0;
			} else {
				$bool = 1;
			}
			$PARMS{"PARMS"}{"scale_" . lc($1) . "_" . $2 . $3} = $bool;
		 } elsif ($in_data =~ /^\s*COU_1-4_SCALE_FACTOR\s+(\d+\.\d+)/) {
			$PARMS{"PARMS"}{"scale_cou"} = $1;
			$PARMS{"PARMS"}{"scale_cou_14"} = $1 if ($PARMS{"PARMS"}{"scale_cou_14"})
		 } elsif ($in_data =~ /^\s*VDW_1-4_SCALE_FACTOR\s+(\d+\.\d+)/) {
			$PARMS{"PARMS"}{"scale_vdw"} = $1;
			$PARMS{"PARMS"}{"scale_vdw_14"} = $1 if ($PARMS{"PARMS"}{"scale_vdw_14"});
			if($PARMS{"PARMS"}{"scale_cou"} ne $PARMS{"PARMS"}{"scale_vdw"}) { 
				$PARMS{"PARMS"}{"same_scale"} = 0;
			} else {
				$PARMS{"PARMS"}{"same_scale"} = 1;
			}
		} elsif ($in_data =~ /^ HYDROGEN_BONDS\s+T/) {
			$use_hb = 1;
			$PARMS{"PARMS"}{"hbond_distance"} = 4.5;
			$PARMS{"PARMS"}{"hbond_angle"} = 60;
			$PARMS{PARMS}{USE_HBOND} = 1;
		} elsif ($use_hb and $in_data =~ /^ H-BOND_LIST_DISTANCE_OFF\s+(\d+\.\d+)/) {
			$PARMS{"PARMS"}{"hbond_distance"} = $1;
		} elsif ($use_hb and $in_data =~ /^ H-BOND_LIST_ANGLE_OFF\s+(\d+\.\d+)/) {
			$PARMS{"PARMS"}{"hbond_angle"} = $1;
		} elsif ($in_data =~ /^ ASSIGN_CHARGE\s+T/) {
			$use_charge  = 1;
		} elsif ($in_data =~ /^ TIP4_OM_DIST\s+(\d+\.\d+)/) {
			$PARMS{"PARMS"}{"tip4_om_dist"} = $1;
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
		} elsif ($in_data =~ /^SEPARATED_STRETCH_STRETCH/) {
			$which_var = 6;
			$crossType = "13BondBond";
		} elsif ($in_data =~ /^TORSION_BEND_BEND/) {
			$which_var = 6;
			$crossType = "AngleAngle";
		} elsif ($in_data =~ /^BEND_BEND/) {
			$which_var = 7;
			$crossType = "AngleAngle";
		} elsif ($in_data =~ /^GENERATOR/) {
			$which_var = 11;
		} elsif ($in_data =~ /^QEq/) {
			$which_var = 10;
		} elsif ($in_data =~ /^\s*(\S+)\s+(.+)$/ and ($which_var == 11)) {
			#UFF Generator
			@tmp = split /\s+/, $2;
			next if ($#tmp < 9);
			$PARMS{PARMS}{UFF}{$1} = (
										{
										"radius"		=> $tmp[0],
										"angle"			=> $tmp[1],
										"zstar"			=> $tmp[2],
										"zeta"			=> $tmp[3],
										"uenerg"		=> $tmp[4],
										"uang"			=> $tmp[5],
										"prd"			=> $tmp[6],
										"cis"			=> $tmp[7],
										"torbar"		=> $tmp[8],
										"elecneg"		=> $tmp[9],
										}
									);
		} elsif ($in_data =~ /^\s*(\S+)\s+(\d+\.\d+\e?\-?\d*)\s+(\d+\.\d+\e?\-?\d*)\s+(\d+\.\d+\e?\-?\d*)/ and ($which_var == 10)) {
			$PARMS{PARMS}{QEq} = 1;
			$PARMS{QEq}{$1} = (
								{
									"CHI"   => $2,
									"ETA"   => $3,
									"GAMMA" => $4,
									"USED"  => 0,
								}
							);
		} elsif ($in_data =~ /^\s*(\S+)\s+(\w+)\s+(\d+\.\d+)\s+(\-?\d+\.\d+)\s+(\d+)\s+(\d+)\s+(\d+)/ and ($which_var == 1)) {
			$type_counter += 1;
			$PARMS{"ATOMTYPES"}{$1} = (
										{
											"TYPEID"		=> $type_counter,
											"ATOM"			=> $2,
											"MASS"			=> $3,
											"CHARGE"		=> $4,
											"NUMBONDS"		=> $5,
											"LONEPAIRS"		=> $7,
											"OTHER"			=> $6,
											"LABEL"			=> $1,
											"USED"			=> 0,
											"USE_CHARGE"	=> $use_charge,
										}
									);
#			print "ATOMTYPES: $1: $type_counter\n";
		} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+.*)/i and $which_var == 2) { #hbonds
			$donors = $hydrogens = $acceptors = "";
			$donors = $1 if(exists($PARMS{"ATOMTYPES"}{$1}) || $1 eq "X");
			$tmp = getLammpsOpts($3,"hbond", $CNV);
			if(! exists($PARMS{ATOMTYPES}{$3}) && $3 ne "X" && ! exists($tmp->{missing})) { #old style hbonds, only donor and acceptor
				$acceptors = $2 if(exists($PARMS{ATOMTYPES}{$2}) || $2 eq "X");
				$hydrogens = "X";
				$vdwType = $3;
				@vdws = getVals("$4 $5");
			} else { #new hbond: donor-hydrogen ... aceptor
				$acceptors = $3 if(exists($PARMS{ATOMTYPES}{$3}) || $3 eq "X");
				$hydrogens = $2 if(exists($PARMS{ATOMTYPES}{$2}) || $2 eq "X");
				$vdwType = $4;
				@vdws = getVals($5);
			}
			next if (! $donors || ! $acceptors || ! $hydrogens);
			$tmp = getHbondList(\%PARMS, $donors, $hydrogens, $acceptors, $elements);
			next if (! $tmp);
			push(@vdws, "4") if($vdwType eq "LJ_12_10"); # append the power for the cosine - 4 - to the end of the data
			for $i (@{ $tmp }) {
				$atom1 = $i->{donor};
				$atom2 = $i->{acceptor};
				if (! exists($PARMS{VDW}{$atom1}{$atom2})) {
					$vdw_counter = 1;
				} else {
					$vdw_counter = searchForHBparm($PARMS{VDW}{$atom1}{$atom2}, $i->{hydrogen});
				}
				$currParm = \%{ $PARMS{"VDW"}{$atom1}{$atom2} };
				$currParm->{$vdw_counter} = (
												{
													"TYPE"		=> $vdwType,
													"Lammps"	=> getLammpsOpts($vdwType,"hbond", $CNV),
													"KEY"		=> "$atom1 $atom2 ",
													"VALS"		=> [$i->{hydrogen},@vdws],
													"ATOM"		=> "$atom1 $atom2 ",
													"USED"		=> 0,
													"IT"		=> "hbond",
													"HATOM"		=> $i->{hydrogen},
												}
											);
				print "";
			}

		} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s*(\S*)\s*(.*)/i and ($which_var == 3 or $which_var == 8)) { #vdw
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
			if (uc($vdwType) =~ /IGNORE/) {
				if (! exists($PARMS{VDW}{$atom1}{$atom2})) {
					$vdw_counter = 1;
				} else {
					$vdw_counter = findDuplicate($PARMS{VDW}{$atom1}{$atom2}, $vdwType);
					$vdw_counter = scalar(keys %{ $PARMS{VDW}{$atom1}{$atom2} }) + 1 if (! $vdw_counter);
				}
				$tmp1 = "vdw";
				$currParm = \%{ $PARMS{"VDW"}{$atom1}{$atom2} };
				$currParm->{$vdw_counter}{USED} = 1;
				$currParm->{$vdw_counter}{IGNORE} = 1;
				$currParm->{$vdw_counter}{VALS} = ();
			} elsif (uc($vdwType) !~ /UFF_GEN/ and defined($4)) {
				next if (! exists($PARMS{"ATOMTYPES"}{$atom1}));
				next if (! exists($PARMS{"ATOMTYPES"}{$atom2}));

				if ($vdwDat =~ /(.*)\# 1\-4 scaling: (.*)/) {
					@vdws = getVals("$1 $2");
				} else {
					@vdws = getVals($vdwDat);
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
#						$vdws[1] = $vdws[1]/2;
					} elsif (uc($vdwType) eq "EXPO_6" and $calcExp6) {
						($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
						@tmp = @vdws;
						@vdws = ();
						#$vdws[0] = 6 * $tmp[1] * exp($tmp[2]) / ($tmp[2] - 6);
						#$vdws[1] = $tmp[0]/$tmp[2];
						#$vdws[2] = $tmp[1] * $tmp[2] * $tmp[0]**6/($tmp[2] - 6);
						$vdws[0] = $tmp[0] * (6/($tmp[2]-6)) * exp($tmp[2]);
						$vdws[1] = $tmp[1]/$tmp[2];
						$vdws[2] = $tmp[1]**6 * $tmp[0] * ($tmp[2]/($tmp[2] - 6));
					} elsif (uc($vdwType) eq "BUCKINGHAM") {
						#($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
						$vdws[1] = 1/$vdws[1];
					} elsif (uc($vdwType) eq "LJ_6_9") {
						($vdws[0],$vdws[1]) = ($vdws[1],$vdws[0]);
						$vdws[1] = $vdws[1] / 1.14471438354;
					} elsif (uc($vdwType) eq "BORN") {
						splice(@vdws,2,0,0);
						$vdws[$#vdws] *= -1;
					}
				}
				if (! exists($PARMS{VDW}{$atom1}{$atom2})) {
					$vdw_counter = 1;
				} else {
					$vdw_counter = findDuplicate($PARMS{VDW}{$atom1}{$atom2}, $vdwType);
					$vdw_counter = scalar(keys %{ $PARMS{VDW}{$atom1}{$atom2} }) + 1 if (! $vdw_counter);
				}
				$tmp1 = "vdw";
				$currParm = \%{ $PARMS{"VDW"}{$atom1}{$atom2} };
				$currParm->{$vdw_counter} = (
												{
													"TYPE"		=> $vdwType,
													"Lammps"	=> getLammpsOpts($vdwType,"vdw", $CNV, $istip4p, $PARMS{PARMS}),
													"KEY"		=> "$atom1 $atom2 ",
													"VALS"		=> [@vdws],
													"ATOM"		=> "$atom1 $atom2 ",
													"USED"		=> 0,
													"IGNORE"	=> 0,
													"IT"		=> "vdw",
												 }
											);
				if ($vdwType eq "YUKAWA") {
					$currParm->{$vdw_counter}{Lammps}{opts} = "$vdws[1] $vdws[2]";
					pop @{ $currParm->{$vdw_counter}{VALS} };
					pop @{ $currParm->{$vdw_counter}{VALS} };
				}
#				print "VDW: $1: $type_id\n"
			}
		} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\w+)\s+(.+)/ and ($which_var == 4)) {
			if (uc($3) !~ /IGNORE|UFF_GEN/) {
				next if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
				next if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
				
				@bonds = getVals($4);
				next if (! @bonds);
#				GetAbs(\@bonds);
				if (! defined($alter) || $alter != 0) {
					$bonds[0] = $bonds[0] / 2; #fix for the 1/2 factor in harmonic eqns between cerius and lammps
					if (uc($3) eq "MORSE") {
						$bonds[0] *= 2;
						($bonds[0], $bonds[2]) = ($bonds[2], $bonds[0]); # kb is now bonds[0]
						($bonds[1], $bonds[2]) = ($bonds[2], $bonds[1]); # r0 is now bonds[2]
						$bonds[1] = sqrt($bonds[1]/(2 * $bonds[0]));
						pop @bonds if ($#bonds > 2);
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
				$currParm->{$bond_counter} = (
												{
													"INDEX"		=> $bond_counter,
													"TYPE"		=> $3,
													"Lammps"	=> getLammpsOpts($3,"bond"),
													"VALS"		=> [@bonds],
													"USED"		=> 0,
													"KEY"		=> "$1 $2 ",
												}
											);

#				print "BOND $bond_counter: $key_code\n";
			}
		} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)/ and ($which_var == 5)) {
			if (uc($4) !~ /IGNORE|UFF_GEN/) {
				next if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
				next if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
				next if (! exists($PARMS{"ATOMTYPES"}{$3}) and lc($3) ne "x");
				
				@angles = getVals($5);
				next if (! @angles);
#				GetAbs(\@angles);
				if ($crossType eq "" && (! defined($alter) || $alter != 0)) {
					$angles[0] = $angles[0]/2; #same fix as bonds
				}
				if ($4 eq "COS_HARMON") {
					if (sin($angles[1] * pi/180) > 0.001) {
						$angles[0] /= sin($angles[1] * pi/180)**2;
					}
				} elsif ($4 eq "R-COSINE") { #bond angle cosine cross term
					$angles[3] = -$angles[3]/sin($angles[2] * pi/180);
					$angles[4] = -$angles[4]/sin($angles[2] * pi/180);
				} elsif ($4 eq "SW") {
					$angles[0] *= 2;
				}
				($atom1, $atom2, $atom3) = ($1, $2, $3);
				($atom1, $atom2, $atom3) = ($3, $2, $1) if (exists($PARMS{ANGLES}{$3}{$2}{$1}) and $4 ne "SW");
				$angle_counter = 1;
				if (exists($PARMS{ANGLES}{$atom1}{$atom2}{$atom3})) {
					$angle_counter = findDuplicate($PARMS{ANGLES}{$atom1}{$atom2}{$atom3}, $4);
					$angle_counter = scalar(keys %{ $PARMS{ANGLES}{$atom1}{$atom2}{$atom3} }) + 1 if (! $angle_counter);
				}
				$currParm = \%{ $PARMS{"ANGLES"}{$atom1}{$atom2}{$atom3} };
				$currParm->{$angle_counter} =  (
												{
													"INDEX"		=> $angle_counter,
													"TYPE"		=> $4,
													"Lammps"	=> getLammpsOpts($4,"angle", $CNV),
													"VALS"		=> [@angles],
													"USED"		=> 0,
													"KEY"		=> "$1 $2 $3 ",
													"CTYPE"		=> $crossType,
												}
												);
												#$currParm->{$angle_counter}{USED} = 1 if ($4 eq "SW");
#				print "ANGLE $angle_counter: $key_code\n";
			}
		} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)/ and ($which_var == 9)) {
			if (uc($4) !~ /IGNORE|UFF_GEN/) {
				next if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
				next if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
				next if (! exists($PARMS{"ATOMTYPES"}{$3}) and lc($3) ne "x");

				@angles = getVals($5);
#			   GetAbs(\@angles);
				if ($crossType eq "" && (! defined($alter) || $alter != 0)) {
					$angles[0] = $angles[0]/2; #same fix as bonds
				}
				($atom1, $atom2, $atom3) = ($1, $2, $3);
				($atom1, $atom2, $atom3) = ($3, $2, $1) if (! exists($PARMS{"ANGLES"}{$1}{$2}{$3}));
				next if (! keys %{ $PARMS{"ANGLES"}{$atom1}{$atom2}{$atom3} });
				for $i (values %{ $PARMS{"ANGLES"}{$atom1}{$atom2}{$atom3} }) {
					next if ($i->{Lammps}{name} !~ /harmonic/);
					$i->{Lammps}{name} = "charmm";
					push @{ $i->{VALS} }, (@angles);
					print "";
				}
#			   print "ANGLE $angle_counter: $key_code\n";
			}
		} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)/ and ($which_var == 6)) {
			$torsion_type = $5;
			$currParm = ();
			($atom1,$atom2,$atom3,$atom4) = ();
			if (uc($torsion_type) !~ /IGNORE|UFF_GEN/) {
				next if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
				next if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
				next if (! exists($PARMS{"ATOMTYPES"}{$3}) and lc($3) ne "x");
				next if (! exists($PARMS{"ATOMTYPES"}{$4}) and lc($4) ne "x");

				$torsion_counter++;
				@torsions = getVals($6);
#				GetAbs(\@torsions);
 
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
				($atom1, $atom2, $atom3, $atom4) = ($1, $2, $3, $4);
				($atom1, $atom2, $atom3, $atom4) = ($4, $3, $2, $1) if (exists($PARMS{TORSIONS}{$4}{$3}{$2}{$1}));
				$torsion_counter = 1;
				$torsion_counter = scalar(keys %{ $PARMS{"TORSIONS"}{$atom1}{$atom2}{$atom3}{$atom4} }) + 1
					if (exists($PARMS{"TORSIONS"}{$atom1}{$atom2}{$atom3}{$atom4}));
				$scale_flag = $PARMS{"PARMS"}{"scale_torsions"};
				$currParm = \%{ $PARMS{"TORSIONS"}{$atom1}{$atom2}{$atom3}{$atom4} };
				$currParm->{$torsion_counter} =  (
													{
														"TYPE"		=> $torsion_type,
														"Lammps"	=> getLammpsOpts($torsion_type,"dihedral", $CNV),
														"INDEX"		=> $torsion_counter,
														"VALS"		=> [@torsions],
														"USED"		=> 0,
														"KEY"		=> "$1 $2 $3 $4 ",
														"NUM"		=> 1,
														"PER"		=> ($#torsions + 1),
														"1_4scale"	=> $PARMS{"PARMS"}{"scale_vdw_14"},
														"CTYPE"		=> $crossType,
														"do_scale"	=> $scale_flag,
													}
												);
#				print "TORSION $torsion_counter: $key_code\n";
			}
		} elsif ($in_data =~ /^\s+(\-?\d+\.\d+)\s+(.+)/ and $which_var == 6) { #multiple torsions
			next if (! defined($atom1));
			@torsions = ();
			@torsions = getVals($2);
#			GetAbs(\@torsions);
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
			$currParm = \%{ $PARMS{"TORSIONS"}{$atom1}{$atom2}{$atom3}{$atom4} };
			push @{ $currParm->{$torsion_counter}{"VALS"} }, $tmp1;
			push @{ $currParm->{$torsion_counter}{"VALS"} }, @torsions;
			$currParm->{$torsion_counter}{NUM}++;
#			print "FOUND MULTIPLE TORSIONS FOR $key_code\n";
		} elsif ($in_data =~ /^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)/ and ($which_var == 7)) {
			$inversion_type = $5;
			if (uc($inversion_type) !~ /IGNORE|UFF_GEN/) {
				next if (! exists($PARMS{"ATOMTYPES"}{$1}) and lc($1) ne "x");
				next if (! exists($PARMS{"ATOMTYPES"}{$2}) and lc($2) ne "x");
				next if (! exists($PARMS{"ATOMTYPES"}{$3}) and lc($3) ne "x");
				next if (! exists($PARMS{"ATOMTYPES"}{$4}) and lc($4) ne "x");

				$inversion_counter++;
				@inversions = getVals($6);
				next if (! @inversions);
#				GetAbs(\@torsions);
 
				if ($crossType eq "" && (! defined($alter) || $alter != 0)) {
					if (lc($inversion_type) eq "it_jikl") {
						if ($inversions[1] == 180) {
							$inversions[1] = -1;
						} else {
							$inversions[1] = 1;
						}
					}
					$inversions[0] = $inversions[0] / 2 if ($inversion_type ne "UMBRELLA"); #1/2 fix again
				}
				$currParm = \%{ $PARMS{"INVERSIONS"}{$1}{$2}{$3}{$4} };
				$currParm->{$inversion_counter} = (
														{
															"TYPE"		=> $inversion_type,
															"Lammps"	=> getLammpsOpts($inversion_type,"inversion", $CNV),
															"INDEX"		=> $inversion_counter,
															"VALS"		=> [@inversions],
															"USED"		=> 0,
															"KEY"		=> "$1 $2 $3 $4 ",
															"CTYPE"		=> $crossType,
														}
												);
#				print "TORSION $torsion_counter: $key_code\n";
			}
		}			
		
	}
	close FORCEFIELD;

	die "ERROR: Force Field file $ff_file is invalid!\n"
		if (! %PARMS);
	
	if ($PARMS{"PARMS"}{"cut_vdw"} > $PARMS{"PARMS"}{"cut_coul"}) {
		$PARMS{"PARMS"}{"cut_vdw"} = $PARMS{"PARMS"}{"cut_coul"};
	}
	$PARMS{"PARMS"}{"type_counter"} = $type_counter;

	return (\%PARMS);
}

sub saveCeriusFF {
	my ($ffData, $save_name, $ELEMENTS) = @_;
	my ($c, $shft_angle, $torsions, $counter, $element, $vdws, $i);
	my ($atom1, $atom2, $atom3, $atom4, @ATMS, $hk, $inversions, @TMP);

	$vdws = "";
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
				$c->{"t0"} = ($c->{"t0"} * 180/pi);
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
						$c->{"p0"} = (($c->{"p0"} * 180/pi));# % 360);
						if ($counter == 0) {
							printf OUTFILE "%11.6f%10.4f%10.4f\n", ($c->{"kp"} * 2), $c->{"n"}, $c->{"p0"};
						} else {
							printf OUTFILE "%50s%11.6f%10.4f%10.4f\n", "", ($c->{"kp"} * 2), $c->{"n"}, $c->{"p0"};
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
	print OUTFILE "END\n\#\nCOULOMBIC\n X		X		   CONST-EPS\nEND\n";
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
	my ($ffs, $step, $alter) = @_;
	my ($i, $PARMS, $counter, $len, $printStr, $ffType);

	$alter = 1 if (! defined($alter));
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
			$PARMS = parseCerius2FF($i, $alter, $PARMS);
		} else {
			print "Loading $i->{FFTYPE} force field $i->{FF}...";
			$len = 43 + length($i->{FF});
			if ($i->{FFTYPE} eq "CERIUS2") {
				$PARMS = &parseCerius2FF($i->{FF},$alter, $PARMS);
			} else {
				$PARMS = &parseMPSimFF($i->{FF}, $alter, $PARMS);
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
	my ($FF, $reaxFF) = @_;
	my ($i, @FFILES, @tmp, $rec, $isCerius2, $ffType, $currFF, @tmp1);
	
	if ($FF =~ /\s+/) {
		@tmp = split /\s+/, $FF;
	} else {
		@tmp = ($FF);
	}
	$ffType = 0;
	for $i (@tmp) {
		$rec = ();
		$i = "${ff_dir}/${i}.ff" if (! -e $i and -e "${ff_dir}/${i}.ff");
		$i = "${ff_dir}/WAT/${i}.ff" if (! -e $i and -e "${ff_dir}/WAT/${i}.ff");
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
			$rec->{FF} = abs_path($i);
			if ($isCerius2) {
				$rec->{FFTYPE} = "CERIUS2";
			} else {
				$rec->{FFTYPE} = "MPSIM";
			}
			if ($#FFILES == -1) {
				$ffType = 1 if ($i =~ /amber|gaff/i);
				$ffType = 2 if ($i =~ /charmm/i);
				$ffType = 4 if ($i =~ /dreid/i);
				$ffType = 6 if ($i =~ /opls/i);
			}
		} elsif (uc($i) eq "GAFF") {
			$rec = (
						{
							"FF"		=> "${ff_dir}/GAFF.ff",
							"FFTYPE"	=> "CERIUS2",
						}
					);
			$ffType = 1 if (! defined($currFF) or $currFF eq "AMBER");
			$currFF = "AMBER";
		} elsif (uc($i) eq "AMBER03") {
			$rec = (
						{
							"FF"		=> "${ff_dir}/AMBER03.ff",
							"FFTYPE"	=> "CERIUS2",
						}
					);
			$ffType = 1;
		} elsif (uc($i) eq "AMBER03_SPC") {
			$rec = (
						{
							"FF"		=> "${ff_dir}/AMBER03_SPC.ff",
							"FFTYPE"	=> "CERIUS2",
						}
					);
			$ffType = 1 if (! defined($currFF) or $currFF eq "AMBER");
			$currFF = "AMBER";
		} elsif (uc($i) eq "AMBER96") {
			$rec = (
						{
							"FF"		=> "${ff_dir}/AMBER96.ff",
							"FFTYPE"	=> "CERIUS2",
						}
					);
			$ffType = 1 if (! defined($currFF) or $currFF eq "AMBER");
			$currFF = "AMBER";
		} elsif (uc($i) eq "AMBER99") {
			$rec = (
						{
							"FF"		=> "${ff_dir}/AMBER99.ff",
							"FFTYPE"	=> "CERIUS2",
						}
					);
			$ffType = 1 if (! defined($currFF) or $currFF eq "AMBER");
			$currFF = "AMBER";
		} elsif (uc($i) eq "AMBER_MBSC0") {
			$rec = (
						{
							"FF"		=> "${ff_dir}/AMBER_MBSC0.ff",
							"FFTYPE"	=> "CERIUS2",
						}
					);
			$ffType = 1 if (! defined($currFF) or $currFF eq "AMBER");
			$currFF = "AMBER";
		} elsif (uc($i) eq "AMBER91") {
			$rec = (
						{
							"FF"		=> "${ff_dir}/AMBER91.ff",
							"FFTYPE"	=> "CERIUS2",
						}
					);
			$ffType = 1 if (! defined($currFF) or $currFF eq "AMBER");
			$currFF = "AMBER";
		} elsif (uc($i) eq "CHARMM_LIPID") {
			$rec = (
						{
							"FF"		=> "${ff_dir}/charmm_par_all27_lipid.ff",
							"FFTYPE"	=> "CERIUS2",
						}
					);
			$ffType = 2 if (! defined($currFF) or $currFF eq "CHARMM");
			$currFF = "CHARMM";
		} elsif (uc($i) eq "CHARMM") {
			$rec = (
						{
							"FF"		=> "${ff_dir}/charmm_par_all27_prot_na.ff",
							"FFTYPE"	=> "CERIUS2",
						}
					);
			$ffType = 2 if (! defined($currFF) or $currFF eq "CHARMM");
			$currFF = "CHARMM";
		} elsif (uc($i)  eq "MESODNA") {
			$rec = (
						{
							"FF"		=> "${ff_dir}/ff/MesoDNA.ff",
							"FFTYPE"	=> "CERIUS2",
						}
					);
			$ffType = 3;
		} elsif (uc($i) eq "DREIDING") {
			$rec = (
						{
							"FF"		=> "${ff_dir}/DREIDING2.21.ff",
							"FFTYPE"	=> "CERIUS2",
						}
					);
			$ffType = 4;
		} elsif (uc($i) eq "REAX") {
			@tmp1 = split /\s+/, $reaxFF;
			$ffType = 5;
			for (@tmp1) {
				next if (! -e $_ or ! -r $_ or ! -T $_);
				$rec = (
						{
							"FF"		=> $_,
							"FFTYPE"	=> "REAX",
						}
					);
				push @FFILES, $rec;
				undef($rec);
			}
		}
		push @FFILES, $rec if ($rec);
	}
	die "ERROR: No valid files found!\n" if (! @FFILES);

	return (\@FFILES, $ffType);
}

sub loadElements {
	my (%ELEMENTS, $indata, $eleNum, $myDir);
	$myDir = "${scripts_dir}/Packages";
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

package Packages::Bulk;

use strict;
use constant q2C => 1.602176487E-19; # q = 1.602176487 x 10^-19 C
use constant a2m => 10E-10; # 1 A = 10^-10 m
use constant kb => 1.3806504E-23; # boltzmann constant m^2 kg s^-2 K^-1
use constant PI => atan2(1,1) * 4;
use constant qa2Cm => 1.602176487E-29; # q x A -> C x m
use constant debye => 2.9979E29; # 1 C m = 2.9979 x 10^29 D 
use constant e0 => 8.854187817E-12; # vacuum permitivity A^2 s^4 kg^-1 m^-3

use Packages::General qw(CrossProduct CenterText GetEquilPoint GetStats);
require Exporter;

our (@ISA, @EXPORT, $VERSION, @EXPORT_OK);
local $Carp::Carplevel = 1;
my ($cpack, $cfile) = caller();
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(CalcMolOrientCorrelation CalcDipole CalcDielectricConstant CalcDiffusionConstant SetOpts
		GetOrientOpt CalcHVap CalcHCap CalcTCompress CalcTExpand GetFluct WriteStats WriteData);
$VERSION = "1.00";

sub numerically { ($a<=>$b); }

sub acos {
    my($x) = $_[0];
    if (abs($x) > 1.0) {
        return 0;
    } else {
        return atan2(sqrt(1 - $x * $x), $x);
    }
}

sub dotProduct {
    my ($vec1, $vec2, $coord) = @_;
    my ($dM, $hKey, $i);

    $coord = "" if (! defined($coord));
    for $i ("X", "Y", "Z") {
        $hKey = "${i}${coord}";
        $dM += $vec1->{$hKey}*$vec2->{$hKey};
    }
    return $dM;
}

sub CalcMolOrientCorrelation {
    my ($atoms, $mols, $data, $orientOpts, $isVariable) = @_;
    my ($i, $count, $corr, $j, $curr, $start, $type, $a1, $a2, $dipole, $junk);
    my (@tmp, $offset, $molMap, %MOLORIENT, $v1, $v2, $tmpMol, $molData, $tot);

    $molMap->{0} = 0;
    $molMap->{"-1"} = -1;
    
    for $type (keys %{ $orientOpts }) {
	$count = $corr = 0;
	if ($type eq "T") {
	    for $i (keys %{ $mols }) {
		$tmpMol->{1} = $mols->{$i};
		if (! keys%{ $data->{START}{$i}{DIPOLE} }) {
		    $data->{START}{$i}{DIPOLE} = CalcDipole($atoms, $tmpMol, undef, -1);
		    $corr++;
		} else {
		    $curr = CalcDipole($atoms, $tmpMol, undef, -1);
		    $start = $data->{START}{$i}{DIPOLE};
		    $corr += dotProduct($curr, $start)/dotProduct($start, $start);
		}
		$count++;
		$tmpMol = ();
	    }
	} elsif ($orientOpts->{$type} =~ /^(\-?\d+)_(\-?\d+)_x_(\-?\d+)_(\-?\d+)$/) { #cross product
	    for $i (keys %{ $mols }) {
		@tmp = sort numerically keys %{ $mols->{$i} };
		for $j (0 .. $#tmp) {
		    $offset = $j + 1;
		    $molMap->{$offset} = $tmp[$j];
		}
		if ($1 == 0 and $2 == -1) {
		    $tmpMol->{1} = $mols->{$i};
		    $junk = CalcDipole($atoms, $tmpMol, undef, -1);
		    for $j ("X", "Y", "Z") {
			$v1->{"${j}COORD"} = $junk->{$j};
		    }
		} else {
		    $a1 = $atoms->{ $molMap->{$1} };
		    $a2 = $atoms->{ $molMap->{$2} };
		    $v1 = getVecs($a1, $a2);
		}

		if ($3 == 0 and $4 == -1) {
		    $tmpMol->{1} = $mols->{$i};
		    $junk = CalcDipole($atoms, $tmpMol, undef, -1);
		    for $j ("X", "Y", "Z") {
			$v2->{"${j}COORD"} = $junk->{$j};
		    }
		} else {
		    $a1 = $atoms->{ $molMap->{$3} };
		    $a2 = $atoms->{ $molMap->{$4} };
		    $v2 = getVecs($a1, $a2);
		}
		
		if (! keys %{ $data->{START}{$i}{$type} }) { # if this is the first time
		    $data->{START}{$i}{$type} = CrossProduct($v1, $v2);
		    $corr++;
		} else {
		    $start = $data->{START}{$i}{$type};
		    $curr = CrossProduct($v1, $v2);
		    $curr = CrossProduct($v1, $v2);	    
		    $corr += dotProduct($curr, $start, "COORD")/dotProduct($start, $start, "COORD");
		}
		$count++;				
	    }
	} elsif ($orientOpts->{$type} =~ /^(\d+)_(\d+)$/) { #regular vector
	    for $i (keys %{ $mols }) {
		@tmp = sort numerically keys %{ $mols->{$i} };
		for $j (0 .. $#tmp) {
		    $offset = $j + 1;
		    $molMap->{$offset} = $tmp[$j];
		}
		$a1 = $atoms->{ $molMap->{$1} };
		$a2 = $atoms->{ $molMap->{$2} };

		if (! keys %{ $data->{START}{$i}{$type} }) { # if this is the first time
		    $data->{START}{$i}{$type} = getVecs($a1, $a2);
		    $corr++;
		} else {
		    $start = $data->{START}{$i}{$type};
		    $curr = getVecs($a1, $a2);
		    $corr += dotProduct($curr, $start, "COORD")/dotProduct($start, $start, "COORD");
		}
		$count++;
	    }		
	}
	
	$corr /= $count;
	$MOLORIENT{$type} = $corr;
    }

    $tot = 0;
    if ($isVariable) { # remove starting values of all molecules not in current tstep
        for (keys %{ $data->{START} }) {
            if (! exists($mols->{$_})) {
                delete $data->{START}{$_};
                $tot++;
            }
        }
    }

    return \%MOLORIENT;
}

sub CalcDipole {
    my ($atoms, $mols, $data, $doMoments) = @_;
    my ($i, %DIPOLE, $dMoment, %MOL_DIPOLE, $j, @tmp, $k, $dp, $factor);

    @tmp = ("X", "Y", "Z");

    $factor = qa2Cm;
    $factor = 1 if ($doMoments == -1);
    for $j (@tmp) {
	$DIPOLE{$j} = 0;
	for $i (keys %{ $mols }) {
	    for $k (keys %{ $mols->{$i} }) {
		$dp = $factor * $atoms->{$k}{CHARGE} * $atoms->{$k}{"${j}COORD"};
		$DIPOLE{$j} += $dp;
		$MOL_DIPOLE{$i}{$j} += $dp;
	    }
	}
    }

    for $i (keys %MOL_DIPOLE) {
	$dMoment += sqrt($MOL_DIPOLE{$i}{X}**2 +  $MOL_DIPOLE{$i}{Y}**2 + $MOL_DIPOLE{$i}{Z}**2);
    }
    $dMoment /= (scalar keys %MOL_DIPOLE);
    $DIPOLE{u} = $dMoment;
    $DIPOLE{u2} = $dMoment * $dMoment;
    $DIPOLE{T} = $DIPOLE{X} + $DIPOLE{Y} + $DIPOLE{Z};
    if ($doMoments == 1) {
	$data->{DIPOLE_MOMENT} += debye * $dMoment;
	$dMoment = $data->{DIPOLE_MOMENT}/$data->{Count};
	return $dMoment;
    } else {
	return \%DIPOLE;
    }
}

sub CalcDielectricConstant {
    my ($atoms, $mols, $data, $volData, $tempData) = @_;
    my ($avgVol, $avgTemp, $factor, $i, $currMol, $fluct, $selfD); 
    my ($tot, $dipole, @list, $dipoleM, $dielectric, $dipoleA);

    $tot = $fluct = $dipoleM = 0;
    $data->{u} = $data->{G} = $data->{Count} = 0 if (! $data);
    $data->{Count}++;
    @list = keys %{ $mols };

    for $i (keys %{ $volData->{tStep} }) {
	$avgVol += $volData->{tStep}{$i};
        $avgTemp += $tempData->{tStep}{$i};
    }
    $avgVol *= a2m**3/$data->{Count}; # average volume in m^3
    $avgTemp /= $data->{Count}; # average temperature
    $factor = (4 * PI)/(3 * kb * e0 * $avgTemp * $avgVol); #4pi/3Ve0kbT - dimensionless since C = A s

    for $i (@list) {
	$tot++;
	$currMol->{1} = $mols->{$i};
	$dipole->{$i} = CalcDipole($atoms, $currMol, undef, 0); 
	$dipoleA->{X} += $dipole->{$i}{X};
        $dipoleA->{Y} += $dipole->{$i}{Y};
        $dipoleA->{Z} += $dipole->{$i}{Z};
    }

    for $i (@list) {
	$selfD = $dipole->{$i}{u2};
	$dipoleM += $selfD;
	$fluct += dotProduct($dipole->{$i}, $dipoleA);
    }
    $data->{u} += $dipoleM/$tot; #average dipole moment squared
    $data->{C} += $tot;
    $data->{G} += $fluct/(($data->{C}/$data->{Count}) * ($data->{u}/$data->{Count}));
    $dielectric->{G} = $data->{G}/$data->{Count};
    $dielectric->{T} = $factor * ($data->{C}/$data->{Count}) * ($data->{u}/$data->{Count}) * $dielectric->{G};
    $dielectric->{T} *= 50; #i can't justify where this factor comes from but seems to be necessary!
    return $dielectric;
}

sub CalcDiffusionConstant {
    my ($atoms, $mols, $data, $delT, $isVariable) = @_;
    my ($i, $j, $k, $d, $dist, $MOL, @tmp, $CENTER, $factor, $tot);

    $d = {
	"r" => 0,
	"T" => 0,
	"x" => 0,
	"y" => 0,
	"z" => 0,
         };
    $tot = 0;
    @tmp = ("X", "Y", "Z");
    for $k (keys %{ $mols }) {
	for $i (keys %{ $mols->{$k} }) {
	    if (! exists($data->{START}{$i})) {
		$data->{START}{$i} = {
					"XCOORD" => $atoms->{$i}{XCOORD}, 
					"YCOORD" => $atoms->{$i}{YCOORD}, 
					"ZCOORD" => $atoms->{$i}{ZCOORD}
				  };
		next;
	    }
	    $tot++;
	    for $j (@tmp) {
		$dist = $data->{START}{$i}{"${j}COORD"} - $atoms->{$i}->{"${j}COORD"};
		$d->{lc($j)} += $dist**2;
	    }
	}
    }
    return $d if (! $tot);

    $factor = 10/6; #convert msd to x 10E-9 m^2/s diffusion constant
    $d->{r} = $d->{x} + $d->{y} + $d->{z};
    for $i (keys %{ $d }) {
	$d->{$i} /= $tot;
    }
    $d->{T} = $d->{r} * $factor;
    $data->{T} = ($d->{T}) if (! exists($data->{T})); #store original displacement
    $d->{T} = ($d->{T} - $data->{T})/$delT; # approximation to slope

    $tot = 0;
    if ($isVariable) { # remove starting values of all molecules not in current tstep
        for (keys %{ $data->{START} }) {
	    if (! exists($atoms->{$_})) {
		delete $data->{START}{$_};
		$tot++;
	    }
        }
    }

    return $d;
}

sub getVecs {
    my ($a1, $a2) = @_;
    my (%vec, $i);

    for $i ("XCOORD", "YCOORD", "ZCOORD") {
        $vec{$i} = $a1->{$i} - $a2->{$i};
    }

    return \%vec;
}

sub GetOrientOpt {
    my ($atoms, $mols, $orientOpts) = @_;
    my (%OPTS, $i, $currOpt, $molMap, $offset, $opt1, $opt2, %OPTLIST, $optName, @tmp, @tmp2);

    @tmp = keys %{ $mols };
    @tmp2 = sort numerically keys %{ $mols->{ $tmp[0] } };

    for $i (@tmp2) {
        $offset = $i - $tmp2[0] + 1;
        $molMap->{$offset} = $i;
    }

    $orientOpts = lc($orientOpts);

    while ($orientOpts =~ /(\S+)/g) {
        $currOpt = $1;
        if ($currOpt =~ /dm/i) {
            $OPTS{T} = "dM";
            $OPTLIST{"dm"} = "0_-1";
	}
        if ($currOpt =~ /(\w+)_x_(\w+)\:(\w+)/) {
            ($opt1, $opt2, $optName) = ($1, $2, $3);
            next if (lc($opt1) eq lc($opt2));
            for ($opt1, $opt2) {
                if ($_ =~ /(\d+)_(\d+)/) {
                    next if ($1 == $2);
                    if (exists($molMap->{$1}) and exists($molMap->{$2}) and $1 != $2) {
                        $OPTLIST{"${1}_${2}"} = 1;
                    }
                }
            }
            if (exists($OPTLIST{$opt1}) and exists($OPTLIST{$opt2})) {
                $opt1 = $OPTLIST{$opt1} if ($OPTLIST{$opt1} ne "1");
                $opt2 = $OPTLIST{$opt2} if ($OPTLIST{$opt2} ne "1");
                $OPTS{$optName} = "${opt1}_x_${opt2}";
            }
        } elsif ($currOpt =~ /(\d+)_(\d+)\:(\w+)/) {
            next if ($1 == $2);
            $optName = $3;
            undef($opt1);
            if (exists($molMap->{$1}) and exists($molMap->{$2}) and $1 != $2) {
                $OPTLIST{"${1}_${2}"} = 1;
                $OPTLIST{$optName} = "${1}_${2}";
                $opt1 = "${1}_${2}";
            }
            if (defined($opt1) and exists($OPTLIST{$opt1})) {
                $OPTS{$optName} = "${1}_${2}";
            }
        }
    }
    return \%OPTS;
}

sub CalcTExpand {
    my ($STATS) = @_;
    my ($count, $factor, $thermalExpand, $results);

    $factor = 10E6/5.77668; # 1/kb; kb in kcal/mol
    $count = $STATS->{counter};

    $results->{vh} = $STATS->{VH}/$count;
    $results->{temp_2} = $STATS->{T_2}/$count;
    $results->{vol} = $STATS->{V}/$count;
    $results->{h} = $STATS->{H}/$count;

    $thermalExpand = ($results->{vh} - ($results->{vol}*$results->{h}))/($results->{temp_2} * $results->{vol});
    $results->{T} = $factor * $thermalExpand; # 10E-4 K-1
    return $results;
}

sub CalcTCompress {
    my ($STATS) = @_;
    my ($count, $factor, $thermalCompress, $results);

    $factor = 10E3/1.380; # V/kb; kb in kj/mol
    $count = $STATS->{counter};

    $results->{vol_2} = $STATS->{V_2}/$count;
    $results->{vol} = $STATS->{V}/$count;
    $results->{temp} = $STATS->{T}/$count;

    $thermalCompress = ($results->{vol_2} - $results->{vol}**2)/($results->{temp} * $results->{vol});
    $results->{T} = $factor * $thermalCompress; # 10E-6 atm-1
    return $results;
}

sub CalcHCap {
    my ($STATS) = @_;
    my ($count, $heatCapacity, $Gas_Constant, $results);

    $Gas_Constant = 1.987/1000; # kcal mol-1 K-1
    $count = $STATS->{counter};

    $results->{h_2} = $STATS->{H_2}/$count; #(kcal mol-1)**2
    $results->{h} = $STATS->{H}/$count;
    $results->{temp} = $STATS->{T}/$count;

    $heatCapacity = ($results->{h_2} - $results->{h}**2)/($Gas_Constant * $results->{temp}**2);
    $results->{T} = $heatCapacity; #kcal mol-1 K-1
    return $results;
}

sub CalcHVap {
    my ($STATS, $totMols) = @_;
    my ($hVap, $counter, $results);

    $counter = $STATS->{counter};

    $results->{U} = $STATS->{E_int}/$counter;
    $results->{vol} = $STATS->{V}/$counter;
    $results->{press} = $STATS->{P}/$counter;
    $results->{temp} = $STATS->{T}/$counter;
    $results->{pV} = ($results->{press} * $results->{vol} * 0.6023)/(9.8692 * 1000); #kJ/mol;
    $results->{rT} = 8.314472 * $results->{temp}/1000; # kJ/mol

    #$hVap = -1 * 1000 * $potEng/($pV);
    #$pV = $rt = 0 if (! $applyCorrection);
    $hVap = (-1*($results->{U}/$totMols) + $results->{rT} - $results->{pV});
    #$hVap = -1*$potEng/800;
    $results->{T} = $hVap;
    #printf "HVap: %.3f (kJ/mol)...", $hVap;
    return $results;
}

sub GetFluct {
    my ($data, $STATS) = @_;

    $STATS->{H} += $data->{poteng};
    $STATS->{H_2} += $data->{poteng}**2;
    $STATS->{E_int} += ($data->{poteng} - $data->{e_bond} - $data->{e_angle}); # kcal/mol
    #$rt += (8.314472 * $data->{$i}{Temp}/1000); # kJ/mol
    $STATS->{T} += $data->{temp};
    $STATS->{T_2} += $data->{temp}**2;
    $STATS->{P} += $data->{press};
    $STATS->{V} += $data->{volume};
    $STATS->{V_2} += $data->{volume}**2;
    $STATS->{VH} += $data->{volume}*$data->{poteng};
    #$STATS->{VH_2} += ($data->{volume}*$data->{poteng})**2;
    $STATS->{counter}++;
}

sub WriteStats {
    my ($stats, $saveName, $isMulti) = @_;
    my ($outStr, $i, $title, $j, $header, @tmp, @fields, $count);

    $outStr = sprintf("%-45s%21s%21s\n", "BULK PROPERTY",
                      CenterText("Calc",21),CenterText("Exp(300K 1atm)",21));

    for $i (1 .. 90) {
        $outStr .= "=";
    }
    $outStr .= "\n";

    $title = "ALL";
    if (defined($isMulti)) {
	delete $stats->{VOLUME};
	delete $stats->{TEMPERATURE};
	delete $stats->{MOL_CORRELATION};
	$outStr = sprintf("%-45s%14s", "BULK PROPERTY",CenterText("Exp**",14));
	@tmp = sort {$a cmp $b} keys %{ $stats };
	$count = 0;
	for $i (@tmp) {
	    $outStr .= sprintf("%14s",CenterText(uc($i),14));
	    $count++;
	}
	$outStr .= "\n";
	$count *= 14;
	$count += 56;
        for $i (1 .. $count) {
            $outStr .= "=";
        }
        $outStr .= "\n";

	@fields = keys %{ $stats->{$tmp[0]} };
	for $i (@fields) {
	    if ($stats->{$tmp[0]}{$i}{UNITS}) {
		$header = sprintf("%-45s", $stats->{$tmp[0]}{$i}{STATS}{T}{CONVERGED} . "${i} (" . $stats->{$tmp[0]}{$i}{UNITS} . ")");
	    } else {
		$header = sprintf("%-45s", $stats->{$tmp[0]}{$i}{STATS}{T}{CONVERGED} . $i);
	    }

	    $outStr .= sprintf("%-45s%8.3f %5.3f", $header, $stats->{$tmp[0]}{$i}{EXP_A}, $stats->{$tmp[0]}{$i}{EXP_D});
	    for $j (@tmp) {
		$outStr .= sprintf("%8.3f %5.3f", $stats->{$j}{$i}{STATS}{T}{AVG}, $stats->{$j}{$i}{STATS}{T}{STDEV});
	    }
	    $outStr .= "\n";
	}
    } else {
	$outStr = sprintf("%-45s%21s%21s\n", "BULK PROPERTY",
                      CenterText("Calc",21),CenterText("Exp**",21));

	for $i (1 .. 90) {
	    $outStr .= "=";
	}
        $outStr .= "\n";

        for $i (keys %{ $stats }) {
	    if ($stats->{$i}{UNITS}) {
		$header = sprintf("%-45s", $stats->{$i}{STATS}{T}{CONVERGED} . "${i} (" . $stats->{$i}{UNITS} . ")");
	    } else {
		$header = sprintf("%-45s", $stats->{$i}{STATS}{T}{CONVERGED} . $i);
	    }
	    $outStr .= $header .
	    sprintf("%8.3f +/- %-8.3f%8.3f +/- %-8.3f\n", $stats->{$i}{STATS}{T}{AVG},
                    $stats->{$i}{STATS}{T}{STDEV}, $stats->{$i}{EXP_A}, $stats->{$i}{EXP_D});
	}
    }
    $outStr .= "note: * indicates converged result not obtained\n\t** experimental values for water at 300K and 1atm\n";
    print $outStr;
    open OUTDATA, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
    print OUTDATA $outStr;
    close OUTDATA;
}

sub WriteData {
    my ($data, $savePrefix, $calcStats) = @_;
    my ($type, $fileName, @index, @headers, $startStats, $num_str);
    my ($i, $j, $k, @tmp, $avg_str, $stdev_str, $count, $sData, $title);

    
    delete $data->{VOLUME};
    delete $data->{TEMPERATURE};
    $calcStats = 0 if (! defined($calcStats) or $calcStats !~ /^(1|yes)$/i);
    $calcStats = 1 if ($calcStats =~ /^(1|yes)$/i);

    for $type (keys %{ $data }) {
        if (! keys %{ $data->{$type}{tStep} }) {
            delete $data->{$type};
            next;
        }
        $count = 0;
        @index = sort numerically keys %{ $data->{$type}{tStep} };
        @headers = ("T", grep {!/T/} keys %{ $data->{$type}{tStep}{$index[0]} });
        $fileName = $savePrefix . "_" . lc($type) . ".dat";
        open DATA, "> $fileName" or die "ERROR: Cannot create $fileName: $!\n";
        printf DATA "#%11s ", "TIME(ps)";
        for $j (@headers) {
            printf DATA "%12s ", $j;
        }
        print DATA "\n";
        for $i (@index) {
            $count++;
            printf DATA "%12.3f ", $i;
            for $j (@headers) {
                printf DATA "%12.5f ", $data->{$type}{tStep}{$i}{$j};
            }
            print DATA "\n";
        }
	close DATA;
	next;
        $avg_str = sprintf("#%11s ", "AVG");
        $stdev_str = sprintf("#%11s ", "STDEV");
        $num_str = sprintf("#%11s ", "num_pts");
        for $j (@headers) {
            $data->{$type}{STATS}{$j}{CONVERGED} = "";

	    if ($calcStats or $type ne "MOL_CORRELATION") {
		($k, $sData) = GetEquilPoint($data->{$type}{tStep}, 0.02, $j);
		if (! $k) {
		    $startStats = sprintf("%.0f", (0.66 * scalar(@index)));
		} else {
		    $startStats = 1;
		}
		$data->{$type}{STATS}{$j} = GetStats($sData, $startStats);
		$data->{$type}{STATS}{$j}{CONVERGED} = "*" if (! $k);
	    } else {
		$data->{$type}{STATS}{$j} = findFirstZero($data->{$type}{tStep});
	    }
            $avg_str .= sprintf("%12.5G ",$data->{$type}{STATS}{$j}{AVG} );
            $stdev_str .= sprintf("%12.8G ", $data->{$type}{STATS}{$j}{STDEV});
            $num_str .= sprintf("%11d%1s ", $data->{$type}{STATS}{$j}{NUM}, $data->{$type}{STATS}{$j}{CONVERGED});
        }
        print DATA "${avg_str}\n${stdev_str}\n${num_str}\n";
        close DATA;
    }
}

sub findFirstZero {
    my ($data) = @_;
    my ($i, $firstZero, $STATS, $pos);

    $STATS = {"AVG" => scalar(@{ $data }), "STDEV" => 0, "NUM" => 0, "CONVERGED" => "*"};
    for $i (0 .. ($#{ $data } - 1)) {
	if ($data->[$i] == 0 or ($data->[$i] * $data->[$i+1]) < 0) {
	    $firstZero = $data->[$i];
	    $pos = $i;
	    last;
	}
    }
    if (defined($firstZero)) {
	$STATS->{AVG} = $firstZero;
	$STATS->{NUM} = $#{ $data } - $pos;
	$STATS->{CONVERGED} = "";
    }

    return $STATS;
}

sub SetOpts {
    my ($BULK) = $_[0];
    $BULK->{DIFFUSION_CONSTANT}{UNITS} = "x 10-9 m2 s-1";
    $BULK->{DIFFUSION_CONSTANT}{EXP_A} = 2.23;
    $BULK->{DIFFUSION_CONSTANT}{EXP_D} = 0.1;

    $BULK->{DENSITY}{UNITS} = "g/cm3";
    $BULK->{DENSITY}{EXP_A} = 0.99953;
    $BULK->{DENSITY}{EXP_D} = 0.0001;

    $BULK->{DIELECTRIC_CONSTANT}{UNITS} = "";
    $BULK->{DIELECTRIC_CONSTANT}{EXP_A} = 78.46;
    $BULK->{DIELECTRIC_CONSTANT}{EXP_D} = 0.04;

    $BULK->{HEAT_VAPORIZATION}{UNITS} = "kcal/mol";
    $BULK->{HEAT_VAPORIZATION}{EXP_A} = 10.5176;
    $BULK->{HEAT_VAPORIZATION}{EXP_D} = 0.03;

    $BULK->{HEAT_CAPACITY}{UNITS} = "cal mol-1 K-1";
    $BULK->{HEAT_CAPACITY}{EXP_A} = 18.004;
    $BULK->{HEAT_CAPACITY}{EXP_D} = 0.0006;

    $BULK->{THERMAL_COMPRESSIBILITY}{UNITS} = "x 10E-6 atm-1";
    $BULK->{THERMAL_COMPRESSIBILITY}{EXP_A} = 45.86;
    $BULK->{THERMAL_COMPRESSIBILITY}{EXP_D} = 0.2;

    $BULK->{THERMAL_EXPANSION}{UNITS} = "x 10E-4 K-1";
    $BULK->{THERMAL_EXPANSION}{EXP_A} = 2.558;
    $BULK->{THERMAL_EXPANSION}{EXP_D} = 0.2;

    $BULK->{MOL_CORRELATION}{UNITS} = "ps";
    $BULK->{MOL_CORRELATION}{EXP_A} = 3.6;
    $BULK->{MOL_CORRELATION}{EXP_D} = 0.2;
}

1;

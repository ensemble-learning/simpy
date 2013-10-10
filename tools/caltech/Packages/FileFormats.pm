package Packages::FileFormats;

use strict;
use Storable qw(dclone);

use Packages::ManipAtoms qw(GetMols);
use Packages::General qw(Trim);
use File::Basename qw(basename);
use Math::Trig qw(acos);

require Exporter;

our (@ISA, @EXPORT, $VERSION, @EXPORT_OK);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(createBGF GetBGFFileInfo addElement addCon addHeader sortByRes createMOL2 GetMSIFileInfo
createHeaders addBoxToHeader GetPDBFileInfo GetMOL2FileInfo GetResAtoms GetBondList createPDB createPQR 
PrintCoord DeleteAtoms AddAtom UpdateBGF GetSystemCharge AddMass GetBGFAtoms GetMass AMBER2MOL2Types 
insertHeaderRemark InsertMol);
$VERSION = "1.00";

sub numerically { ($a<=>$b); }

sub has_field {
   my ($hstruct, $key) = @_;
   my ($found, $i);

   $found = 0;
   for $i (keys %{ $hstruct }) {
	if(exists($hstruct->{$i}{$key})) {
	    $found = 1;
	    last;
	}
    }

    return $found;
}

sub H2others {
   my ($H) = $_[0];
   my ($crystx,$a,$b,$c,$la,$lb,$lc,$alpha,$beta,$gamma);

   my $fac=180.0/acos(-1.0);
   $a->[0]=$H->[0][0]; $a->[1]=$H->[0][1]; $a->[2]=$H->[0][2];
   $b->[0]=0;          $b->[1]=$H->[1][1]; $b->[2]=$H->[1][2];
   $c->[0]=0;          $c->[1]=0;          $c->[2]=$H->[2][2];
   $la = sqrt($a->[0]*$a->[0] + $a->[1]*$a->[1] + $a->[2]*$a->[2]);
   $lb = sqrt($b->[1]*$b->[1] + $b->[2]*$b->[2]); 
   $lc = $c->[2];
   $alpha = $fac * acos( $b->[2] / $lb );
   $beta =  $fac * acos( $a->[2] / $la );
   $gamma = $fac * acos( ($a->[1]*$b->[1] + $a->[2]*$b->[2])/$la/$lb );

   $crystx = sprintf("CRYSTX %11.5f%11.5f%11.5f%11.5f%11.5f%11.5f",$la,$lb,$lc,$alpha,$beta,$gamma);
   return $crystx;
}

sub GetMSIFileInfo {
   my ($msiName, $saveHeaders) = @_;
   my (%ATOMS, %BONDS, @HEADERS, %AtomMap, $mode, $i, $j);  
   my ($rec, $index, $atom1, $atom2, $counter, %DATA);
   my ($bond_order, $disp, $cell, $tmp, $coord1, $coord2);

   $HEADERS[0] = "BGVER 322";
   $HEADERS[1] = "";
   open MSIFILE, $msiName or die "ERROR: Cannot open $msiName: $!\n";
   while (<MSIFILE>) {
	chomp;
	$_ =~ s/\#//;
	if ($_ =~ /Model/) {
	    $mode = 1; #headers
	} elsif ($_ =~ /^\s*\((\d+) Atom/) {
	    $mode = 2; #atoms
	    $rec = ();
            $counter = $1;
	} elsif ($_ =~ /^\s*\(\d+ Bond/) {
	    $mode = 3; #bonds
            undef($atom1);
            undef($atom2);
	    undef($bond_order);
	    undef($disp);
	} elsif ($_ =~ /^\s+\)$/) {
            if ($mode == 2 and $rec and defined($rec->{INDEX})) {
                $rec->{RESNAME} = "RES";
                $rec->{RESNUM} = 444;
		$rec->{CHARGE} = 0 if (! exists($rec->{CHARGE}));
		$rec->{LABEL} = "ATOM";
		$rec->{LONEPAIRS} = 0;
		$rec->{ATMNAME} = substr($rec->{FFTYPE},0,2) if (! exists($rec->{ATMNAME}));
                for $i (keys %{ $rec }) {
                    $DATA{$counter}{$i} = $rec->{$i};
                }
            }
            $rec = ();
            if ($mode == 3 and defined($atom1) and defined($atom2)) {
		$DATA{$atom1}{tmp}{$atom2}{BOND} = 1;
		$DATA{$atom2}{tmp}{$atom1}{BOND} = 1;
                $DATA{$atom1}{NUMBONDS}++;
                $DATA{$atom2}{NUMBONDS}++;
		$DATA{$atom1}{tmp}{$atom2}{ORDER} = $bond_order if(defined($bond_order));
		$DATA{$atom2}{tmp}{$atom1}{ORDER} = $bond_order if(defined($bond_order));
		if(defined($disp)) {
		    $DATA{$atom1}{tmp}{$atom2}{DISPX} = $disp->[0] if($disp->[0] != 0);
		    $DATA{$atom1}{tmp}{$atom2}{DISPY} = $disp->[1] if($disp->[1] != 0);
		    $DATA{$atom1}{tmp}{$atom2}{DISPZ} = $disp->[2] if($disp->[2] != 0);
		    $DATA{$atom2}{tmp}{$atom1}{DISPX} = $disp->[0] if($disp->[0] != 0);
		    $DATA{$atom2}{tmp}{$atom1}{DISPY} = $disp->[1] if($disp->[1] != 0);
		    $DATA{$atom2}{tmp}{$atom1}{DISPZ} = $disp->[2] if($disp->[2] != 0);
		}
            }
	} elsif ($mode == 1 and $_=~ /Label \"(.+)\"/) {
	    $HEADERS[1] = "DESCRP $1";
	} elsif ($mode == 1 and $_ =~ /A C \"Save-comment \d+\" \"(.+)/) {
	    push @HEADERS, "REMARK $1";
	} elsif ($mode == 1 and $_ =~ /PeriodicType (\d+)/) {
	    push @HEADERS, "PERIOD $1";
	} elsif ($mode == 1 and $_ =~ /SpaceGroup \"(.*)\"/) {
	    push @HEADERS, "SGNAME $1";
	} elsif ($mode == 1 and $_ =~ /(A.|B.|C.) \((\-?\d+\.?\d*E?\-?\d*) (\-?\d+\.?\d*E?\-?\d*) (\-?\d+\.?\d*E?\-?\d*)\)/) {
	    push @{ $cell }, [$2,$3,$4];
	} elsif ($mode == 2 and $_ =~ /FFType \"(.+)\"/) {
	    $rec->{FFTYPE} = $1;
	} elsif ($mode == 2 and $_ =~ /XYZ\s*\((\-?\d+\.?\d*E?\-?\d*)\s*(\-?\d+\.?\d*E?\-?\d*)\s*(\-?\d+\.?\d*E?\-?\d*)/) {
	    $rec->{XCOORD} = $1;
	    $rec->{YCOORD} = $2;
	    $rec->{ZCOORD} = $3;
	} elsif ($mode == 2 and $_ =~ /Id\s*(\d+)/) {
	    $index = $1;
	    $rec->{INDEX} = $1;
	} elsif ($mode == 2 and $_ =~ /Label\s*\"(.+)\"/) {
	    $rec->{ATMNAME} = $1;
	} elsif ($mode == 2 and $_ =~ /PerParent (\d+)/) {
	    $AtomMap{$counter} = $1;
	} elsif ($mode == 2 and $_ =~ /Charge (\-?\d+\.?\d*)/) {
	    $rec->{CHARGE} = $1;
	} elsif ($mode == 3 and $_ =~ /Order (\d+)/) {
	    $bond_order = $1;
	} elsif ($mode == 3 and $_ =~ /PerTrans \((\-?\d+) (\-?\d+) (\-?\d+)\)/) {
	    @{ $disp } = ($1,$2,$3);
	} elsif ($mode == 3 and $_ =~ /Atom1 (\d+)/) {
	    $atom1 = $1;
	} elsif ($mode == 3 and $_ =~ /Atom2 (\d+)/) {
	    $atom2 = $1;
	}
    }
    close MSIFILE;
    die "ERROR: $msiName does not contain any valid information!\n" if (! %DATA);

    if(defined($cell)) {
	$cell  = H2others($cell);
	push @HEADERS, $cell;
    }

    for $i (keys %DATA) {
	$atom1 = $i;
	$atom1 = $AtomMap{$i} if(exists($AtomMap{$i}));
	for $j (keys %{ $DATA{$i}{tmp} }) {
	    $atom2 = $j;
	    $atom2 = $AtomMap{$j} if(exists($AtomMap{$j}));
	    next if(!exists($AtomMap{$i}) and ! exists($AtomMap{$j}));
	    for (keys %{ $DATA{$i}{tmp}{$j} }) {
		$DATA{$atom1}{tmp}{$atom2}{$_} = $DATA{$i}{tmp}{$j}{$_};
	    }
	    delete $DATA{$i}{tmp}{$j};
	}
	delete $DATA{$i} if (exists($AtomMap{$i}));
    }

    for $i (keys %DATA) {
	$atom1 = $DATA{$i}{INDEX};
	$bond_order = has_field($DATA{$i}{tmp}, "ORDER");
	$disp->[0] = has_field($DATA{$i}{tmp}, "DISPX");
	$disp->[1] = has_field($DATA{$i}{tmp}, "DISPY");
	$disp->[2] = has_field($DATA{$i}{tmp}, "DISPZ");
	for $j (sort numerically keys %{ $DATA{$i}{tmp} }) {
	    $atom2 = $DATA{$j}{INDEX};
	    push @{ $BONDS{$atom1} }, $atom2;
	    if($bond_order) {
		push @{ $DATA{$i}{ORDER} }, $DATA{$i}{tmp}{$j}{ORDER} if(exists($DATA{$i}{tmp}{$j}{ORDER}));
		push @{ $DATA{$i}{ORDER} }, 1 if(!exists($DATA{$i}{tmp}{$j}{ORDER}));
	    }
	    if($disp->[0]) {
		push @{ $DATA{$i}{DISPX} }, $DATA{$i}{tmp}{$j}{DISPX} if(exists($DATA{$i}{tmp}{$j}{DISPX}));
		push @{ $DATA{$i}{DISPX} }, 0 if(!exists($DATA{$i}{tmp}{$j}{DISPX}));
	    }
	    if($disp->[1]) {
		push @{ $DATA{$i}{DISPY} }, $DATA{$i}{tmp}{$j}{DISPY} if(exists($DATA{$i}{tmp}{$j}{DISPY}));
		push @{ $DATA{$i}{DISPY} }, 0 if(!exists($DATA{$i}{tmp}{$j}{DISPY}));
	    }
	    if($disp->[2]) {
		push @{ $DATA{$i}{DISPZ} }, $DATA{$i}{tmp}{$j}{DISPZ} if(exists($DATA{$i}{tmp}{$j}{DISPZ}));
		push @{ $DATA{$i}{DISPZ} }, 0 if(!exists($DATA{$i}{tmp}{$j}{DISPZ}));
	    }
	}
	delete $DATA{$i}{tmp};
	for $j (keys %{ $DATA{$i} }) {
	    $ATOMS{$atom1}{$j} = $DATA{$i}{$j};
	}
    }
 
    for $atom1 (keys %ATOMS) {
        for $i ("X","Y","Z") {
	    next if (!exists($ATOMS{$atom1}{"DISP${i}"}));
	    for $j (0 .. $#{ $ATOMS{$atom1}{"DISP${i}"} }) {
	        $tmp = $ATOMS{$atom1}{"DISP${i}"}[$j];
		$tmp = -$tmp if ($tmp < 0);
		next if($tmp == 0);
		$coord1 = $ATOMS{$atom1}{"${i}COORD"};
		$coord2 = $ATOMS{ $BONDS{$atom1}[$j] }{"${i}COORD"};
		$ATOMS{$atom1}{"DISP${i}"}[$j] = $tmp if($coord1>=$coord2);
		$ATOMS{$atom1}{"DISP${i}"}[$j] = -$tmp if($coord1<$coord2);
	    }
	}
    }

    if (defined($saveHeaders) and $saveHeaders =~ /^1/) {
	return (\%ATOMS, \%BONDS, \@HEADERS);
    } else {
	return (\%ATOMS, \%BONDS);
    }
}
sub createBGF {
    my ($AtmData, $CON, $file_name, $saveRadii) = @_;
    my ($atom, $Atm_no, $bond, $label, $i, $j, $x, $y, $z, @order, $tmp, $fields, $val);

    if (! defined($saveRadii)) {
	$saveRadii = 0;
    }

    open OUTDATA, "> $file_name" or die "Cannot create $file_name: $!\n";
    if (exists($AtmData->{"HEADER"})) {
	for (@{ $AtmData->{"HEADER"} }) {
	    print OUTDATA $_ . "\n";
	}
	delete $AtmData->{"HEADER"};
    }
    print OUTDATA "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5,f10.5)\n";

    for $Atm_no (sort numerically keys %{ $AtmData }) {
	$atom = \%{ $AtmData->{$Atm_no} };
	if (! exists($atom->{"LABEL"})) {
	    $label = "ATOM";
	} else {
	    $label = $atom->{"LABEL"};
	}
        $atom->{"NUMBONDS"} = 0 if (! exists($atom->{"NUMBONDS"}));
	$atom->{"NUMBONDS"} = scalar(@{ $CON->{$Atm_no} }) if (defined($CON->{$Atm_no}));
        $atom->{"LONEPAIRS"} = 0 if (! exists($atom->{"LONEPARIS"}));
	#if (exists($atom->{ORDER})) {
	    #$atom->{"NUMBONDS"} = 0;
	    #for $i (@{ $atom->{ORDER} }) {
		#$atom->{"NUMBONDS"} += $i;
	    #}
	#}
        #if (length($atom->{"ATMNAME"}) > 4) {
	    #$atom->{"ATMNAME"} = substr($atom->{"ATMNAME"},0,4);
        #}
        for $i ("RESONANCE", "OCCUPANCY", "RADII") {
	    $atom->{$i} = 0 if (! exists($atom->{$i}));
	}
	$atom->{RESNAME} = substr($atom->{RESNAME},0,3) if (length($atom->{RESNAME}) > 3);
        $atom->{CHAIN} = "X" if (! defined($atom->{CHAIN}));
	($x, $y, $z) = ($atom->{"XCOORD"}, $atom->{"YCOORD"}, $atom->{"ZCOORD"});
	$x -= $atom->{OFFSET}{XCOORD} if (exists($atom->{OFFSET}{XCOORD}));
        $y -= $atom->{OFFSET}{YCOORD} if (exists($atom->{OFFSET}{YCOORD}));
        $z -= $atom->{OFFSET}{ZCOORD} if (exists($atom->{OFFSET}{ZCOORD}));

	printf OUTDATA "%-6s %5d %-5s %3s %1s %5d%10.5f%10.5f%10.5f %-5s%3d%2d %8.5f",
	$label, $Atm_no, $atom->{"ATMNAME"}, $atom->{"RESNAME"}, $atom->{CHAIN}, $atom->{"RESNUM"}, 
	$x, $y, $z, $atom->{"FFTYPE"}, $atom->{"NUMBONDS"}, $atom->{"LONEPAIRS"}, $atom->{"CHARGE"}; 
	if ($saveRadii) {
	    printf OUTDATA "%2d%4d%8.3f", $atom->{"RESONANCE"}, $atom->{"OCCUPANCY"}, $atom->{"RADII"};
	}
	print OUTDATA "\n";
    }

    #print OUTDATA "FORMAT CONECT (a6,12i6)\nFORMAT ORDER (a6,i6,13f6.3)\n";
    print OUTDATA "FORMAT CONECT (a6,12i6)\n";
    
    for $atom (sort numerically keys %{ $AtmData }) {
	$tmp = ();
	$fields = ();
	$fields->{CONECT} = 1;
	for $i (0 .. $#{ $CON->{$atom} }) {
	    $tmp->[$i]{CONECT} = sprintf("%6d",$CON->{$atom}[$i]);
	    for $j ("ORDER", "DISPX", "DISPY", "DISPZ") {
		if (exists($AtmData->{$atom}{$j})) {
		    $val = 0;
		    $val = $AtmData->{$atom}{$j}[$i] if($i <= $#{ $AtmData->{$atom}{$j} } and defined($AtmData->{$atom}{$j}[$i]));
		    $tmp->[$i]{$j} = sprintf("%6d",$val);
		    $fields->{$j} = 1;
		}
	    }
	}
	#sort hashes of array
	if (! $tmp) {
	    printf OUTDATA "%-6s%6d\n","CONECT",$atom;
	    next;
	}
	@order = sort{$a->{CONECT}<=>$b->{CONECT}} @{ $tmp };
	for $j ("CONECT", "ORDER", "DISPX", "DISPY", "DISPZ") {
	    next if(! exists($fields->{$j}));
	    printf OUTDATA "%-6s%6d", $j, $atom;
	    for $i (@order) {
		printf OUTDATA $i->{$j};
	    }
	    print OUTDATA "\n";
	}
	$tmp = ();
    }
    
    print OUTDATA "END\n";
    close OUTDATA;

}

sub GetBGFFileInfo {

    my ($in_file, $saveheaders, $verbose) = ($_[0],$_[1], $_[2]);
    my (%ATOMS, %BONDS, @CON, $counter, $in_data, $atm_patern, %BBOX, @ORDER, $nb);
    my ($calcBox, @HEADERS, $chain, $atomC, $numHetAtms, $fftype, $radii, $molid);
    my (@tmp, $i, $j, $k, $conStr, $atm_patern2);

    if (defined($saveheaders) and $saveheaders =~ /^1/) {
       $saveheaders = 1;
    } else {
       $saveheaders = 0;
    }

    $numHetAtms = 0;

    $atm_patern = '^ATOM\s*(\d+)\s+(.{5})\S?\s+(\S+)\s\w?\s*(\d+)\s*' . 
	'(\-?\d+\.\d{5})\s*(\-?\d+\.\d{5})\s*(\-?\d+\.\d{5})\s*(\S+)\s+' . 
	'(\d+)\s+(\d+)\s+(\-?\d+\.\d+)\s*(\d+\s+\d+\s+\d+\.\d+)?';
    $atm_patern2 = '^ATOM\s*(\d+)\s+(\S+)\s+(\S+)\s\w?\s*(\d+)\s*' . 
	'(\-?\d+\.\d{5})\s*(\-?\d+\.\d{5})\s*(\-?\d+\.\d{5})\s*(\S+)\s+' . 
	'(\d+)\s+(\d+)\s+(\-?\d+\.\d+)\s*(\d+\s+\d+\s+\d+\.\d+)?';
#    $atm_patern = '^ATOM\s*(\d+)\s+(.{5})\S?\s+(\S+)\s\w?\s*(\d+)' . 
#	'(.{10})(.{10})(.{10})\s*(\S+)\s+' . 
#	'(\d+)\s+(\d+)\s+(\-?\d+\.\d+)\s*(\d+\s+\d+\s+\d+\.\d+)?';
    open BGFFILE, $in_file or die "Cannot open BGFfile $in_file: $!\n";
    $molid = 1;
    while (<BGFFILE>) {
	chomp;
	$in_data = $_;
	$in_data =~ s/^HETATM/ATOM  /;
	if ($in_data =~ /^XTLGRF|BIOGRF|DESCRP|REMARK|REMARK|FORCEFIELD|PERIOD|AXES|SGNAME|CELLS/) {
	    push @HEADERS, $in_data;
	} elsif ($in_data =~ /$atm_patern2/ or $in_data =~ /$atm_patern/) {
	    $atomC = Trim($1);
	    $fftype = $8;
	    die "ERROR: Invalid BGF format for atom $atomC. Aborting.\nLine is $_\n" if (! defined($fftype));
	    $ATOMS{$atomC} = (
			    {
				"INDEX"       => $atomC,
				"ATMNAME"     => $2,
				"RESNAME"     => $3,
				"RESNUM"      => $4,
				"XCOORD"      => $5,
				"YCOORD"      => $6,
				"ZCOORD"      => $7,
				"FFTYPE"      => $fftype,
				"NUMBONDS"    => $9,
				"LONEPAIRS"   => $10,
				"CHARGE"      => $11,
				"MOLECULEID"  => $molid,
                             }
                           );
	    $radii = 0;
	    $BONDS{$atomC} = ();
	    if ($12) {
		($ATOMS{$atomC}{"OCCUPANCY"},
		$ATOMS{$atomC}{"RESONANCE"},
		$ATOMS{$atomC}{"RADII"}) = split /\s+/,$12;
		$radii = $ATOMS{$atomC}{"RADII"};
	    }
            if ($in_data =~ /^ATOM\s+\d+\s+\S+\s+\S+\s+(\w)\s+/) {
                $ATOMS{$atomC}{"CHAIN"} = $1;
            } else {
                $ATOMS{$atomC}{"CHAIN"} = "A";
            }
            if ($_ =~ /^HETATM/) {
                $ATOMS{$atomC}{"LABEL"} = "HETATM";
		$numHetAtms++;
            } else {
                $ATOMS{$atomC}{"LABEL"} = "ATOM";
            }

        } elsif ($in_data =~ /^CRYSTX\s+(.+)/) {
            print "Using Box Information from file..." if (! defined($verbose) || $verbose eq "1");
	    push @HEADERS, $in_data;
        
	} elsif ($in_data =~ /^CONECT(.+)$/) {
	    $i = $1;
	    @CON = ();
	    while ($i =~ /(.{6})/g) {
		$j = $1;
		next if ($j !~ /\S+/);
		$j =~ s/\s+//g;
		push @CON, $j;
	    }
	    #for $i (0 .. $#CON) {
		#$CON[$i] =~ s/\s+//g;
	    #}
	    $atomC = shift @CON;
	    next if (! exists($ATOMS{$atomC}));
	    #@CON = sort  @CON;
	    for $counter (@CON) {
		next if (! exists($ATOMS{$counter}));
	        push @{ $BONDS{$atomC} }, $counter;
	    }
	    $ATOMS{$atomC}{NUMBONDS} = scalar(@CON);
	} elsif ($_ =~ /^(ORDER|DISP.{1})\s+(.+)$/) {
	    @ORDER = split /\s+/, $2;
	    $atomC = shift @ORDER;
	    $nb = 0;
	    for $counter (@ORDER) {
		push @{ $ATOMS{$atomC}{$1} }, $counter;
		$nb += $counter;
	    }
	    $ATOMS{$atomC}{NUMBONDS} = $nb if ($_ =~ /ORDER/);
	} elsif ($_ =~ /^TER|ENDMDL/) {
	    $molid++;
	}
    }
    close BGFFILE;

    die "ERROR: $in_file does not contain any ATOM/CONNECTION information!\n"
    	if (! %ATOMS or (! %BONDS and $atomC > $numHetAtms));

    #&GetMols(\%ATOMS, \%BONDS);
    if (! $saveheaders) {
	return (\%ATOMS, \%BONDS);
    } else {
        return (\%ATOMS, \%BONDS, \@HEADERS);
    }
}

sub StoreExtrema {
    my ($x, $y, $z, $radii, $par) = @_;
    my ($counter, $tester);

    @{ $tester } = keys %{ $par };
    if ($#{ $tester } == -1) {
	for $counter ("X", "Y", "Z") {
	    $par->{$counter} = (
				{
				    "lo" => 99999.999,
				    "hi" => -99999.999,
				}
				);
	}
    } else {
	$par->{"X"}{"hi"} = ($x + $radii)
	    if (($x + $radii) > $par->{"X"}{"hi"});
	$par->{"X"}{"lo"} = ($x - $radii)
	    if (($x - $radii) < $par->{"X"}{"lo"});
	$par->{"Y"}{"hi"} = ($y + $radii)
	    if (($y + $radii) > $par->{"Y"}{"hi"});
	$par->{"Y"}{"lo"} = ($y - $radii)
	    if (($y - $radii) < $par->{"Y"}{"lo"});
	$par->{"Z"}{"hi"} = ($z + $radii)
	    if (($z + $radii) > $par->{"Z"}{"hi"});
	$par->{"Z"}{"lo"} = ($z - $radii)
	    if (($z - $radii) < $par->{"Z"}{"lo"});
    }
}

sub addHeader {
    my ($ATMS, $HDATA) = @_;

    for (@{ $HDATA }) {
	push @{ $ATMS->{"HEADER"} }, $_;
    }
}

sub addCon {
    my ($ATOMS, $CONNECTIONS) = @_;
    my ($atom, $bond);
    for $atom (keys %{ $CONNECTIONS }) {
	for $bond (@{ $CONNECTIONS->{$atom} }) {
	    $ATOMS->{$atom}{"CONNECT"}{$bond} = "";
	    $ATOMS->{$bond}{"CONNECT"}{$atom} = "";
	}
    }
}
    
sub addElement {
    my ($DATA, $field) = @_;
    my ($element, $hashKey, $eleChar);

    for $hashKey (keys %{ $DATA }) {
	if ($DATA->{$hashKey}{$field} =~ /([A-Za-z]+)/) {
	    $eleChar = $1;
	} else {
	    $eleChar = $DATA->{$hashKey}{$field};
	}
	
	$DATA->{$hashKey}{"ELEMENT"} = substr($eleChar, 0, 1);
	$DATA->{$hashKey}{"INDEX"} = $hashKey;
    }
}

sub createHeaders {
    my ($BBOX, $title) = @_;
    my (@HEADER, %BOX, $muser, $mdate);

    $muser = $ENV{USER};
    $mdate = scalar(localtime);

    # get the box information
    if ($BBOX) {
	%BOX = (
		"X" => {
		    "ANGLE" => 90,
		    "VAL"   => $BBOX->{XCOORD}{hi} - $BBOX->{XCOORD}{lo},
		},
		"Y" => {
		    "ANGLE" => 90,
		    "VAL"   => $BBOX->{YCOORD}{hi} - $BBOX->{YCOORD}{lo},
		},
		"Z" => {
		    "ANGLE" => 90,
		    "VAL"   => $BBOX->{ZCOORD}{hi} - $BBOX->{ZCOORD}{lo},
		},
		);
    }

    $HEADER[0] = "BIOGRF 200";
    $HEADER[1] = "DESCRP $title";
    $HEADER[2] = "REMARK Created by $muser on $mdate"; 
    $HEADER[3] = "FORCEFIELD AMBER";
    if (%BOX) {
	$HEADER[4] = "PERIOD 111";
	$HEADER[5] = "AXES   ZYX";
	$HEADER[6] = "SGNAME P 1                  1    1";
	$HEADER[7] = sprintf("CRYSTX %11.5f%11.5f%11.5f%11.5f%11.5f%11.5f", 
			     $BOX{"X"}{"VAL"}, $BOX{"Y"}{"VAL"}, $BOX{"Z"}{"VAL"},
			     $BOX{"X"}{"ANGLE"}, $BOX{"Y"}{"ANGLE"}, $BOX{"Z"}{"ANGLE"});
	$HEADER[8] = sprintf("CELLS %5d%5d%5d%5d%5d%5d", -1, 1, -1, 1, -1, 1);
    }

    return \@HEADER;
}

sub addBoxToHeader {
    my ($HEADER, $BOX) = @_;
    my ($rec, $i, $j, $inc, $angles,$ang);

    $rec  = "PERIOD 111\nAXES   ZYX\nSGNAME P 1                  1    1\nCRYSTX";
    for ("X", "Y", "Z") {
	$rec .= sprintf("%11.5f", ($BOX->{$_}{"hi"} - $BOX->{$_}{"lo"}));
	$ang = 90;
	$ang = $BOX->{$_}{angle} if(exists($BOX->{$_}{angle}));
	$angles .= sprintf("%11.5f",$ang);
    }
    $rec .=  "$angles\nCELLS    -1    1   -1    1   -1    1";

    $i = 0;
    while ($i <= $#{ $HEADER }) {
	$inc = 1;
	for $j ("PERIOD", "AXES", "SGNAME", "CRYSTX", "CELLS") {
	    if ($HEADER->[$i] =~ /^$j/) {
		splice @{ $HEADER }, $i, 1;
		$inc = 0;
	    }
	}
	$i++ if ($inc);
    }
    push @{ $HEADER }, $rec;
}

sub GetPDBFileInfo {

    my ($infile, $include_solvent, $res_nm, $atom_no, $hasRadii) = @_;
    my ($inText, %FData, $rec, $isvalid, %BONDS);
    $isvalid = 0;
    open INFILE, $infile or die "Cannot open $infile: $!\n";
    while (<INFILE>) {
	chomp;
	$inText = $_;
	if ($inText =~ /^(ATOM  |HETATM)(.{5})(.{6})\s*(\S+)\s+\w*\s+(\d+)\s+(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*(\-?\d*\.*\d*)\s*(\d*\.*\d*)/) {
	    $res_nm = $4;
	    $atom_no = $2;
	    $rec = (
		    {
			"HEADER"   => $1,
			"LABEL"    => $3,
			"ATMNAME"  => $3,
			"RES_NAME" => $res_nm,
			"RESNAME"  => $res_nm,
			"RES_ID"   => $5,
			"RESNUM"   => $5,
			"XCOORD"   => $6,
			"YCOORD"   => $7,
			"ZCOORD"   => $8,
			"INDEX"    => $2,
			"ATMNUM"   => $2,
		    }
		    );
            $atom_no =~ s/\s+//g;
            next if ($atom_no !~ /^\d+$/);
	    
	    $rec->{"CHARGE"} = $9 if ($9);
	    $rec->{"RADII"} = $10 if ($10);		
	    if ($res_nm =~ /WAT|NA|MG/i) {
                $rec->{"SOLUTE"} = 0;
		if ($include_solvent) {
		    $FData{$atom_no} = $rec;
		}
	    } else {
                $rec->{"SOLUTE"} = 1;
		$FData{$atom_no} = $rec;
	    }

	    $isvalid = 1;
	    if (defined($hasRadii) and $hasRadii eq "1") {
		if ($inText =~ /(\-?\d+\.\d+)\s*(\-?\d+\.\d+)\s*/g) {
		    $FData{$atom_no}{"CHARGE"} = $1;
		    $FData{$atom_no}{"RADII"} = $2;
		} else {
		    $isvalid = 0;
		}
	    }
	} elsif ($inText =~ /^CONECT(.{5})(.+)$/) {
	     $atom_no = $1+0;
	     $rec = $2;
	     $BONDS{$atom_no} = ();
	     while ($rec =~ /(.{5})/g) {
		push @{ $BONDS{$atom_no} }, $1;
	     }
	}
    }
    close INFILE;

    die "Invalid PDB file $infile\n"
	if (! $isvalid);

    return (\%FData, \%BONDS);
}

sub GetMOL2FileInfo {
    my ($in_file, $calcBox) = @_;
    my ($isAtm, $isBond, $atm_patern, $bond_patern, $isHeader); 
    my (%ATOMS, $description, $order, %CONN, $tmp);

    $isAtm = $isBond = $isHeader = 0;
    
    $atm_patern = '^\s*(\d+)\s+(\S+)\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)' . 
	'\s+(\-?\d+\.\d+)\s+(\S+)\s+(\d+)\s+(\S*)\s*(\-?\d+\.\d+)';
    $bond_patern = '^\s+\d+\s+(\d+)\s+(\d+)\s+(\w+)';

    open MOL2FILE, $in_file or die "Cannot open MOL2 file $in_file: $!\n";
    while (<MOL2FILE>) {
	chomp;
	if ($_ =~ /@<TRIPOS>MOLECULE/) {
	    $isHeader = 1;
	    $isAtm = $isBond = 0;
	} elsif ($_ =~ /^@<TRIPOS>ATOM/) {
	    $isAtm = 1;
	    $isBond = $isHeader = 0;
	} elsif ($_ =~ /^@<TRIPOS>BOND/) {
	    $isBond = 1;
	    $isAtm = $isHeader = 0;
	} elsif ($isHeader and $_ =~ /(\S+)/) {
	    $description .= "$1 ";
	} elsif ($_ =~ /$atm_patern/ and $isAtm) {
	    $ATOMS{$1} = (
			  {
			      "ATMNAME"     => $2,
			      "XCOORD"      => $3,
			      "YCOORD"      => $4,
			      "ZCOORD"      => $5,
			      "FFTYPE"      => $6,
			      "LONEPAIRS"   => 0,
			      "CHARGE"      => $9,
			      "RESNAME"     => $8,
			      "RESNUM"      => $7,
			      "INDEX"       => $1,
			  }
			  );
	    $tmp  = $1;
	    $ATOMS{$tmp}{RESNAME} =~ s/\d+$//;
	} elsif ($_ =~ /$bond_patern/ and $isBond) {
	    if($3 eq "ar") {
		$order = 2;
	    } else {
		$order = $3;
	    }
	    push @{ $CONN{$1} }, $2;
	    push @{ $CONN{$2} }, $1;
	    push @{ $ATOMS{$1}{"CONNECT"} }, $2;
	    push @{ $ATOMS{$2}{"CONNECT"} }, $1;
	    push @{ $ATOMS{$1}{"ORDER"} }, $order;
	    push @{ $ATOMS{$2}{"ORDER"} }, $order;
	
	    $ATOMS{$1}{"NUMBONDS"} = ($#{ $ATOMS{$1}{"CONNECT"} } + 1);
	    $ATOMS{$2}{"NUMBONDS"} = ($#{ $ATOMS{$2}{"CONNECT"} } + 1);
	}
    }
    close MOL2FILE;

    my ($molname) = basename($in_file);
    $molname =~ s/\.\w+$//;
    push @{ $ATOMS{"HEADER"} }, "BIOGRF 200";
    push @{ $ATOMS{"HEADER"} }, "DESCRP $molname";
    push @{ $ATOMS{"HEADER"} }, "REMARK $description";
    push @{ $ATOMS{"HEADER"} }, "REMARK BGF file created by mol2bgf";
    push @{ $ATOMS{"HEADER"} }, "FORCEFIELD DREIDING";
    return (\%ATOMS, \%CONN);
}

sub createMOL2 {
    my ($ATOMS, $BONDS, $fName, $saveResInfo) = @_;
    my ($atom, $bond, @atmList, $bList, $aList, $resID); 
    my ($resName, $RESLIST, $i, $j, $order, $title);

    $title = substr(basename($fName),0,3);
    open MOL2, "> $fName" || die "ERROR: Cannot create MOL2 file $fName: $!\n";
    $i = $j = 0;
    $saveResInfo = 0 if (! defined($saveResInfo));
    @atmList = sort numerically keys %{ $ATOMS };
    $bList = "";
    for $atom (@atmList) {
	$j++;
	if (! $saveResInfo) {
	    $resID = 1;
	    $resName = "<1>";
	} else {
	    $resID = $ATOMS->{$atom}{RESNUM};
	    $resName = $ATOMS->{$atom}{RESNAME};
	}
	$aList .= sprintf("%7d %-4s %13.4f%10.4f%10.4f %-4s%7d%4s%16.5f\n",
			  $atom, $ATOMS->{$atom}{ATMNAME}, $ATOMS->{$atom}{XCOORD},
			  $ATOMS->{$atom}{YCOORD}, $ATOMS->{$atom}{ZCOORD},
			  $ATOMS->{$atom}{FFTYPE}, $resID, $resName, $ATOMS->{$atom}{CHARGE});
	$RESLIST->{ $ATOMS->{$atom}{RESNAME} } = 1;
	for $bond (0 .. $#{ $BONDS->{$atom} }) {
	    next if ($atom > $BONDS->{$atom}[$bond]);
	    $order = 1;
	    if (exists($ATOMS->{$atom}{ORDER})) {
		$order = $ATOMS->{$atom}{ORDER}[$bond];
	    }
	    $i++;
	    $bList .= sprintf("%6d %4d %4d %-5s\n", $i, $atom, $BONDS->{$atom}[$bond], $order);
	}
    }
    print MOL2 "\@<TRIPOS>MOLECULE\n$title\n";
    print MOL2 "$j $i 0 0 0\n";
    print MOL2 "SMALL\nUSER_CHARGES\n\n\@<TRIPOS>ATOM\n";
    print MOL2 "${aList}@<TRIPOS>BOND\n${bList}";
    close MOL2;    
}

sub IsInteger {
    my ($testStr) = $_[0];
    my ($returnStr);

    $returnStr = 1;
    if  ($testStr !~ /^\d+$/) {
	$returnStr = 0;
    }	

    return $returnStr;
}

sub GetResAtoms {
    my ($atom, $BONDS, $bondList) = @_;
    my ($tmp, $bond, $i);
    
    $bondList->{$atom} = 1;
    for $bond (@{ $BONDS->{$atom} }) {
	next
	    if (exists($bondList->{$bond}));
	$bondList->{$bond} = 1;
	$tmp = GetResAtoms($bond, $BONDS, $bondList);
	for $i (keys %{ $tmp }) {
	    $bondList->{$i} = 1;
	}
    }
    
    return $bondList;
}

sub PrintCoord {
    my ($atom) = $_[0];
    my ($returnStr);
    
    $returnStr = "(";
    
    for ("XCOORD", "YCOORD", "ZCOORD") {
	$returnStr .= sprintf("%.3f, ", $atom->{$_});
    }
    
    $returnStr = substr($returnStr, 0, -2) . ")";
   
    return $returnStr;
}

sub DeleteAtoms {
    my ($atmList, $atoms, $bonds) = @_;
    my ($atomC, $bondC, $i);
    
    for $atomC (keys %{ $atmList }) {
	delete $atoms->{$atomC};
	delete $bonds->{$atomC};
	#for $i (keys %{ $bonds }) {
	#    for $bondC (0 .. $#{ $bonds->{$i} }) {
	#	if ($bonds->{$i}[$bondC] and $bonds->{$i}[$bondC] == $atomC) {
	#	    delete $bonds->{$i}[$bondC];
	#	}
	#    }
	#}
    }
}

sub AddAtom {
    my ($currAtom, $atomIndex, $ATOMS, $BONDS) = @_;
    
    if (exists($ATOMS->{$atomIndex})) {
	$atomIndex = getLastAtom($ATOMS) + 1;
    }
    
    %{ $ATOMS->{$atomIndex} } = %{ $currAtom };
    @{ $BONDS->{$atomIndex} } = (); 
}

sub InsertMol {
    my ($molAtoms, $molBonds, $atoms, $bonds, $offset) = @_;
    my ($i, $j);

    $offset = scalar(keys %{ $atoms }) if (exists($atoms->{$offset+1}));
	
    for $i (keys %{ $molAtoms }) {
	$atoms->{$i+$offset} = dclone($molAtoms->{$i});
	$bonds->{$i+$offset} = ();
	for $j (0 .. $#{ $molBonds->{$i} }) {
	    $bonds->{$i+$offset}[$j] = ($molBonds->{$i}[$j]+$offset);
	}
    }
}

sub getLastAtom {
    my ($ATOMS) = $_[0];
    my (@atmList) = sort numerically keys %{ $ATOMS };
    
    return $atmList[$#atmList];
}

sub UpdateBGF {
    my ($ATOMS, $BONDS) = @_;
    my (%atmList, %bondList, $index, $atomC, $bondC);
    
    $index = 1;
    for $atomC (sort numerically keys %{ $ATOMS }) {
	%{ $atmList{$index} } = %{ $ATOMS->{$atomC} };
	$atmList{$index}{"INDEX"} = $atomC;
	$ATOMS->{$atomC}{"INDEX"} = $index;
	$index++;
    }
    
    for $atomC (keys %{ $BONDS }) {
	for $bondC (@{ $BONDS->{$atomC} }) {
            if (exists($ATOMS->{$bondC})) {
                push @{ $bondList{ $ATOMS->{$atomC}{"INDEX"} } }, $ATOMS->{$bondC}{"INDEX"};
            }
	}
    }
    
    $ATOMS = ();
    $BONDS = ();
    
    return (\%atmList, \%bondList);
}

sub GetSystemCharge {
    my ($ATOMS) = $_[0];
    my ($atom, $totCharge);

    $totCharge = 0;
    for $atom (keys %{ $ATOMS }) {
	$totCharge += $ATOMS->{$atom}{"CHARGE"};
    }
    
    return $totCharge;
}


sub AddMass {
    my ($atoms, $parms) = @_;
    my ($i, $fftype, $mass);
    
    for $i (keys %{ $atoms }) {
        $fftype = $atoms->{$i}{"FFTYPE"};
        die "ERROR: Atom $atoms->{$i}{ATMNAME} does not have a force field type\n"
            if (! defined($fftype));
        $mass = $parms->{"ATOMTYPES"}{$fftype}{"MASS"};
        die "ERROR: Atom $atoms->{$i}{ATMNAME} has an invalid force field type $fftype\n"
            if (! defined($mass));
        $atoms->{$i}{"MASS"} = $mass;
    }
}

sub GetMass {
    my ($parms) = $_[0];
    my ($i, $MASSES);

    for $i (keys %{ $parms->{"ATOMTYPES"} }) {
        $MASSES->{$i} = $parms->{"ATOMTYPES"}{$i}{"MASS"};
    }
    return $MASSES;
}

sub GetBGFAtoms {
    my ($atmList, $atoms, $bonds) = @_;
    my ($field, $index, $atomC, %BGF, %CONS, %TMP, $i);

    $atomC = 1;

    for $index (sort numerically keys %{ $atmList }) {
        next if (! exists($atoms->{$index}) or ! $atoms->{$index});
	$BGF{$atomC} = \%{ $atoms->{$index} };
	$atoms->{$index}{"FINDEX"} = $atomC;
	$BGF{$atomC}{"OINDEX"} = $index;
	$atomC++;
    }

    for $atomC (keys %BGF) {
	@{ $CONS{$atomC} } = ();
	for $index (@{ $bonds->{$BGF{$atomC}{"OINDEX"} } }) {
	    $i = $atoms->{$index}{"FINDEX"};
	    if ($i and exists($BGF{$i})) {
		push @{ $CONS{$atomC} }, $i;
	    }
	}
    }

    return (\%BGF, \%CONS, $atmList);

}

sub getAtmList {
    my ($atoms, $selectList, $field) = @_;
    my ($i, $atom, %sorted, %atmList);

    for $i (keys %{ $atoms }) {
	$sorted{ uc($atoms->{$i}{$field}) }{$i} = 1;
    }

    for $i (keys %{ $selectList }) {
	$i = uc($i);
	if (exists($sorted{$i})) {
	    for $atom (keys %{ $sorted{$i} }) {
		$atmList{$atom} = 1;
	    }
	}
    }

    return \%atmList;
}

sub sortByRes {
    my ($ATOMS, $select) = @_;
    my ($i, %RES, $resNum);

    %{ $select } = %{ $ATOMS } if (! defined($select));
    for $i (keys %{ $ATOMS }) {
	next if (! exists($select->{$i}));
	$resNum = $ATOMS->{$i}{RESNUM};
	$RES{$resNum}{ATOMS}{$i} = 1;
    }

    return \%RES;
}

sub GetBondList {
    my ($ATOMS, $BONDS, $SELECT) = @_;
    my ($atomC, $atom, $resAtms, $resC, %RESINFO, $i);

    $resC = 1;
    for $atomC (keys %{ $ATOMS }) {
	$atom = \%{ $ATOMS->{$atomC} };
	next if (keys %{ $SELECT } and ! exists($SELECT->{$atomC}));
	next if (exists($atom->{"bondpointer"}));
	$resAtms = GetResAtoms($atomC, $BONDS, ());
	$atom->{"bondpointer"} = $resC;
	$RESINFO{$resC}{$atomC} = 1;
	for $i (keys %{ $resAtms }) {
	    next if (keys %{ $SELECT } and ! exists($SELECT->{$i}));
	    $ATOMS->{$i}{"bondpointer"} = $resC;
	    $RESINFO{$resC}{$i} = 1;
	}
	$resC++;
    }

    return \%RESINFO;
}

sub createPDB {
    my ($atoms, $bonds, $saveName, $terList) = @_;
    my ($atmName, $resName, $i, $bondList, $j);
    my ($fmt) = '%-6s%5d%4s  %3s  %4d %11.3f%8.3f%8.3f';
    open PDBFILE, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
    for $i (sort numerically keys %{ $atoms }) {
        $atmName = $atoms->{$i}{"ATMNAME"};
        $resName = $atoms->{$i}{"RESNAME"};
	if (length($atmName) > 4) {
	    $atmName =~ /(\S+)/;
	    $atmName = $1;
	    $atmName =~ /^(.{4})/;
	    $atmName = $1;
	}
        $resName = substr($resName, 0,2) . substr($resName,-1,1) if (length($resName) > 3);
	print PDBFILE "TER\n" if(defined($terList) and exists($terList->{$i}));
        printf PDBFILE "$fmt\n", "ATOM", $i, $atmName, $resName, $atoms->{$i}{RESNUM},
        $atoms->{$i}{XCOORD}, $atoms->{$i}{YCOORD}, $atoms->{$i}{ZCOORD};
	$bondList .= sprintf("%-6s%5d", "CONECT", $i);
	for $j (@{ $bonds->{$i} }) {
	    $bondList .= sprintf("%5d",$j);
	}
	$bondList .= "\n";
    }

    print PDBFILE $bondList if (keys %{ $bonds });

    close PDBFILE;
}

sub createPQR {
    my ($atoms, $bonds, $saveName,$headers) = @_;
    my ($atmName, $resName, $i, $bondList, $j, @tmp);
    my ($fmt) = '%-6s%5d %4s %3s %5d %11.3f%8.3f%8.3f%8.3f%8.3f';

    open PQRFILE, "> $saveName" or die "ERROR: Cannot create $saveName: $!\n";
    if(defined($headers)) {
	for $i (@{ $headers }) {
	    if ($i =~ /CRYSTX\s+(.+)$/) {
		@tmp = split /\s+/, $1;
		printf PQRFILE "CRTST1 %8.3f %8.3f %8.3f %6.2f %6.2f %6.2f P 1           1\n",
		    $tmp[0],$tmp[1],$tmp[2],$tmp[3],$tmp[4],$tmp[5];
		last;
	    }
	}
    }
    for $i (sort numerically keys %{ $atoms }) {
	die "ERROR: Radii for atom $i not found.. Aborting\n" if (! exists($atoms->{$i}{RADII}));
        $atmName = $atoms->{$i}{"ATMNAME"};
        $resName = $atoms->{$i}{"RESNAME"};
	$atmName =~ s/\s//g if (length($atmName) > 4);
        if (length($atmName) > 4) {
            $atmName =~ /(.{4})$/;
            $atmName = $1;
        }
        $resName = substr($resName, 0,2) . substr($resName,-1,1) if (length($resName) > 3);
        printf PQRFILE "$fmt\n", "ATOM", $i, $atmName, $resName, $atoms->{$i}{RESNUM},
        $atoms->{$i}{XCOORD}, $atoms->{$i}{YCOORD}, $atoms->{$i}{ZCOORD}, 
	$atoms->{$i}{CHARGE}, $atoms->{$i}{RADII};
        $bondList .= sprintf("%-6s%5d", "CONECT", $i);
        for $j (@{ $bonds->{$i} }) {
            $bondList .= sprintf("%5d",$j);
        }
        $bondList .= "\n";
    }

    print PQRFILE "CONECT\n$bondList" if (keys %{ $bonds });

    close PQRFILE;
}

sub AMBER2MOL2Types {
    my ($atoms) = $_[0];
    my ($i, %ffmap, $ffType, $convertFile);
    
    $convertFile = "/ul/tpascal/scripts/dat/amber2mol2Types.perldata";
    die "ERROR: Cannot read $convertFile: $!\n" 
	if (! -e $convertFile or ! -r $convertFile);
    scalar eval `cat $convertFile`;
    for $i (keys %{ $atoms }) {
	$ffType = uc($atoms->{$i}{FFTYPE});
	die "ERROR: atom $i ($atoms->{$i}{ATMNAME}) has a atomtype $ffType with no match!\n" 
	    if (! exists($ffmap{$ffType}));
	$atoms->{$i}{FFTYPE} = $ffmap{$ffType};
    }
}

sub insertHeaderRemark {
    my ($headers, $insertRmk) = @_;
    my ($i, $end, $start);

    $headers = createHeaders(undef,"default.bgf") if (! $headers);
    $start = $end = 0;
    for $i (0 .. $#{ $headers }) {
	if ($headers->[$i] =~ /^DESCRP/) {
	   $start = $i;
	} elsif ($headers->[$i] =~ /^REMARK/) {
	   $start = $i - 1 if (! $start);
	   $end = $i;
	}
    }

    return () if (! $start and ! $end);
    $end++;
    splice @{ $headers }, $end, 0, $insertRmk;
}

1;

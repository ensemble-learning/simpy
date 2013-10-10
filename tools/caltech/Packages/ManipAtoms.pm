package Packages::ManipAtoms;
require Exporter;
use strict;
use Storable qw(dclone);
use Cwd;
use Packages::General qw(CoM GetBondLength);
no warnings 'recursion';

our (@ISA, @EXPORT, @EXPORT_OK, $VERSION);

@ISA = qw(Exporter);
$VERSION = "1.00";
@EXPORT = ();
@EXPORT_OK = qw(ImageAtoms ScaleAtoms UnwrapAtoms FindElement CenterSystem GetAtmList 
		RotateAbout GetSolvent SplitAtomsByMol GetAtmData GetMols MapOrigin 
		GetAbituraryRotationMatrix TransAtoms GroupAtomsByField MakeFieldSequential
		BuildAtomSelectionString SelectAtoms ReimageAtoms);

sub numerically { ($a<=>$b); }

sub MapOrigin {
    my ($box, $com) = @_;
    my ($boxCenter, $offset);

    for (keys %{ $box }) {
        $boxCenter->{$_} = ($box->{$_}{hi} - $box->{$_}{lo})/2;
        $offset->{$_} = $boxCenter->{$_} - $com->{$_ . "COORD"};
        $box->{$_}{hi} -= $offset->{$_};
        $box->{$_}{lo} -= $offset->{$_};
    }
}

sub CenterSystem {
    my ($atoms, $offset, $atomList) = @_;
    my ($i, $dim);

    for $i (keys %{ $atoms }) {
	next if (keys %{ $atomList } && ! exists($atomList->{$i}));
	for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	    $offset->{$dim} = $offset->{substr($dim,0,1)}{lo} if (! exists($offset->{$dim}));
	    next if (! $offset->{$dim});
	    $atoms->{$i}{$dim} -= $offset->{$dim};
	}
    }
}

sub FindElement {
    my ($fftype, $parms, $elements) = @_;
    my ($elementName, $elementNum, $i);

    die "ERROR: $fftype is not a valid force field type!\n" if (! exists($parms->{$fftype}));
    $elementName = $parms->{$fftype}{ATOM};
    $elementNum = 0;

    for $i (keys %{ $elements }) {
        if (uc($elements->{$i}{SYMBOL}) eq uc($elementName)) {
            $elementNum = $i;
            last;
        }
    }

    #die "ERROR: Element $elementName is not a valid element!\n" if (! $elementNum);

    return ($elementName, $elementNum);
}

sub ImageAtoms {
    my ($atoms, $CoM, $box) = @_;
    my ($MOL, $i, $OFFSET, $dim, $index, $CENTER);

    %{ $CENTER } = %{ $CoM };
    for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	$OFFSET->{$dim} = 0;
	while ($CENTER->{$dim} > $box->{$dim}{hi}) {
	    $OFFSET->{$dim}++;
	    $CENTER->{$dim} -= $box->{$dim}{len};
        }
	while ($CENTER->{$dim} < $box->{$dim}{lo}) {
	    $OFFSET->{$dim}--;
	    $CENTER->{$dim} += $box->{$dim}{len};
	}
    }
    
    for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	$index = $dim;
	$index =~ s/COORD/INDEX/;
	for $i (keys %{ $atoms }) {
	    $atoms->{$i}{$dim} -= ($OFFSET->{$dim} * $box->{$dim}{len});
	    $atoms->{$i}{$index} = $OFFSET->{$dim};
	    if ($atoms->{$i}{$dim} < $box->{$dim}{lo}) {
		$atoms->{$i}{$index}--;
	    } elsif ($atoms->{$i}{$dim} > $box->{$dim}{hi}) {
		$atoms->{$i}{$index}++;
	    }
	    $atoms->{$i}{IMAGED} = 1;
	}
    }

    return $CENTER;
}

sub ScaleAtoms {
    my ($atoms, $box) = @_;
    my ($i, $coord, $dim, $index, $pos);

    for $dim ("XCOORD", "YCOORD", "ZCOORD") {
	$index = $dim;
	$index =~ s/COORD/INDEX/;
	for $i (keys %{ $atoms }) {
	    $pos = $atoms->{$i}{$dim};
	    $pos /= $box->{$dim}{len};
	    if ($pos =~ /(\d+)(\.\d+)/) {
		if ($1 > 0) {
		    $atoms->{$i}{$index} = $1 + 1;
		} else {
		    $atoms->{$i}{$index} = $1;
		}
		$atoms->{$i}{$dim} = sprintf("%.6f ", $pos);
	    }
	}
    }
}

sub UnwrapAtoms {
    my ($ATOMS, $BOX, $scaled) = @_;
    my ($atomC, $atom, $coord, $index, $pos);

    for $atomC (keys %{ $ATOMS }) {
	$atom = \%{ $ATOMS->{$atomC} };
	for $coord ("XCOORD", "YCOORD", "ZCOORD") {
	    $index = $coord;
	    $index =~ s/COORD/INDEX/;
	    next if(!exists($atom->{$index}));
	    $index = $atom->{$index};
	    if (! defined($scaled) or $scaled == 1) {
	        $atom->{$coord} *= $BOX->{$coord}{len};
	        $atom->{$coord} += $BOX->{$coord}{lo};
	    }
	    $atom->{$coord} += ($index * $BOX->{$coord}{len});
	}
    }
}

sub GetAtmList {
    my ($select, $ATOMS) = @_;
    my ($i, %LIST, $field, $val, $rec, $operator, $excluded, $tmp);
    my ($count);

    for $field (keys %{ $select }) {
	$tmp = "";
	for $val (keys %{ $select->{$field} }) {
	    $excluded = 0;
	    $operator = "";
	    if ($val eq "*") {
		$operator = ">";
		$val = "0 and \$ATOMS->{\$i}{INDEX} < 999999";
	    } elsif ($val =~ /^\d+/) { #integer
		$operator = "==";
	    } elsif ($val =~ /^\!(\d+)/) {
		$operator = "!=";
		$val = $1;
		$excluded = 1;
	    } elsif ($val =~ /^(>|<|=|\+|\-)(\d+\.*\d*)/) {
		$operator = $1;
		$val = $2;
	    } elsif ($val =~ /^\!(>|<|=)(\d+\.*\d*)/) {
		$excluded = 1;
		if ($1 eq ">") {
		    $operator = "<";
		} elsif ($1 eq "<") {
		    $operator = ">";
		} else {
		    $operator = "!=";
		}
		$val = $2;
	    } elsif ($val =~ /^(\w+)$/) {
		$operator = "eq";
		$val = "\"${1}\"";
	    } elsif ($val =~ /^\!(\w+)$/) {
		$excluded = 1;
		$operator = "ne";
		$val = "\"${1}\"";
	    } elsif ($val =~ /^(\S+)/) {
		$operator = "=~";
		$val = "/$1/";
	    }
	    next if (! $operator);
	    if ($tmp) {
		if (! $excluded) {
		    $tmp .= "or \$ATOMS->{\$i}{${field}} $operator $val ";
		} else {
		    $tmp .= "and \$ATOMS->{\$i}{${field}} $operator $val ";
		}
	    } else {
		$tmp = "{${field}} $operator $val ";
	    }
	}
        if ($rec) {
            $rec .= "and (\$ATOMS->{\$i}$tmp) ";
        } else {
            $rec = "$tmp";
        }
    }

    $count = 0;
    for $i (keys %{ $ATOMS }) {
	$excluded = 0;
	if (! eval('$ATOMS->{$i}' . $rec)) {
	    $excluded = 1;
	}
	if (! $excluded) {
	    $count++;
	    $LIST{$i} = $count;
	}
    }
    
    return \%LIST;
}

sub GetSolvent {
    my ($atoms, $solvType) = @_;
    my (%SOLVOPTS, %SOLVATMS, $i);

    %SOLVOPTS = (
			"WATER" => {
					"MOLSIZE" => 3,
					"RESNAME" => "WAT|HOH|RES",
					"FFTYPE"  => "OW|HW|OT|HT|HO|OH|H_|O_3"
				    }
		);

    return () if (! defined($SOLVOPTS{uc($solvType)}));

    for $i (keys %{ $atoms }) {
	if ($atoms->{$i}{MOLSIZE} == $SOLVOPTS{$solvType}{MOLSIZE}) {
	    if ($atoms->{$i}{RESNAME} =~ /$SOLVOPTS{$solvType}{RESNAME}/ or
		$atoms->{$i}{FFTYPE} =~ /$SOLVOPTS{$solvType}{FFTYPE}/) {
		    $atoms->{$i}{IS_SOLVENT} = 1;
		    $SOLVATMS{$i} = 1;
	    }
	}
    }

    return \%SOLVATMS;
}

sub getMolList {
    my ($atoms, $bonds, $atomID, $mol, $used) = @_;
    my ($i);

    $mol->{MEMBERS}{$atomID} = 1;
    $mol->{MOLSIZE}++;
    for $i ("MOLECULE", "MOLECULEID", "MOLSIZE") {
	delete $atoms->{$atomID}{$i};
    }
    $atoms->{$atomID}{MOLECULE} = $mol;
    $atoms->{$atomID}{MOLECULEID} = \$atoms->{$atomID}{MOLECULE}{INDEX};
    $atoms->{$atomID}{MOLSIZE} = \$atoms->{$atomID}{MOLECULE}{MOLSIZE};
    $used->{$atomID} = 1;
    for $i (@{ $bonds->{$atomID} }) {
	next if (exists($used->{$i}));
	$used->{$i} = 1;
	&getMolList($atoms, $bonds, $i, $mol, $used);
     }
}

sub GetMols {
    my ($atoms, $bonds, $select) = @_;
    my ($i, $j, $MOLS, $counter, $msize, $USED);

    $select = $atoms if (! defined($select));
    $i = 1;
    while (! exists($atoms->{$i})) {
	$i++;
    }

    $counter = 0;
    for $i (sort numerically keys %{ $atoms }) {
	if (! exists($USED->{$i}) and exists($select->{$i})) {
	    $counter++;
	    $MOLS->{$counter} = ();
	    $MOLS->{$counter}{INDEX} = $counter;
	    $MOLS->{$counter}{MOLSIZE} = 0;
	    &getMolList($atoms, $bonds, $i, $MOLS->{$counter}, $USED);
	}
    }

    return $MOLS;
}

sub GetMols_old {
    my ($ATOMS, $BONDS) = @_;
    my ($atomC, @tmp, $rec, $i, $min, $molNum);

    for $atomC (keys %{ $ATOMS }) {
	delete $ATOMS->{$atomC}{MOLECULE};
	delete $ATOMS->{$atomC}{MOLECULEID};
    }
    
    @tmp = sort numerically keys %{ $ATOMS };
    for $atomC (@tmp) {
	$rec = ();
	$rec->{MEMBERS}{$atomC} = 1;
	if ($ATOMS->{$atomC}{NUMBONDS} == 0 || ! $BONDS->{$atomC}) { #ions
	    $rec->{INDEX} = $molNum;
	    $ATOMS->{$atomC}{MOLECULE} = \%{ $rec };
	    $molNum++;
	} else {
	    $min = $atomC;
	    for $i (@{ $BONDS->{$atomC} })  {
		$rec->{MEMBERS}{$i} = 1;
		$min = $i if ($i < $min);
	    }
	    if ($min < $atomC) { # found a member which has a lower index, so merge this rec with it, and update
		for $i (keys %{ $rec->{MEMBERS} }) {
		    $ATOMS->{$min}{MOLECULE}{MEMBERS}{$i} = 1;
		    next if ($min == $i or $i == $atomC);
		    if (exists($ATOMS->{$i}{MOLECULE}{MEMBERS})) {
			for (keys %{ $ATOMS->{$i}{MOLECULE}{MEMBERS} }) {
			    $ATOMS->{$min}{MOLECULE}{MEMBERS}{$_} = 1;
			}
		    }
		}
		$rec = \%{ $ATOMS->{$min}{MOLECULE} };
	    } else {
		$rec->{INDEX} = $molNum;
		$molNum++;
	    }

	    for $i (keys %{ $rec->{MEMBERS} })  {
		next if ($i == $min);
		$ATOMS->{$i}{MOLECULE} = \%{ $rec };
	    }
	}
    }

    $molNum = 1;
    for $atomC (@tmp) {
	if (! exists($ATOMS->{$atomC}{MOLECULE}{INDEX}) || $ATOMS->{$atomC}{NUMBONDS} == 0 || ! $BONDS->{$atomC}) {
	    $ATOMS->{$atomC}{MOLECULE}{INDEX} = $molNum;
	    $molNum++;
	}
	$ATOMS->{$atomC}{MOLECULEID} = $ATOMS->{$atomC}{MOLECULE}{INDEX};
	$ATOMS->{$atomC}{MOLECULE}{SIZE} = scalar(keys %{ $ATOMS->{$atomC}{MOLECULE}{MEMBERS} });
	$ATOMS->{$atomC}{MOLSIZE} = $ATOMS->{$atomC}{MOLECULE}{SIZE};
    }
}

sub SplitAtomsByMol {
    my ($atoms, $selList) = @_;
    my ($i, %molList, $counter, $j, @tmp);

    $selList = $atoms if (! defined($selList));
    @tmp = sort { ($a<=>$b) } keys %{ $selList };

    for $i (@tmp) {
	next if (exists($atoms->{$i}{IS_SPLIT}));
	for $j (keys %{ $atoms->{$i}{MOLECULE}{MEMBERS} }) {
	    next if (! exists($selList->{$j}));
	    $atoms->{$j}{IS_SPLIT} = 1;
	    $counter = $atoms->{$j}{MOLECULEID};
	    $molList{$counter}{$j} = $atoms->{$i}{MOLECULE}{MEMBERS}{$j};
	}
    }

    for $i (@tmp) {
        delete $atoms->{$i}{IS_SPLIT}
    }

    return \%molList;
}

sub GetAtmData {
    my ($allAtoms, $atomList) = @_;
    my (%ATOMS, $i);
                                                                                                                                      
    for $i (keys %{ $atomList }) {
        $ATOMS{$i} = \%{ $allAtoms->{$i} };
    }
    return \%ATOMS;
}

sub vecLen {
    my ($v) = $_[0];
    my ($vLen);

    for ("XCOORD", "YCOORD", "ZCOORD") {
	$vLen += $v->{$_}*$v->{$_};
    }
    return sqrt($vLen);
}

sub GetAbituraryRotationMatrix {
    my ($v, $s, $c) = @_;
    my ($rM,$vx,$vy,$vz,$scale,$tmp);
    $tmp = vecLen($v);
    $scale = 1/$tmp;
    ($vx,$vy,$vz) = ($v->{XCOORD}*$scale,$v->{YCOORD}*$scale,$v->{ZCOORD}*$scale);
    $rM = [
           [ 1+(1-$c)*($vx*$vx-1), (1-$c)*$vx*$vy-$vz*$s, (1-$c)*$vx*$vz+$vy*$s],
           [(1-$c)*$vx*$vy-$vz*$s,  1+(1-$c)*($vy*$vy-1), (1-$c)*$vy*$vz-$vx*$s],
           [(1-$c)*$vx*$vz-$vy*$s, (1-$c)*$vy*$vz+$vx*$s,  1+(1-$c)*($vz*$vz-1)]
          ];
    return $rM;
}

sub RotateAbout {
    my ($atoms, $rotM) = @_;
    my ($orig, $i, $j, $count);

    for $i (keys %{ $atoms }) {
        for $j ("XCOORD","YCOORD","ZCOORD") {
            $orig->{$i}{$j} = $atoms->{$i}{$j};
        }
    }

    for $i (keys %{ $atoms }) {
        $atoms->{$i}{XCOORD} =  $orig->{$i}{XCOORD}*$rotM->[0][0]+
                                $orig->{$i}{YCOORD}*$rotM->[0][1]+
                                $orig->{$i}{ZCOORD}*$rotM->[0][2];
        $atoms->{$i}{YCOORD} =  $orig->{$i}{XCOORD}*$rotM->[1][0]+
                                $orig->{$i}{YCOORD}*$rotM->[1][1]+
                                $orig->{$i}{ZCOORD}*$rotM->[1][2];
        $atoms->{$i}{ZCOORD} =  $orig->{$i}{XCOORD}*$rotM->[2][0]+
                                $orig->{$i}{YCOORD}*$rotM->[2][1]+
                                $orig->{$i}{ZCOORD}*$rotM->[2][2];
    }
}

sub TransAtoms {
    my ($atoms, $com, $d) = @_;
    my ($i, $j);

    $d = 1 if (! defined($d));
    for $i ("XCOORD", "YCOORD", "ZCOORD") {
	for $j (keys %{ $atoms }) {
	    $atoms->{$j}{$i} += $d*$com->{$i};
	}
    }
}

sub updateAtomIndex {
    my ($atoms, $bonds) = @_;
    my ($i, $j, $index, $newAtoms, $newBonds, @fields, @list);

    $newBonds = ();
    @list = keys %{ $atoms };
    @fields = grep {!/MOL/} keys %{ $atoms->{ $list[0] } };
    for $i (@list) {
        $index = $atoms->{$i}{INDEX};
	for $j (@fields) {
	    $newAtoms->{$index}{$j} = $atoms->{$i}{$j};
	}
        $newBonds->{$index} = ();
        next if (! exists($newAtoms->{$index}{BONDS}) || ! $newAtoms->{$index}{BONDS});
        for $j (0 .. $#{ $newAtoms->{$index}{BONDS} }) {
            $newBonds->{$index}[$j] = $newAtoms->{$index}{BONDS}[$j]{INDEX};
        }
        delete $newAtoms->{$index}{BONDS};
	delete $atoms->{$i};
	delete $bonds->{$i};
    }
    @list = keys %{ $atoms };
    for $i (@list) {
	delete $atoms->{$i};
	delete $bonds->{$i};
    }
    for $i (keys %{ $newAtoms }) {
	%{ $atoms->{$i} } = %{ $newAtoms->{$i} };
	$bonds->{$i} = ();
	@{ $bonds->{$i}  } = @{ $newBonds->{$i} } if($newBonds->{$i});
    }
}

sub GroupAtomsByField {
    my ($atoms, $bonds, $mode, $select, $mode2, $reverse_opt) = @_;
    my ($i, $j, $FIELD, $index, @sorted, $sort_val); 
    my ($resid, $resNum, $MOLS, $tmp, $hasSelect);

    $reverse_opt = 0 if (! defined($reverse_opt));
    $hasSelect = 1;
    if(! defined($select)) {
      $hasSelect = 0;
      $select = $atoms;
    }
    $mode2 = "INDEX" if (! defined($mode2));

    if ($mode =~ /MOL/ or $mode2 =~ /MOL/) {
	$MOLS = &GetMols($atoms, $bonds, $select);
    }

    for $i (keys %{ $atoms }) {
	$sort_val = $atoms->{$i}{$mode2};
	$sort_val = ${ $atoms->{$i}{$mode2} } if ($mode2 =~ /MOL/);
	if ($mode =~ /MOL/ and exists($select->{$i})) {
	    $FIELD->{ ${ $atoms->{$i}{$mode} } }{$i} = $sort_val;
	} elsif (exists($select->{$i})) {
	    $FIELD->{ $atoms->{$i}{$mode} }{$i} = $sort_val;
	}
	next if (! exists($bonds->{$i}) || ! $bonds->{$i});
	$index = 0;
	for $j (@{ $bonds->{$i} }) {
	    $atoms->{$i}{BONDS}[$index] = \%{ $atoms->{$j} };
	    $index++;
	}
    }

    if ($mode =~ /MOL|RESNUM|NUMBONDS/) {
	@sorted = sort {$a<=>$b} keys %{ $FIELD };
    } else {
	@sorted = sort {$a cmp $b } keys %{ $FIELD };
    }

    $index = $resid = $resNum = 0;
    if($hasSelect) {
	@{ $tmp } = sort {$a<=>$b} keys %{ $FIELD->{$sorted[0]} };
	$index = $atoms->{$tmp->[0]}{INDEX}-1;
	$resid = $resNum = $atoms->{$tmp->[0]}{RESNUM};
    }
    for $i (0 .. $#sorted) {
	@{ $tmp } = sort { $FIELD->{$sorted[$i]}{$a} cmp $FIELD->{$sorted[$i]}{$b} } keys %{ $FIELD->{$sorted[$i]} };
	@{ $tmp } = reverse sort @{ $tmp } if ($reverse_opt);
	for $j (@{ $tmp }) {
	    $index++;
	    $atoms->{$j}{INDEX} = $index;
	    if (! $resNum || $atoms->{$j}{RESNUM} != $resNum) {
		$resid++;
		$resNum = $atoms->{$j}{RESNUM};
	    }
	    $atoms->{$j}{RESNUM} = $resid;
	    if ($mode =~ /MOLECULEID/) {
		$atoms->{$j}{CHAIN} = "X";
		$atoms->{$j}{CHAIN} = chr(64+ ${ $atoms->{$j}{MOLECULEID} }) 
		    if($atoms->{$j}{MOLECULEID} < 10);
	    } elsif ($mode =~ /MOLSIZE/) {
		$atoms->{$j}{CHAIN} = chr(64+ $i+1);
	    }
	}
    }
    $FIELD = ();
    &updateAtomIndex(\%{ $atoms },\%{ $bonds });
}

sub MakeFieldSequential {
    my ($atoms, $field) = @_;
    my ($i, $atomsByField, @list, $isnumeric, $index, $j);

    @list = keys %{ $atoms };
    $isnumeric = 1;
    die "ERROR: field $field does not exists in MakeFieldSequential\n"
	if (! exists($atoms->{$list[0]}{$field}));
    for $i (@list ) {
	push @{ $atomsByField->{ $atoms->{$i}{$field} } }, $atoms->{$i};
	$isnumeric=0 if ($atoms->{$i}{$field} !~ /^\-?\d+\.?\d*$/);
    }
    return 0 if (scalar(keys %{ $atomsByField} ) == 1);
    @list = sort numerically keys %{ $atomsByField } if($isnumeric);
    @list = sort { ($a cmp $b); } keys %{ $atomsByField } if(!$isnumeric);
    $index = shift @list;
    $index++;
    for $i (@list) {
	for $j (@{ $atomsByField->{$i} }) {
	    $j->{$field} = $index;
	}
	$index++;
    }
}

sub BuildAtomSelectionString {
    my ($atomSel) = $_[0];
    my ($fields, $fl, $ud, $sel, $atm_no);

    $fields = (
	       {
		   "INDEX"       => 0,
		   "ATMNAME"     => 1,
		   "RESNAME"     => 1,
		   "CHAIN"       => 1,
		   "RESNUM"      => 0,
		   "XCOORD"      => 0,
		   "YCOORD"      => 0,
		   "ZCOORD"      => 0,
		   "FFTYPE"      => 1,
		   "NUMBONDS"    => 0,
		   "LONEPAIRS"   => 0,
		   "CHARGE"      => 0,
		   "MOLECULEID"  => 0,
		   "MOLSIZE"     => 0,
	       }
	      );
   $fl = "(" . join("|",keys %{ $fields }) . ")";
   $sel = $atomSel;
   while ($atomSel =~ /$fl/gi) {
	$ud = '$atoms->{$i}{' .uc  $1 . '}';
	$ud = '${ $atoms->{$i}{' .uc  $1 . '} }' if (uc($1) eq "MOLSIZE" or uc($1) eq "MOLECULEID");
	$sel =~ s/$1/$ud/;
   }
    while ($sel =~ /dist\((\w+)\)/) { 
        $atm_no = $1;
        $sel =~ s/dist\(\w+\)/GetBondLength\(\$atoms->\{\$i\}\,\$atoms->\{$atm_no\},\$box\)/;
    }
   return $sel;
}

sub SelectAtoms {
    my ($selectionStr, $atoms, $box) = @_;
    my ($atomList, $i, $err);

    #first try to see if expresion is valid
    $i = 1;
    #try
    $SIG{__WARN__} = sub {  };
    eval($selectionStr);
    if($@) {
	$err = $@;
	$err =~ s/^.* line \d+, near .../near \"/;
	die "ERROR: Invalid atom selection ${err}";
    }
    #now select atoms
    for $i (keys %{ $atoms }) {
	if (eval($selectionStr)) {
	    $atomList->{$i} = 1;
	}
    }
    # make sure at least 1 atom record selection
    die "ERROR: No atoms matched selection!\n" if (! $atomList);
    return $atomList;
}

sub ReimageAtoms {
    my ($atoms, $bonds, $mols, $box, $selection) = @_;
    my ($i, $j, $k, $l, $dist, $factor, $altered, $com, $curr);

    $altered = ();

    for $i (keys %{ $mols }) {
	for $j (keys %{ $mols->{$i}{MEMBERS} }) {
	    next if(!exists($selection->{$j}));
	    for $k (@{ $bonds->{$j} }) {
		next if($j > $k);
		for $l ("X", "Y", "Z") {
		    $dist = getdist($atoms->{$j}, $atoms->{$k}, $l);
		    $factor = 1;
		    $factor = -1 if ($atoms->{$k}{"${l}COORD"} < $atoms->{$j}{"${l}COORD"});
		    while($dist>4) {
			$dist -= $box->{$l};
			$atoms->{$k}{"${l}COORD"} -= $factor * $box->{$l};
		    }
		}
	    }
	}
	$curr = GetAtmData($atoms, $mols->{$i}{MEMBERS});
	$com = CoM($curr);
	for $l ("X", "Y", "Z") {
	    next if ($com->{"${l}COORD"} > 0 and $com->{"${l}COORD"} < $box->{$l});
	    $factor = 0;
	    if ($com->{"${l}COORD"} < 0) {
		while($com->{"${l}COORD"}<0) {
		    $factor += $box->{$l};
		    $com->{"${l}COORD"} += $box->{$l};
		}
	    } else {
		while($com->{"${l}COORD"}>$box->{$l}) {
		    $factor -= $box->{$l};
		    $com->{"${l}COORD"} -= $box->{$l};
		}
	    }
	    for $j (keys %{ $mols->{$i}{MEMBERS} }) {
		$atoms->{$j}{"${l}COORD"} += $factor;
	    }
	}
    }
}

sub getdist {
    my ($atom1, $atom2, $dim) = @_;

    return (sqrt(($atom1->{"${dim}COORD"}-$atom2->{"${dim}COORD"})**2));
}

1;

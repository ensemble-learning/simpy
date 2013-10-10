package Packages::FSearch;
use Exporter ();

our @ISA = qw(Exporter);
our @EXPORT = qw();
our @EXPORT_OK = qw(&GetResiduals &SetResiduals &StoreCoords &GetGridATOMS);

$VERSION = 1.0;
use strict;
use Packages::ManipAtoms qw(GetAtmData);
use Packages::General qw(CoM);

sub numerically { ($a<=>$b); }

sub GetResiduals {
    my ($arrayIndex, $distArray, $residuals, $validRes, $tot, $boxLen) = @_;
    my ($i, $baseVal, $offset, $j, $atomC, $tmp, $count);
                                                                                                                             
    $baseVal = $distArray->[$arrayIndex][0];
    $count = 0;
                                                                                                                             
    for $i (1 .. $tot) {
        $tmp->[$i] = $validRes->[$i];
        $validRes->[$i] = 0;
    }
                                                                                                                             
    for ($j = ($arrayIndex + 1); $j < $tot; $j++) {
        $atomC = $distArray->[$j][1];
        next if (! $tmp->[$atomC]);
        $offset = $distArray->[$j][0] - $baseVal;
        $offset -= $boxLen while ($offset > $boxLen);
        $residuals->[$atomC] -= ($offset * $offset);
        last if ($residuals->[$atomC] < 0);
        $validRes->[$atomC] = $atomC;
        $count++;
    }
                                                                                                                             
    for ($j = ($arrayIndex - 1); $j > -1; $j--) {
        $atomC = $distArray->[$j][1];
        next if (! $tmp->[$atomC]);
        $offset = ($baseVal - $distArray->[$j][0]);
        $offset -= $boxLen while ($offset > $boxLen);
        $residuals->[$atomC] -= ($offset * $offset);
        last if ($residuals->[$atomC] < 0);
        $validRes->[$atomC] = $atomC;
        $count++;
    }
                                                                                                                             
    return $count;
}

sub SetResiduals {
    my ($solvAtoms, $distMax, $skip, $tot, $residuals, $vRes) = @_;
    my ($i);
                                                                                                                             
    for $i (1 .. $tot) {
        $residuals->[$i] = $distMax;
        $vRes->[$i] = $i;
        $vRes->[$i] = 0 if (! $solvAtoms->[$i]);
    }
    $vRes->[$skip] = 0;
}
                                                                                                                             
sub StoreCoords {
    my ($atoms, $box) = @_;
    my (%cSORTED, $DATA, $i, $d, %atomMap, $index, @arryIndex, $tot);
                                                                                                                             
    @arryIndex =keys %{ $atoms };
    $tot = $#arryIndex;
    for $i (@arryIndex) {
        $d = $atoms->{$i}{XCOORD};
        $d += 0.00001 while (exists($DATA->{X}{$d}));
        $DATA->{X}{$d} = $i;
                                                                                                                             
        $d = $atoms->{$i}{YCOORD};
        $d += 0.00001 while (exists($DATA->{Y}{$d}));
        $DATA->{Y}{$d} = $i;
                                                                                                                             
        $d = $atoms->{$i}{ZCOORD};
        $d += 0.00001 while (exists($DATA->{Z}{$d}));
        $DATA->{Z}{$d} = $i;
    }
                                                                                                                             
    for $i ("X", "Y", "Z") {
        $index = 0;
        $#{ $cSORTED{$i} } = $tot;
        for $d (sort numerically keys %{ $DATA->{$i} }) {
            $cSORTED{$i}[$index][0] = $d;
            $cSORTED{$i}[$index][1] = $DATA->{$i}{$d};
            $atoms->{ $DATA->{$i}{$d} }{SORT}{$i} = $index;
            $index++;
        }
    }
                                                                                                                             
    return \%cSORTED;
}

sub GetGridATOMS {
    my ($allAtoms, $soluMols, $solvMols) = @_;
    my ($i, $count, %atomCOMS, @solu, @solv, $MOL, %SOLVSOLU, @tmp);
                                                                                                                             
    $count = 0;
    $solu[0] = $solv[0] = 0;
    for $i (keys %{ $soluMols }) {
        $count++;
        $MOL = GetAtmData($allAtoms, $soluMols->{$i});
        $atomCOMS{$count} = CoM($MOL);
        $atomCOMS{$count}{IS_SOLUTE} = 1;
        $atomCOMS{$count}{MOLECULEID} = $i;
        $atomCOMS{$count}{RESNAME} = "SOL";
        $solu[$count] = $count;
        @tmp  = keys %{ $soluMols->{$i} };
        if ($allAtoms->{ $tmp[0] }{IS_SOLVENT}) {
            $SOLVSOLU{ $allAtoms->{ $tmp[0] }{SOLVENTID} } = 1;
            $atomCOMS{$count}{RESNAME} = "WAT";
            $solv[$count] = $count;
        } else {
            $solv[$count] = 0;
        }
    }
                                                                                                                             
    for $i (keys %{ $solvMols }) {
        next if (exists($SOLVSOLU{$i}));
        $count++;
        $MOL = GetAtmData($allAtoms, $solvMols->{$i});
        $atomCOMS{$count} = CoM($MOL);
        $atomCOMS{$count}{IS_SOLUTE} = 1;
        $atomCOMS{$count}{MOLECULEID} = $i;
        $atomCOMS{$count}{RESNAME} = "WAT";
        $solv[$count] = $count;
        $solu[$count] = 0;
    }
                                                                                                                             
    return (\%atomCOMS, \@solu, \@solv);
}

1;

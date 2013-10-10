package Packages::REPLICATE;

require Exporter;
use strict;
use Cwd;
use Packages::General qw(PrintProgress CombineMols GetTime);

use constant PI => atan2(1,1) * 4;

our (@EXPORT_OK, @ISA, @EXPORT, $VERSION);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(ReplicateCell GetBoxVol);
$VERSION = "1.00";

sub transMol {
    my ($unit, $box, $dim, $disp) = @_;
    my ($i,%ATOMS, $j, $offset);
    
    for $i (keys %{ $unit } ) {
	%{ $ATOMS{$i} } = %{ $unit->{$i} };
    }
    for $i (keys %{ $unit }) {
	for $j ("X","Y","Z") {
	    $offset = $box->{DISP_V}{$j} * $disp;
	    $ATOMS{$i}{"${j}COORD"} += $offset;
	}
    }
    return \%ATOMS;
}

sub invertMol {
    my ($atoms, $dim) = @_;
    my ($i);

    for $i (keys %{ $atoms }) {
	$atoms->{$i}{"${dim}COORD"} *= -1;
    }
}

sub GetBoxVol {
	my ($box) = $_[0];

	&getBoxDisplacementTensor($box);

}

sub getBoxDisplacementTensor {
    my ($box) = $_[0];
    my ($lx,$ly,$lz,$a,$b,$c,$cos_alpha,$alpha,$beta,$gamma);
    my ($ax,$bx,$cx,$by,$cy,$cz);

    $lx = $box->{X}{len};
    $ly = $box->{Y}{len};
    $lz = $box->{Z}{len};
    $alpha = $box->{X}{angle}*PI/180;
    $beta  = $box->{Y}{angle}*PI/180;
    $gamma = $box->{Z}{angle}*PI/180;

    # 3 x 3 displacement tensor
    #(a b c) =  (ax bx cx)
    #		(0  by cy)
    #		(0  0  cz)
    #ax = lx
    #bx = ly cos(gamma)
    #cx = lz cos(beta)
    #cy = ly*lz*cos(alpha) - bx*cx/by
    #cz = sqrt(lz*lz - cx*cx - cy*cy)
    $ax = $lx; $bx = $ly*cos($gamma); $cx = $lz*cos($beta);
    $by = $ly*sin($gamma); $cy = ($ly*$lz*cos($alpha)-$bx*$cx)/$by;
    $cz = sqrt($lz*$lz - $cx*$cx - $cy*$cy);
    $box->{X}{DISP_V}{X} = $ax; $box->{X}{DISP_V}{Y} = 0;   $box->{X}{DISP_V}{Z} = 0;
    $box->{Y}{DISP_V}{X} = $bx; $box->{Y}{DISP_V}{Y} = $by; $box->{Y}{DISP_V}{Z} = 0;
    $box->{Z}{DISP_V}{X} = $cx; $box->{Z}{DISP_V}{Y} = $cy; $box->{Z}{DISP_V}{Z} = $cz;
}

sub ReplicateCell {
    my ($atoms, $bonds, $box, $dims, $centerMol, $createNewbonds, $str, $updateResNum, $invertDim) = @_;
    my ($i, $j, $cellAtoms, $cellBonds, $unitAtoms, $pbcbonds); 
    my ($strLen, $offset, $tot, $count, $total, $start);
    
    &getBoxDisplacementTensor($box);
    print "${str}calculating time remaining\r";
    $start = time();
    if ($centerMol) {
	for $i ("X", "Y", "Z") {
	    $j = -1 * int(($dims->{$i} -1)/2);
	    for ($j .. -1) {
		$atoms = transMol($atoms, $box->{$i}, $i, $_);
	    }
	}
    }

    $total = 1;
    for $i ("X", "Y", "Z") {
	$total *= $dims->{$i};
    }

    for $i ("X", "Y", "Z") {
	$unitAtoms = ();
	$cellBonds = ();
	$pbcbonds = ();
	$pbcbonds = getPBCbonds($atoms, $bonds, $box) if ($createNewbonds);
	$tot = 0;
	$offset = 0;
	for $j (keys %{ $atoms }) {
	    %{ $unitAtoms->{$j} } = %{ $atoms->{$j} };
	    $cellBonds->{$j} = ();
	    @{ $cellBonds->{$j} } = @{ $bonds->{$j} } if($bonds->{$j});
	    $tot++;
	}
	for $j (1 .. ($dims->{$i} - 1)) {
	    $count++;
	    $strLen = PrintProgress($count, $total, $start, $str);
	    $offset += $tot;
	    $cellAtoms = transMol($unitAtoms, $box->{$i}, $i, $j);
	    &invertMol($cellAtoms, $invertDim) if (defined($invertDim) and ($count % 2 == 1));
	    ($atoms, $bonds) = CombineMols($atoms, $cellAtoms, $bonds, $cellBonds, $updateResNum);
	    updatePBCbonds($atoms, $bonds, $pbcbonds->{$i}, $offset, $tot, $i) if(exists($pbcbonds->{$i}));
	}
	$box->{$i}{hi} = ($box->{$i}{len} * $dims->{$i});
	$box->{$i}{lo} = 0;
	$box->{$i}{len} = ($box->{$i}{len} * $dims->{$i});
    }
    $tot = GetTime(time() - $start);
    printf "$str%-${strLen}s\n", "${tot}s elapsed...Done";
    return ($atoms, $bonds, $box);
}

sub updatePBCbonds {
    my ($atoms, $bonds, $bondlist, $atomOffset, $tot, $dim) = @_;
    my ($i, $j, $atom1, $atom2, $atom3, $atom4, @list, $pos);

    for $i (keys %{ $bondlist }) {
	@list = keys %{ $bondlist->{$i} };
	for $j (0 .. $#list) {
	    $pos = $bondlist->{$i}{ $list[$j] };
	    $atom1 = $i;
	    $atom2 = $list[$j];
	    $atom3 = $atom1 + $atomOffset;
	    $atom4 = $atom2 + $tot;
	    #for atom1 , delete bond to atom2 and form bond to atom4
	    $bonds->{$atom1}[$pos->[0]] = $atom4;
            #for atom2, delete bond to atom1 and form bond to atom3
	    $bonds->{$atom2}[$pos->[1]] = $atom3;
	    delete $atoms->{$atom2}{"DISP${dim}"};
	    #$atoms->{$atom2}{"DISP${dim}"}[$pos->[1]] = 0;
            #for atom3, delete bond to atom4 and form bond to atom2
            $bonds->{$atom3}[$pos->[0]] = $atom2;
	    delete $atoms->{$atom3}{"DISP${dim}"};
            #$atoms->{$atom3}{"DISP${dim}"}[$pos->[0]] = 0;
            #updateBond($atoms->{$atom3}, $bonds->{$atom3}, $atom4, $atom2);
            #delete $atoms->{$atom3}{"DISP${dim}"};
            #for atom4, delete bond to atom3 and form bond to atom1
            $bonds->{$atom4}[$pos->[1]] = $atom1;
            #updateBond($atoms->{$atom4}, $bonds->{$atom4}, $atom3, $atom1);
            #now update bondlist
            delete $bondlist->{$i}{$atom2};
            $bondlist->{$i}{$atom4}[0] = $pos->[0];
	    $bondlist->{$i}{$atom4}[1] = $pos->[1];
	}
    }
    print "";
}

sub getPBCbonds {
    my ($atoms, $bonds, $box) = @_;
    my ($i, $j, $k, $l, $dist, $atom1, $atom2, $bondlist, $sign,$flag2);

    #search for multiple bond entries to same atom and image flags
    for $i (keys %{ $bonds }) {
	$atom1 = $i;
        for $j (0 .. $#{ $bonds->{$i} }) {
            for $k ("X", "Y", "Z") {
		$sign = 0;
		$sign = $atoms->{$i}{"DISP${k}"}[$j] if (exists($atoms->{$i}{"DISP${k}"}));
		next if($sign>-1);
		#find symmetric entry
		undef($flag2);
		$atom2 =  $bonds->{$i}[$j];
		for $l (0 .. $#{ $bonds->{$atom2} }) {
		    if ($bonds->{$atom2}[$l] == $atom1) {
			if (exists($atoms->{$atom2}{"DISP${k}"}) and $atoms->{$atom2}{"DISP${k}"}[$l] == -1*$sign) {
			    $flag2 = $l;
			    last;
			}
		    }
		}
		next if (! defined($flag2));
		$bondlist->{$k}{$i}{$bonds->{$i}[$j]} = ([$j, $flag2]);
            }
        }
    }

    return $bondlist if(defined($bondlist));
    #search for bonds by distance
    for $i (keys %{ $bonds }) {
        $atom1 = \%{ $atoms->{$i} };
	for $j (0 .. $#{ $bonds->{$i} }) {
            $atom2 = \%{ $atoms->{$bonds->{$i}[$j]} };
	    for $k ("X","Y","Z") {
		$dist = $atom1->{"${k}COORD"} - $atom2->{"${k}COORD"};
		next if ($dist == 0);
		$sign = $dist/abs($dist);
		if ($sign < 0) {
		    $dist *= -1;
		}
		if ($dist > 5 and $dist > $box->{$k}{len}/3) { #search for long bonds
		    $bondlist->{$k}{$bonds->{$i}[$j]}{$i} = ([$j,$sign]);
		}
	    }
	}
    }

    return $bondlist;
}
1;

package Packages::Superimpose;
use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = ();
our @EXPORT_OK = qw(SuperimposeAtoms);

sub numerically { ($a<=>$b); }

sub calcKearsleyMatrix {
    my ($ref, $mov) = @_;
    my (@km, $i, $j, @refMtchAtms, @movMtchAtms);
    my ($Xm, $Xp, $Ym, $Yp, $Zm, $Zp);
    my ($Xm2, $Xp2, $Ym2, $Yp2, $Zm2, $Zp2);
    
    for $i (0 .. 3) {
        for $j (0 .. 3) {
            $km[$i][$j] = 0;
        }
    }
    
    $j = 0;
    for $i (sort numerically keys %{ $ref }) {
        $refMtchAtms[$j][0] = $ref->{$i}{"XCOORD"};
        $refMtchAtms[$j][1] = $ref->{$i}{"YCOORD"};
        $refMtchAtms[$j][2] = $ref->{$i}{"ZCOORD"};
        $j++;
    }
    
    $j = 0;
    for $i (sort numerically keys %{ $mov }) {
        $movMtchAtms[$j][0] = $mov->{$i}{"XCOORD"};
        $movMtchAtms[$j][1] = $mov->{$i}{"YCOORD"};
        $movMtchAtms[$j][2] = $mov->{$i}{"ZCOORD"};
        $j++;
    }
    
    die "ERROR: Unequal number of atoms in ref and match\n" if ($#refMtchAtms != $#movMtchAtms);
    
    $j--;
    for $i (0 .. $j) {
	$Xm = -1*($movMtchAtms[$i][0] - $refMtchAtms[$i][0]);
	$Ym = -1*($movMtchAtms[$i][1] - $refMtchAtms[$i][1]);
	$Zm = -1*($movMtchAtms[$i][2] - $refMtchAtms[$i][2]);

	$Xp = ($refMtchAtms[$i][0] + $movMtchAtms[$i][0]);
	$Yp = ($refMtchAtms[$i][1] + $movMtchAtms[$i][1]);
	$Zp = ($refMtchAtms[$i][2] + $movMtchAtms[$i][2]);

	$Xm2 = $Xm*$Xm;
	$Ym2 = $Ym*$Ym;
	$Zm2 = $Zm*$Zm;

	$Xp2 = $Xp*$Xp;
	$Yp2 = $Yp*$Yp;
	$Zp2 = $Zp*$Zp;

	$km[0][0] = $km[0][0] +  ($Xm2 + $Ym2 + $Zm2);
	$km[0][1] = $km[0][1] +  ($Yp*$Zm - $Ym*$Zp);
	$km[0][2] = $km[0][2] +  ($Xm*$Zp - $Xp*$Zm);
	$km[0][3] = $km[0][3] +  ($Xp*$Ym - $Xm*$Yp);
	$km[1][0] = $km[1][0] +  ($Yp*$Zm - $Ym*$Zp);
	$km[1][1] = $km[1][1] +  ($Yp2 + $Zp2 + $Xm2);
	$km[1][2] = $km[1][2] +  ($Xm*$Ym - $Xp*$Yp);
	$km[1][3] = $km[1][3] +  ($Xm*$Zm - $Xp*$Zp);
	$km[2][0] = $km[2][0] +  ($Xm*$Zp - $Xp*$Zm);
	$km[2][1] = $km[2][1] +  ($Xm*$Ym - $Xp*$Yp);
	$km[2][2] = $km[2][2] +  ($Xp2 + $Zp2 + $Ym2);
	$km[2][3] = $km[2][3] +  ($Ym*$Zm - $Yp*$Zp);
	$km[3][0] = $km[3][0] +  ($Xp*$Ym - $Xm*$Yp);
	$km[3][1] = $km[3][1] +  ($Xm*$Zm - $Xp*$Zp);
	$km[3][2] = $km[3][2] +  ($Ym*$Zm - $Yp*$Zp);
	$km[3][3] = $km[3][3] +  ($Xp2 + $Yp2 + $Zm2);        
    }
    
    return \@km;
}

sub findDiagonalizingEuler {
    my ($evec) = $_[0];
    my ($e0, $e1, $e2, $e3);
    my ($ee0, $ee1, $ee2, $ee3);
    my (@Rot_Mat);
    
    $e0 = $evec->[0][3];
    $e1 = $evec->[1][3];
    $e2 = $evec->[2][3];
    $e3 = $evec->[3][3];

    $ee0 = $e0*$e0;
    $ee1 = $e1*$e1;
    $ee2 = $e2*$e2;
    $ee3 = $e3*$e3;


    $Rot_Mat[0][0] = $ee0 + $ee1 - $ee2 - $ee3;
    $Rot_Mat[0][1] = 2*($e1*$e2 + $e0*$e3);
    $Rot_Mat[0][2] = 2*($e1*$e3 - $e0*$e2);
    $Rot_Mat[1][0] = 2*($e1*$e2 - $e0*$e3);
    $Rot_Mat[1][1] = $ee0 - $ee1 + $ee2 - $ee3;
    $Rot_Mat[1][2] = 2*($e2*$e3 + $e0*$e1);
    $Rot_Mat[2][0] = 2*($e1*$e3 + $e0*$e2);
    $Rot_Mat[2][1] = 2*($e2*$e3 - $e0*$e1);
    $Rot_Mat[2][2] = $ee0 - $ee1 - $ee2 + $ee3;

    return \@Rot_Mat;   
}


sub SuperimposeAtoms {
    my ($refMol, $movMol, $masses) = @_;
    my ($ref_com, $mov_com, $kearsleyMatrix, $kearsleyEigenVec, $RotMat);
    my ($start, $end);
    
    $ref_com = CoM($refMol, $masses);
    $mov_com = CoM($movMol, $masses);
    #translate both molecules to origin
    &translateAllAtoms($refMol, $ref_com, -1);
    &translateAllAtoms($movMol, $mov_com, -1);
    
    $kearsleyMatrix = calcKearsleyMatrix($refMol, $movMol);
    #printMatrix($kearsleyMatrix);
    $kearsleyEigenVec = diagMatrix($kearsleyMatrix, 4);
    #printMatrix($kearsleyEigenVec);
    $RotMat = findDiagonalizingEuler($kearsleyEigenVec);
    #printMatrix($RotMat);
    &rotateMolecule($movMol, $RotMat);
    #translate systems back to com
    &translateAllAtoms($refMol, $ref_com, 1);
    &translateAllAtoms($movMol, $ref_com, 1);
}

sub printMatrix {
    my ($matrix) = $_[0];
    my ($i, $j, $outString);
    
    $outString = "\n{";
    for $i (@{ $matrix }) {
        $outString .= "{";
        for $j (@{ $i }) {
            $outString .= "$j,";
        }
        chop $outString;
        $outString .= "},";
    }
    chop $outString;
    $outString .= "}\n";
    
    print $outString;
}

sub rotateMolecule {
    my ($atoms, $rotMatrix) = @_;
    my (@angles, $PI, $i, $j, $atom, $xC, $yC, $zC, $dim);

    for $i (keys %{ $atoms }) {
        $atom = $atoms->{$i};
        $xC = $atom->{"XCOORD"};
        $yC = $atom->{"YCOORD"};
        $zC = $atom->{"ZCOORD"};
        $j = 0;
        for $dim ("XCOORD", "YCOORD", "ZCOORD") {
            $atom->{$dim} = $xC * $rotMatrix->[$j][0] + $yC * $rotMatrix->[$j][1] + $zC * $rotMatrix->[$j][2];
            $j++;
        }
    }
            
}

sub diagMatrix {
    my ($matrixA, $size) = @_;
    my ($new_matrix, $V, $E, $eigenVec);
    
    $new_matrix = Math::MatrixReal->new_from_rows([\@{ $matrixA->[0]}, \@{ $matrixA->[1]}, \@{ $matrixA->[2]}, \@{ $matrixA->[3]}]);
    #print "Input Matrix:\n";
    #print $new_matrix->as_yacas( ( format => "%.8f", align => "l",name => "A" ) );
    
    ($E, $V) = $new_matrix->sym_diagonalize();
    #print "\nEigenvalues: \n";
    #print $E->as_yacas( ( format => "%.5G", align => "l",name => "A" ) );    
    #print "\nEigenvectors: \n";
    #print $V->as_yacas( ( format => "%.5G", align => "l",name => "A" ) );
    $eigenVec = eigenSort($E, $V, $size);
    #print "\nBest Eigenvector: \n";
    #printMatrix($eigenVec);
    return $eigenVec;
}

sub eigenSort {
    my ($eVal, $eVec, $size) = @_;
    my ($i, $rec, %eigenVecs, $j, @tmp, @ret);
    
    for $i (1 .. $size) {
        $rec = ();
        for $j (1 .. $size) {
            push @{ $rec }, $eVec->element($i,$j);
        }
        $eigenVecs{ $eVal->element($i,1) } = $rec;
    }
   
    @tmp = sort numerically keys %eigenVecs;
    
    for $i (@tmp) {
        push @ret, $eigenVecs{$i};
    }
    
    return \@ret;
}

sub CoM {
    my ($atoms, $masses) = @_;
    my ($i, %COM, $dim, $totalMass, $mass);
    
    $totalMass = 0;
    for $i (keys %{ $atoms }) {
	die "ERROR: Atom $i does not have a type!\n" if (! keys %{ $atoms->{$i} } or ! exists($atoms->{$i}{"FFTYPE"}));
	if (exists($atoms->{$i}{MASS})) {
	    $mass = $atoms->{$i}{MASS};
	} elsif (defined($masses) and exists($masses->{ $atoms->{$i}{"FFTYPE"} })) {
            $mass = $masses->{ $atoms->{$i}{"FFTYPE"} };
	} else {
	    $mass = 1;
	}
        for $dim ("XCOORD", "YCOORD", "ZCOORD") {
            $COM{$dim} += $atoms->{$i}{$dim} * $mass;
        }
        $totalMass += $mass;
    }
    
    for $dim ("XCOORD", "YCOORD", "ZCOORD") {
        $COM{$dim} /= $totalMass;
    }
   
   return \%COM;
}

sub translateAllAtoms {
    my ($atoms, $transVec, $direction) = @_;
    my ($i, $dim);
    
    for $i (keys %{ $atoms }) {
        for $dim (keys %{ $transVec }) {
            $atoms->{$i}{$dim} += $direction * $transVec->{$dim};
        }
    }
}

1;

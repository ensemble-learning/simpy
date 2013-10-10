package Math::Polynomial::Solve;

require 5.005;

use Math::Complex;
use Carp;
use Exporter;
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION $epsilon $use_hessenberg);
use vars qw($debug_cubic  $debug_quartic  $debug_poly);
use strict;

@ISA = qw(Exporter);

#
# Export only on request.
#
@EXPORT_OK = qw(
	poly_roots
	linear_roots
	quadratic_roots
	cubic_roots
	quartic_roots
	get_hessenberg
	set_hessenberg
	GetLeastSquaresPlane
);

$VERSION = '2.10';
$use_hessenberg = 0;	# Set to 1 to force its use regardless of polynomial degree.
$debug_cubic = 0;
$debug_quartic = 0;
$debug_poly = 0;

#
# Set up the epsilon variable, the value that in the floating-point math
# of the computer is the smallest value a variable can have before it is
# indistinguishable from zero when adding it to one.
#

BEGIN
{
	my $epsilon2 = 0.25;
	$epsilon = 0.5;

	while (1.0 + $epsilon2 > 1.0)
	{
		$epsilon = $epsilon2;
		$epsilon2 /= 2.0;
	}

	print "\$Math::Polynomial::Solve::epsilon = ", $epsilon, "\n" if ($debug_poly);
}

#
# Get/Set the flag that tells the module to use the QR Hessenberg
# method regardless of the degree of the polynomial.
#
sub get_hessenberg
{
	return $use_hessenberg;
}

sub set_hessenberg
{
	$use_hessenberg = ($_[0])? 1: 0;
}

#
# @x = linear_roots($a, $b)
#
#
sub linear_roots(@)
{
	my($a, $b) = @_;

	if (abs($a) < $epsilon)
	{
		carp "The coefficient of the highest power must not be zero!\n";
		return (undef);
	}

	return wantarray? (-$b/$a): -$b/$a;
}

#
# @x = quadratic_roots($a, $b, $c)
#
#
sub quadratic_roots(@)
{
	my($a, $b, $c) = @_;

	if (abs($a) < $epsilon)
	{
		carp "The coefficient of the highest power must not be zero!\n";
		return (undef, undef);
	}

	return (0, -$b/$a) if (abs($c) < $epsilon);

	my $dis_sqrt = sqrt($b*$b - $a * 4 * $c);

	$dis_sqrt = -$dis_sqrt if ($b < $epsilon);	# Avoid catastrophic cancellation.

	my $xt = ($b + $dis_sqrt)/-2;

	return ($xt/$a, $c/$xt);
}

#
# @x = cubic_roots($a, $b, $c, $d)
#
#
sub cubic_roots(@)
{
	my($a, $b, $c, $d) = @_;
	my @x;

	if (abs($a) < $epsilon)
	{
		carp "The coefficient of the highest power must not be zero!\n";
		return @x;
	}

	return (0, quadratic_roots($a, $b, $c)) if (abs($d) < $epsilon);

	my $xN = -$b/3/$a;
	my $yN = $d + $xN * ($c + $xN * ($b + $a * $xN));

	my $two_a = 2 * $a;
	my $delta_sq = ($b * $b - 3 * $a * $c)/(9 * $a * $a);
	my $h_sq = 4/9 * ($b * $b - 3 * $a * $c) * $delta_sq**2;
	my $dis = $yN * $yN - $h_sq;
	my $twothirds_pi = (2 * pi)/3;

	#
	#              Debug
	# print "two_a      = $two_a\n" if ($debug_cubic);
	# print "delta_sq   = $delta_sq\n" if ($debug_cubic);
	# print "h_sq       = $h_sq\n" if ($debug_cubic);
	# print "dis        = $dis\n" if ($debug_cubic);
	#
	if ($dis > $epsilon)
	{
		print "Cubic branch 1: (\$dis >  0) " if ($debug_cubic);

		#
		# One real root, two complex roots.
		#
		my $dis_sqrt = sqrt($dis);
		my $r_p = $yN - $dis_sqrt;
		my $r_q = $yN + $dis_sqrt;
		my $p = cbrt( abs($r_p)/$two_a );
		my $q = cbrt( abs($r_q)/$two_a );

		$p = -$p if ($r_p > 0);
		$q = -$q if ($r_q > 0);

		$x[0] = $xN + $p + $q;
		$x[1] = $xN + $p * exp($twothirds_pi * i)
			    + $q * exp(-$twothirds_pi * i);
		$x[2] = ~$x[1];
	}
	elsif ($dis < -$epsilon)
	{
		print "Cubic branch 2: (\$dis <  0) " if ($debug_cubic);

		#
		# Three distinct real roots.
		#
		my $theta = acos(-$yN/sqrt($h_sq))/3;
		my $delta = sqrt($b * $b - 3 * $a * $c)/(3 * $a);
		my $two_d = 2 * $delta;

		@x = ($xN + $two_d * cos($theta),
			$xN + $two_d * cos($twothirds_pi - $theta),
			$xN + $two_d * cos($twothirds_pi + $theta));
	}
	else
	{
		print "Cubic branch 3: (\$dis == 0) " if ($debug_cubic);

		#
		# abs($dis) <= $epsilon (effectively zero).
		#
		# Three real roots (two or three equal).
		#
		my $delta = cbrt($yN/$two_a);

		@x = ($xN + $delta, $xN + $delta, $xN - 2 * $delta);
	}

	return @x;
}

#
# @x = quartic_roots($a, $b, $c, $d, $e)
#
#
sub quartic_roots(@)
{
	my($a,$b,$c,$d,$e) = @_;
	my @x = ();

	if (abs($a) < $epsilon)
	{
		carp "Coefficient of highest power must not be zero!\n";
		return @x;
	}

	return (0, cubic_roots($a, $b, $c, $d)) if (abs($e) < $epsilon);

	#
	# First step:  Divide by the leading coefficient.
	#
	$b /= $a;
	$c /= $a;
	$d /= $a;
	$e /= $a;

	#
	# Second step: simplify the equation to the
	# "resolvant cubic"  y**4 + fy**2 + gy + h.
	#
	# (This is done by setting x = y - b/4).
	#

	my $b4 = $b/4;

	#
	# The f, g, and h values are:
	#
	my $f = $c -
		6 * $b4 * $b4;
	my $g = $d +
		2 * $b4 * (-$c + 4 * $b4 * $b4);
	my $h = $e +
		$b4 * (-$d + $b4 * ($c - 3 * $b4 * $b4));


	if (abs($h) < $epsilon)
	{
		print "Quartic branch 1: (\$h == 0)\n" if ($debug_quartic);

		#
		# Special case: h == 0.  We have a cubic times y.
		#
		@x = (0, cubic_roots(1, 0, $f, $g));
	}
	elsif (abs($g) < $epsilon)
	{
		print "Quartic branch 2: (\$g == 0)\n" if ($debug_quartic);

		#
		# Another special case: g == 0.  We have a quadratic
		# with y-squared.
		#
		my($p, $q) = quadratic_roots(1, $f, $h);
		$p = sqrt($p);
		$q = sqrt($q);
		@x = ($p, -$p, $q, -$q);
	}
	else
	{
		print "Quartic branch 3: Ferrari's method.\n" if ($debug_quartic);

		#
		# Special cases don't apply, so continue on with Ferrari's
		# method.  This involves setting up the resolvant cubic
		# as the product of two quadratics.
		#
		# After setting up conditions that guarantee that the
		# coefficients come out right (including the zero value
		# for the third-power term), we wind up with a 6th
		# degree polynomial with, fortunately, only even-powered
		# terms.  In other words, a cubic with z = y**2.
		#
		# Take a root of that equation, and get the
		# quadratics from it.
		#
		my($z);
		($z, undef, undef) = cubic_roots(1, 2*$f, $f*$f - 4*$h, -$g*$g);

		print STDERR "\n\tz = $z\n" if ($debug_quartic);

		my $alpha = sqrt($z);
		my $rho = $g/$alpha;
		my $beta = ($f + $z - $rho)/2;
		my $gamma = ($f + $z + $rho)/2;

		@x = quadratic_roots(1, $alpha, $beta);
		push @x, quadratic_roots(1, -$alpha, $gamma);
	}

	return ($x[0] - $b4, $x[1] - $b4, $x[2] - $b4, $x[3] - $b4);
}


# Perl code to find roots of a polynomial translated by Nick Ing-Simmons
# <Nick@Ing-Simmons.net> from FORTRAN code by Hiroshi Murakami.
# From the netlib archive: http://netlib.bell-labs.com/netlib/search.html
# In particular http://netlib.bell-labs.com/netlib/opt/companion.tgz

#       BASE is the base of the floating point representation on the machine.
#       It is 16 for base 16 float : for example, IBM system 360/370.
#       It is 2  for base  2 float : for example, IEEE float.

my $MAX_ITERATIONS = 60;
sub BASE ()    { 2 }
sub BASESQR () { BASE * BASE }

#
# $matrix_ref = build_companion(@coefficients);
#
# Build the Companion Matrix of the N degree polynomial.  Return a
# reference to the N by N matrix.
#
sub build_companion(@)
{
	my @coefficients = @_;
	my $n = $#coefficients;
	my @h;			# The matrix.

	#
	# First step:  Divide by the leading coefficient.
	#
	my $cn = shift @coefficients;

	foreach my $c (@coefficients)
	{
		$c /= $cn;
	}

	#
	# Why would we be calling this for a linear equation?
	# Who knows, but if we are, then we can skip all the
	# complicated looping.
	#
	if ($n == 1)
	{
		$h[1][1] = -$coefficients[0];
		return \@h;
	}

	#
	# Next: set up the diagonal matrix.
	#
	for my $i (1 .. $n)
	{
		for my $j (1 .. $n)
		{
			$h[$i][$j] = 0.0;
		}
	}

	for my $i (2 .. $n)
	{
		$h[$i][$i - 1] = 1.0;
	}

	#
	# And put in the coefficients.
	#
	for my $i (1 .. $n)
	{
		$h[$i][$n] = - (pop @coefficients);
	}

	#
	# Now we balance the matrix.
	#
	# Balancing the unsymmetric matrix A.
	# Perl code translated by Nick Ing-Simmons <Nick@Ing-Simmons.net>
	# from FORTRAN code by Hiroshi Murakami.
	#
	#  The fortran code is based on the Algol code "balance" from paper:
	#  "Balancing a Matrix for Calculation of Eigenvalues and Eigenvectors"
	#  by B. N. Parlett and C. Reinsch, Numer. Math. 13, 293-304(1969).
	#
	#  Note: The only non-zero elements of the companion matrix are touched.
	#
	my $noconv = 1;
	while ($noconv)
	{
		$noconv = 0;
		for my $i (1 .. $n)
		{
			#
			# Touch only non-zero elements of companion.
			#
			my $c;
			if ($i != $n)
			{
				$c = abs($h[$i + 1][$i]);
			}
			else
			{
				$c = 0.0;
				for my $j (1 .. $n - 1)
				{
					$c += abs($h[$j][$n]);
				}
			}

			my $r;
			if ($i == 1)
			{
				$r = abs($h[1][$n]);
			}
			elsif ($i != $n)
			{
				$r = abs($h[$i][$i - 1]) + abs($h[$i][$n]);
			}
			else
			{
				$r = abs($h[$i][$i - 1]);
			}

			next if ($c == 0.0 || $r == 0.0);

			my $g = $r / BASE;
			my $f = 1.0;
			my $s = $c + $r;
			while ( $c < $g )
			{
				$f = $f * BASE;
				$c = $c * BASESQR;
			}

			$g = $r * BASE;
			while ($c >= $g)
			{
				$f = $f / BASE;
				$c = $c / BASESQR;
			}

			if (($c + $r) < 0.95 * $s * $f)
			{
				$g = 1.0 / $f;
				$noconv = 1;

				#C Generic code.
				#C   do $j=1,$n
				#C	 $h($i,$j)=$h($i,$j)*$g
				#C   enddo
				#C   do $j=1,$n
				#C	 $h($j,$i)=$h($j,$i)*$f
				#C   enddo
				#C begin specific code. Touch only non-zero elements of companion.
				if ($i == 1)
				{
					$h[1][$n] *= $g;
				}
				else
				{
					$h[$i][$i - 1] *= $g;
					$h[$i][$n] *= $g;
				}
				if ($i != $n)
				{
					$h[$i + 1][$i] *= $f;
				}
				else
				{
					for my $j (1 .. $n)
					{
						$h[$j][$i] *= $f;
					}
				}
			}
		}	# for $i
	}	# while $noconv

	return \@h;
}

#
# @roots = hqr_eigen_hessenberg($matrix_ref)
#
# Finds the eigenvalues of a real upper Hessenberg matrix,
# H, stored in the array $h(1:n,1:n).  Returns a list
# of real and/or complex numbers.
#
sub hqr_eigen_hessenberg($)
{
	my $ref = shift;
	my @h   = @$ref;
	my $n   = $#h;

	#
	#
	# Eigenvalue Computation by the Householder QR method for the Real Hessenberg matrix.
	# Perl code translated by Nick Ing-Simmons <Nick@Ing-Simmons.net>
	# from FORTRAN code by Hiroshi Murakami.
	# The fortran code is based on the Algol code "hqr" from the paper:
	#   "The QR Algorithm for Real Hessenberg Matrices"
	#   by R. S. Martin, G. Peters and J. H. Wilkinson,
	#   Numer. Math. 14, 219-231(1970).
	#
	#
	my($p, $q, $r);
	my($w, $x, $y);
	my($s, $z );
	my $t = 0.0;

	my @w;

	ROOT:
	while ($n > 0)
	{
		my $its = 0;
		my $na  = $n - 1;

		while ($its < $MAX_ITERATIONS)
		{
			#
			# Look for single small sub-diagonal element;
			#
			my $l;
			for ($l = $n; $l >= 2; $l--)
			{
				last if ( abs( $h[$l][ $l - 1 ] ) <= $epsilon *
					( abs( $h[ $l - 1 ][ $l - 1 ] ) + abs( $h[$l][$l] ) ) );
			}

			$x = $h[$n][$n];

			if ($l == $n)
			{
				#
				# One (real) root found.
				#
				push @w, $x + $t;
				$n--;
				next ROOT;
			}

			$y = $h[$na][$na];
			$w = $h[$n][$na] * $h[$na][$n];

			if ($l == $na)
			{
				$p = ( $y - $x ) / 2;
				$q = $p * $p + $w;
				$y = sqrt( abs($q) );
				$x = $x + $t;

				if ($q > 0.0)
				{
					#
					# Real pair.
					#
					$y = -$y if ( $p < 0.0 );
					$y += $p;
					push @w, $x - $w / $y;
					push @w, $x + $y;
				}
				else
				{
					#
					# Complex or twin pair.
					#
					push @w, $x + $p - $y * i;
					push @w, $x + $p + $y * i;
				}

				$n -= 2;
				next ROOT;
			}

			croak "Too many iterations ($its) at n=$n\n" if ($its == $MAX_ITERATIONS);

			if ($its && $its % 10 == 0)
			{
				#
				# Form exceptional shift.
				#
				carp "exceptional shift \@ $its" if ($debug_poly);

				$t = $t + $x;
				for (my $i = 1; $i <= $n; $i++)
				{
					$h[$i][$i] -= $x;
				}

				$s = abs($h[$n][$na]) + abs($h[$na][$n - 2]);
				$y = 0.75 * $s;
				$x = $y;
				$w = -0.4375 * $s * $s;
			}

			$its++;

			#
			# Look for two consecutive small sub-diagonal elements.
			#
			my $m;
			for ($m = $n - 2 ; $m >= $l ; $m--)
			{
				$z = $h[$m][$m];
				$r = $x - $z;
				$s = $y - $z;
				$p = ($r * $s - $w) / $h[$m + 1][$m] + $h[$m][$m + 1];
				$q = $h[$m + 1][$m + 1] - $z - $r - $s;
				$r = $h[$m + 2][$m + 1];

				$s = abs($p) + abs($q) + abs($r);
				$p = $p / $s;
				$q = $q / $s;
				$r = $r / $s;

				last if ($m == $l);
				last if (
					abs($h[$m][$m - 1]) * (abs($q) + abs($r)) <=
					$epsilon * abs($p) * (
						  abs($h[$m - 1][$m - 1]) +
						  abs($z) +
						  abs($h[$m + 1][$m + 1])
					)
				  );
			}

			for (my $i = $m + 2; $i <= $n; $i++)
			{
				$h[$i][$i - 2] = 0.0;
			}
			for (my $i = $m + 3; $i <= $n; $i++)
			{
				$h[$i][$i - 3] = 0.0;
			}

			#
			# Double QR step involving rows $l to $n and
			# columns $m to $n.
			#
			for (my $k = $m; $k <= $na; $k++)
			{
				my $notlast = ($k != $na);
				if ($k != $m)
				{
					$p = $h[$k][$k - 1];
					$q = $h[$k + 1][$k - 1];
					$r = ($notlast)? $h[$k + 2][$k - 1]: 0.0;

					$x = abs($p) + abs($q) + abs($r);
					next if ( $x == 0.0 );

					$p = $p / $x;
					$q = $q / $x;
					$r = $r / $x;
				}

				$s = sqrt($p * $p + $q * $q + $r * $r);
				$s = -$s if ($p < 0.0);

				if ($k != $m)
				{
					$h[$k][ $k - 1 ] = -$s * $x;
				}
				elsif ($l != $m)
				{
					$h[$k][$k - 1] = -$h[$k][$k - 1];
				}

				$p += $s;
				$x = $p / $s;
				$y = $q / $s;
				$z = $r / $s;
				$q /= $p;
				$r /= $p;

				#
				# Row modification.
				#
				for (my $j = $k; $j <= $n; $j++)
				{
					$p = $h[$k][$j] + $q * $h[$k + 1][$j];

					if ($notlast)
					{
						$p = $p + $r * $h[ $k + 2 ][$j];
						$h[ $k + 2 ][$j] -= $p * $z;
					}

					$h[ $k + 1 ][$j] -= $p * $y;
					$h[$k][$j] -= $p * $x;
				}

				my $j = $k + 3;
				$j = $n if $j > $n;

				#
				# Column modification.
				#
				for (my $i = $l; $i <= $j; $i++)
				{
					$p = $x * $h[$i][$k] + $y * $h[$i][$k + 1];

					if ($notlast)
					{
						$p += $z * $h[$i][$k + 2];
						$h[$i][$k + 2] -= $p * $r;
					}

					$h[$i][$k + 1] -= $p * $q;
					$h[$i][$k] -= $p;
				}
			}	# for $k
		}	# while $its
	}	# while $n
	return @w;
}

#
# A debugging aid.
#
sub show_matrix
{
	my $a = shift;
	my @rows = @$a;
	shift (@rows);
	foreach my $row (@rows)
	{
		warn join ( ',', @{$row}[ 1 .. @$a - 1 ] ) . "\n";
	}
}

#
# @x = poly_roots(@coefficients)
#
# Coefficients are fed in highest degree first.  Equation 5x**5 + 4x**4 + 2x + 8
# would be fed in with @x = poly_roots(5, 4, 0, 0, 2, 8);
#
sub poly_roots(@)
{
	my(@coefficients) = @_;
	my(@x) = ();

	#
	# Check for zero coefficients in the higher-powered terms.
	#
	shift @coefficients while ($#coefficients and abs($coefficients[0]) < $epsilon);

	if (@coefficients == 0)
	{
		carp "All coefficients are zero\n";
		return (0);
	}

	#
	# How about zero coefficients in the low terms?
	#
	while (@coefficients and abs($coefficients[$#coefficients]) < $epsilon)
	{
		push @x, 0;
		pop @coefficients;
	}

	#
	# If the remaining polynomial is a quintic or higher, or
	# if use_hessenberg is set, continue with the matrix
	# calculation.
	#
	if (($use_hessenberg and $#coefficients > 0) or $#coefficients > 4)
	{
		my $matrix_ref = build_companion(@coefficients);

		if ($debug_poly)
		{
			carp "Balanced Companion:\n";
			show_matrix($matrix_ref);
		}

		#
		# QR iterations from the matrix.
		#
		push @x, hqr_eigen_hessenberg($matrix_ref);
	}
	elsif ($#coefficients == 4)
	{
		push @x, quartic_roots(@coefficients);
	}
	elsif ($#coefficients == 3)
	{
		push @x, cubic_roots(@coefficients);
	}
	elsif ($#coefficients == 2)
	{
		push @x, quadratic_roots(@coefficients);
	}
	elsif ($#coefficients == 1)
	{
		push @x, linear_roots(@coefficients);
	}

	return @x;
}
sub GetLeastSquaresPlane {
    my ($helix, $CoG) = @_;
    my ($dSumXX, $dSumYY, $dSumZZ, $dSumXY, $dSumXZ, $dSumYZ);
    my ($x1, $y1, $z1, $x2, $y2, $z2, $group, $DATA, $i, $n);
    my ($o, $p, $q, $root, $dnorm, $a, $b, $c, @tmp, $smallest);

    $DATA = ();
    @tmp  = keys %{ $helix->{ATOMS} };
    $n = $#tmp + 1;
    for $group (1 .. $n) {
	for $i ("XCOORD", "YCOORD", "ZCOORD") { # calculate displacement from center of geometry
	    $DATA->[($group - 1)]{$i} = $helix->{ATOMS}{$tmp[$group - 1]}{$i} - $CoG->{$i};
	}
    }

    if ($n == 3) {
	$x1 = $DATA->[1]{XCOORD} - $DATA->[0]{XCOORD};
	$y1 = $DATA->[1]{YCOORD} - $DATA->[0]{YCOORD};
	$z1 = $DATA->[1]{ZCOORD} - $DATA->[0]{ZCOORD};
	$x2 = $DATA->[2]{XCOORD} - $DATA->[1]{XCOORD};
	$y2 = $DATA->[2]{YCOORD} - $DATA->[1]{YCOORD};
	$z2 = $DATA->[2]{ZCOORD} - $DATA->[1]{ZCOORD};
	
	$a = $y1 * $z2 - $z1 * $y2;
	$b = $z1 * $x2 - $x1 * $z2;
	$c = $x1 * $y2 - $y1 * $x2;	
    } else {
	# Calc co-variance
	$dSumXX = $dSumYY = $dSumZZ = $dSumXY = $dSumXZ = $dSumYZ = 0;
	for $i (0 .. ($n - 1)) {
	    $dSumXX += $DATA->[$i]{XCOORD} * $DATA->[$i]{XCOORD};
	    $dSumYY += $DATA->[$i]{YCOORD} * $DATA->[$i]{YCOORD};
	    $dSumZZ += $DATA->[$i]{ZCOORD} * $DATA->[$i]{ZCOORD};

	    $dSumXY += $DATA->[$i]{XCOORD} * $DATA->[$i]{YCOORD};
	    $dSumXZ += $DATA->[$i]{XCOORD} * $DATA->[$i]{ZCOORD};
	    $dSumYZ += $DATA->[$i]{YCOORD} * $DATA->[$i]{ZCOORD};
	}
	
	# Calc coeff. for -l^3 + o * l^2 + p * l + q = 0
	$o = $dSumXX + $dSumYY + $dSumZZ;
	$p = $dSumXY**2 + $dSumXZ**2 + $dSumYZ**2 - 
	    ($dSumXX * $dSumYY + $dSumXX * $dSumZZ + $dSumYY * $dSumZZ);
	$q = $dSumXX * $dSumYY * $dSumZZ + 2.0 * $dSumXY * $dSumXZ * $dSumYZ -
	    ($dSumXX * $dSumYZ * $dSumYZ + $dSumYY * $dSumXZ * $dSumXZ + $dSumZZ * $dSumXY * $dSumXY);

	# solve the cubic equation
	@tmp = cubic_roots(-1.0, $o, $p, $q);
	$root = 0;
	for (@tmp) {
	    if ($_ != 0) {
		$root = $_ if ($_ < $root || ! $root);
	    }
	}

	#calc determinantes
	$a = ($dSumYY - $root) * $dSumXZ - $dSumXY * $dSumYZ;
	$b = ($dSumXX - $root) * $dSumYZ - $dSumXY * $dSumXZ;
	$c = $dSumXY * $dSumXY - ($dSumYY - $root) * ($dSumXX - $root);
    }

    $dnorm = 1/sqrt($a*$a + $b*$b + $c*$c); # normalize
    $a *= $dnorm;
    $b *= $dnorm;
    $c *= $dnorm;

    $helix->{NORMAL}{a} = $a;
    $helix->{NORMAL}{b} = $b;
    $helix->{NORMAL}{c} = $c;
}

1;
__END__

=head1 NAME

Math::Polynomial::Solve - Find the roots of polynomial equations.

=head1 SYNOPSIS

  use Math::Complex;  # The roots may be complex numbers.
  use Math::Polynomial::Solve qw(poly_roots);

  my @x = poly_roots(@coefficients);

or

  use Math::Complex;  # The roots may be complex numbers.
  use Math::Polynomial::Solve qw(poly_roots get_hessenberg set_hessenberg);

  #
  # Force the use of the matrix method.
  #
  set_hessenberg(1);
  my @x = poly_roots(@coefficients);

or

  use Math::Complex;  # The roots may be complex numbers.
  use Math::Polynomial::Solve
	qw(linear_roots quadratic_roots cubic_roots quartic_roots);

  # Find the roots of ax + b
  my @x1 = linear_roots($a, $b);

  # Find the roots of ax**2 + bx +c
  my @x2 = quadratic_roots($a, $b, $c);

  # Find the roots of ax**3 + bx**2 +cx + d
  my @x3 = cubic_roots($a, $b, $c, $d);

  # Find the roots of ax**4 + bx**3 +cx**2 + dx + e
  my @x4 = quartic_roots($a, $b, $c, $d, $e);

=head1 DESCRIPTION

This package supplies a set of functions that find the roots of
polynomials. Polynomials up to the quartic may be solved directly by
numerical formulae. Polynomials of fifth and higher powers will be
solved by an iterative method, as there are no general solutions
for fifth and higher powers.

The linear, quadratic, cubic, and quartic *_roots() functions all expect
to have a non-zero value for the $a term.

If the constant term is zero then the first value returned in the list
of answers will always be zero, for all functions.

=head2 set_hessenberg()

Sets or removes the condition that forces the use of the Hessenberg matrix
regardless of the polynomial's degree.  A non-zero argument forces the
use of the matrix method, a zero removes it.

=head2 get_hessenberg()

Returns 1 or 0 depending upon whether the Hessenberg matrix method is
always in use or not.

=head2 poly_roots()

A generic function that may call one of the other root-finding functions,
or a polynomial solving method using a Hessenberg matrix, depending on the
degree of the polynomial.  You may force it to use the matrix method regardless
of the degree of the polynomial by calling C<set_hessenberg(1)>.
Otherwise it will use the specialized root functions for polynomials of
degree 1 to 4.

Unlike the other root-finding functions, it will check for coefficients
of zero for the highest power, and 'step down' the degree of the
polynomial to the appropriate case. Additionally, it will check for
coefficients of zero for the lowest power terms, and add zeros to its
root list before calling one of the root-finding functions. Thus
it is possible to solve a polynomial of degree higher than 4 without
using the matrix method, as long as it meets these rather specialized
conditions.

=head2 linear_roots()

Here for completeness's sake more than anything else. Returns the
solution for

    ax + b = 0

by returning C<-b/a>. This may be in either a scalar or an array context.

=head2 quadratic_roots()

Gives the roots of the quadratic equation

    ax**2 + bx + c = 0

using the well-known quadratic formula. Returns a two-element list.

=head2 cubic_roots()

Gives the roots of the cubic equation

    ax**3 + bx**2 + cx + d = 0

by the method described by R. W. D. Nickalls (see the Acknowledgments
section below). Returns a three-element list. The first element will
always be real. The next two values will either be both real or both
complex numbers.

=head2 quartic_roots()

Gives the roots of the quartic equation

    ax**4 + bx**3 + cx**2 + dx + e = 0

using Ferrari's method (see the Acknowledgments section below). Returns
a four-element list. The first two elements will be either
both real or both complex. The next two elements will also be alike in
type.

=head2 EXPORT

There are no default exports. The functions may be named in an export list.

=head1 Acknowledgments

=head2 The cubic

The cubic is solved by the method described by R. W. D. Nickalls, "A New
Approach to solving the cubic: Cardan's solution revealed," The
Mathematical Gazette, 77, 354-359, 1993. This article is available on
several different web sites, including
L<http://www.2dcurves.com/cubic/cubic.html> and
L<http://www.m-a.org.uk/resources/periodicals/online_articles_keyword_index/>.
There is also a nice discussion of his paper at
L<http://www.sosmath.com/algebra/factor/fac111/fac111.html>.

The paper is also downloadable from CPAN by choosing the module
Math-Polynomial-Solve-Doc.  This module is a documents-only package.

Dr. Nickalls was kind enough to send me his article, with notes and
revisions, and directed me to a Matlab script that was based on that
article, written by Herman Bruyninckx, of the Dept. Mechanical Eng.,
Div. PMA, Katholieke Universiteit Leuven, Belgium. This function is an
almost direct translation of that script, and I owe Herman Bruyninckx
for creating it in the first place. 

Dick Nickalls, dicknickalls@compuserve.com

Herman Bruyninckx, Herman.Bruyninckx@mech.kuleuven.ac.be,
has web page at L<http://www.mech.kuleuven.ac.be/~bruyninc>.
His matlab cubic solver is at
L<http://people.mech.kuleuven.ac.be/~bruyninc/matlab/cubic.m>.

=head2 The quartic

The method for quartic solution is Ferrari's, as described in the web
page Karl's Calculus Tutor at L<http://www.karlscalculus.org/quartic.html>.
I also made use of some short cuts mentioned in web page Ask Dr. Math FAQ,
at L<http://forum.swarthmore.edu/dr.math/faq/faq.cubic.equations.html>.

=head2 Quintic and higher.

Back when this module could only solve polynomials of degrees 1 through 4,
Matz Kindahl, the author of Math::Polynomial, suggested the poly_roots()
function. Later on, Nick Ing-Simmons, who was working on a perl binding
to the GNU Scientific Library, sent a perl translation of Hiroshi
Murakami's Fortran implementation of the QR Hessenberg algorithm, and it
fit very well into the poly_roots() function. Quintics and higher degree
polynomials can now be solved, albeit through numeric analysis methods.

Hiroshi Murakami's Fortran routines were at
L<http://netlib.bell-labs.com/netlib/>, but do not seem to be available
from that source anymore.

He referenced the following articles:

R. S. Martin, G. Peters and J. H. Wilkinson, "The QR Algorithm for Real Hessenberg
Matrices", Numer. Math. 14, 219-231(1970).

B. N. Parlett and C. Reinsch, "Balancing a Matrix for Calculation of Eigenvalues
and Eigenvectors", Numer. Math. 13, 293-304(1969).

Alan Edelman and H. Murakami, "Polynomial Roots from Companion Matrix
Eigenvalues", Math. Comp., v64,#210, pp.763-776(1995).

For starting out, you may want to read

Numerical Recipes in C, by William Press, Brian P. Flannery, Saul A. Teukolsky,
and William T. Vetterling, Cambridge University Press.

=head1 SEE ALSO

Forsythe, George E., Michael A. Malcolm, and Cleve B. Moler (1977),
Computer Methods for Mathematical Computations, Prentice-Hall.

=head1 AUTHOR

John M. Gamble may be found at B<jgamble@ripco.com>

=cut

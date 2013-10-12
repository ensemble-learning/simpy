package Math::Quaternion;

use 5.004;
use strict;
use warnings;
use Carp;
use Math::Trig; # What?!? Where's acos()? You can't have cos and not acos!

require Exporter;
use overload
	'+' => \&plus,
	'-' => \&minus,
	'bool' => sub { 1; }, # So we can do if ($foo=Math::Quaternion->new) { .. }
	'""' => \&stringify,
	'*' => \&multiply,
	'~' => \&conjugate,
	'abs' => \&modulus,
	'neg' => \&negate,
	'**' => \&power,
	'exp' => \&exp,
	'log' => \&log,
	;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use Math::Quaternion ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.


our %EXPORT_TAGS = ( 'all' => [ qw(
unit
conjugate
inverse
normalize
modulus
isreal
multiply
dot
plus
minus
power
negate
squarednorm
scale
rotation
rotation_angle
rotation_axis
rotate_vector
matrix4x4
matrix3x3
matrix4x4andinverse
stringify
slerp
exp
log
) ],
);

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);

our $VERSION = '0.02';


# Preloaded methods go here.


# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Math::Quaternion - Perl class to represent quaternions

=head1 SYNOPSIS

 use Math::Quaternion qw(slerp);
 my $q = Math::Quaternion->new;  # Make a new unit quaternion
 
 # Make a rotation about the axis (0,1,0)
 my $q2 = Math::Quaternion->new({axis=>[0,1,0],angle=>0.1});
 my @v = (1,2,3); # A vector.
 my @vrotated = $q2->rotate_vector(@v); # Rotate @v about (0,1,0).
 
 my $q3 = Math::Quaternion::rotation(0.7,2,1,4); # A different rotation.
 my $q4 = slerp($q2,$q3,0.5);                   # Interpolated rotation.
 my @vinterp = $q4->rotate_vector(@v);


=head1 DESCRIPTION

This package lets you create and manipulate quaternions. A
quaternion is a mathematical object developed as a kind of
generalization of complex numbers, usually represented by an array
of four real numbers, and is often used to represent rotations in
three-dimensional space.

See, for example, L<http://mathworld.wolfram.com/Quaternion.html> for
more details on the mathematics of quaternions.

Quaternions can be added, subtracted, and scaled just like complex
numbers or vectors -- they can also be multiplied, but quaternion
multiplication DOES NOT COMMUTE. That is to say, if you have
quaternions $q1 and $q2, then in general $q1*$q2 != $q2*$q1. This is
related to their use in representing rotations, which also do not
commute.

If you just want to represent rotations and don't care about the
internal mathematical details, this should be all you need to know:

All quaternions have a quantity called the "norm",  similar to the
length of a vector. A quaternion with norm equal to 1 is called a
"unit quaternion". All quaternions which represent rotations are
unit quaternions.

If you call new() without any arguments, it will give you a unit
quaternion which represents no rotation:
 
   $q = Math::Quaternion->new;

You can make a quaternion which represents a rotation of a given
angle (in radians) about a given axis:

   $qrot = Math::Quaternion->new({ axis => 0.1, angle => [ 2,3,4]});

Say you have two rotations, $q1 and $q2, and you want to make a
quaternion representing a rotation of $q1 followed by $q2. Then, you
do:

  $q3 = $q2 * $q1;   # Rotate by $q1, followed by $q2.

Remember that this is NOT the same as $q1 * $q2, which will reverse
the order of the rotations.

If you perform many iterated quaternion operations, the result may
not quite be a unit quaternion due to numerical inaccuracies. You
can make sure any quaternion has unit length, by doing:

  $unitquat = $anyquat->normalize;

If you have a rotation quaternion, and you want to find the 3x3
matrix which represents the corresponding rotation, then:

  @matrix = $q->matrix3x3;

Similarly, you can generate a 4x4 matrix of the sort you'd pass to
OpenGL:

  @glmatrix = $q->matrix4x4;

If you have a vector representing a direction, and you want to
rotate the vector by a quaternion $q:

  my @vector = (0,0,1);  # Vector pointing in the Z direction. 
  
  my @newvec = $q->rotate_vector(@vector); # New direction.

Say you're using quaternions to represent the orientation of a
camera, and you have two quaternions: one to represent a
starting orientation, and another to represent a finishing
position. If you want to find all the quaternions representing
the orientations in between, allowing your camera to move
smoothly from start to finish, use the slerp() routine:

  use Math::Quaternion qw(slerp);
  
  my ($qstart, $qend) = ... ;
   
  # Set $tween to 9 points between start and end, exclusive.
  
  for my $t (1..9) {
    my $tween = slerp($qstart,$qend,0.1*$t); 
    ...
  }


=head1 METHODS

=over 1

=item B<new>

 my $q = Math::Quaternion->new;          # Make a new unit quaternion.
 my $q2 = Math::Quaternion->new(1,2,3,4);# Make a specific quaternion.
 my $q3 = Math::Quaternion->new($q2);    # Copy an existing quaternion.
 my $q4 = Math::Quaternion->new(5.6);    # Make the quaternion (5.6,0,0,0)
 my $q5 = Math::Quaternion->new(7,8,9);  # Make the quaternion (0,7,8,9)
  
 my $q6 = Math::Quaternion->new({ # Make a quaternion corresponding
       axis => [ 1,2,3],          # to a rotation of 0.2 radians
       angle => 0.2,              # about the vector (1,2,3).
 });
 
 my $q7 = Math::Quaternion->new({ # Make a quaternion which would
       'v1' => [ 0,1,2],            # rotate the vector (0,1,2) onto
       'v2' => [ -1,2,0],           # the vector (-1,2,0).
 });

If no parameters are given, a unit quaternion is returned.  If one
non-reference parameter is given, a "scalar" quaternion is returned.
If one parameter is given and it is a reference to a quaternion or
an array of four numbers, the corresponding quaternion object is
returned.  If three parameters are given, a "vector" quaternion is
returned.  If four parameters are given, the corresponding
quaternion is returned.

Rotation quaternions may also be created by passing a hashref with
the axis and angle of rotation, or by specifying two vectors
specifying start and finish directions. Bear in mind that the latter
method will take the shortest path between the two vectors, ignoring
the "roll" angle.

=cut

sub new {
	my $class = shift;

	my $arr=undef;

	if (0==@_) {
		# No arguments, default to unit quaternion.
		$arr = [ 1,0,0,0];
	} elsif (1==@_) {
		# One argument: if it's not a reference, construct
		# a "scalar quaternion" (x 0 0 0).
		my $arg = $_[0];
		my $reftype = ref($arg);

		if (!$reftype) {
			$arr = [ $arg,0,0,0];
		} else {
			# We've been passed a reference. If it's an array
			# ref, then construct a quaternion out of the
			# corresponding array.
			if ("ARRAY" eq $reftype) {
				return Math::Quaternion->new(@$arg);
			} elsif ("Math::Quaternion" eq $reftype) {
				# If it's a reference to another quaternion,
				# copy it.
				return Math::Quaternion->new(@$arg);
			} elsif ("HASH" eq $reftype) {
				# Hashref.
				my %hash = %$arg;
				if (defined($hash{'axis'})) {
					# Construct a rotation.
					return rotation(
						$hash{'angle'},
						@{$hash{'axis'}}
					);
				} elsif (defined($hash{'v2'})) {
					return rotation(
						$hash{'v1'},$hash{'v2'}
					);
				}
			}
			croak("Don't understand arguments to new()");

		}
	} elsif (3==@_) {
		# Three arguments: construct a quaternion to represent
		# the corresponding vector.
		$arr = [ 0, @_[0,1,2] ];
	} elsif (4==@_) {
		# Four arguments: just slot the numbers right in.
		$arr = [ @_[0,1,2,3] ];
	} else {
		croak("Don't understand arguments passed to new()");
	}
		

	bless $arr, $class;

}

=item B<unit>

Returns a unit quaternion.

 my $u = Math::Quaternion->unit; # Returns the quaternion (1,0,0,0).

=cut

sub unit {
	my $class = shift;

	bless [ 1,0,0,0 ], $class;
}

=item B<conjugate>

Returns the conjugate of its argument.

 my $q = Math::Quaternion->new(1,2,3,4);
 my $p = $q->conjugate;              # (1,-2,-3,-4)

=cut

sub conjugate {
	my $q=shift;

	return Math::Quaternion->new(
		  $q->[0],
		- $q->[1],
		- $q->[2],
		- $q->[3],
	);
}

=item B<inverse>

Returns the inverse of its argument.

 my $q = Math::Quaternion->new(1,2,3,4);
 my $qi = $q->inverse;

=cut

sub inverse {
	my $q = shift;

	return scale(conjugate($q),1.0/squarednorm($q));

}


=item B<normalize>

Returns its argument, normalized to unit norm.

  my $q = Math::Quaternion->new(1,2,3,4);
  my $qn = $q->normalize;

=cut

sub normalize {
	my $q = shift;
	return scale($q,1.0/sqrt(squarednorm($q)));
}

=item B<modulus>

Returns the modulus of its argument, defined as the 
square root of the scalar obtained by multiplying the quaternion
by its conjugate.

 my $q = Math::Quaternion->new(1,2,3,4);
 print $q->modulus;

=cut

sub modulus {
	my $q = shift;
	return sqrt(squarednorm($q));
}

=item B<isreal>

Returns 1 if the given quaternion is real ,ie has no quaternion
part, or else 0.

 my $q1 = Math::Quaternion->new(1,2,3,4);
 my $q2 = Math::Quaternion->new(5,0,0,0);
 print $q1->isreal; # 1;
 print $q2->isreal; # 0;

=cut

sub isreal {
	my $q = shift;
	my ($q0,$q1,$q2,$q3)=@$q;

	if ( (0.0==$q1) && (0.0==$q2) && (0.0==$q3) ) {
		return 1;
	} else {
		return 0;
	}
}

=item B<multiply>

Performs a quaternion multiplication of its two arguments.
If one of the arguments is a scalar, then performs a scalar
multiplication instead.

 my $q1 = Math::Quaternion->new(1,2,3,4);
 my $q2 = Math::Quaternion->new(5,6,7,8);
 my $q3 = Math::Quaternion::multiply($q1,$q2);         # (-60 12 30 24)
 my $q4 = Math::Quaternion::multiply($q1,$q1->inverse); # (1 0 0 0) 

=cut

sub multiply {
	my ($a,$b,$reversed) = @_;
	($a,$b) = ($b,$a) if $reversed;

	if (!ref $a) { return scale($b,$a); }
	if (!ref $b) { return scale($a,$b); }

	my $q = new Math::Quaternion;

	$q->[0] = $a->[0] * $b->[0] 
		- $a->[1]*$b->[1]
		- $a->[2]*$b->[2]
		- $a->[3]*$b->[3];
	
	$q->[1] = $a->[0] * $b->[1]
		+ $b->[0] * $a->[1]
		+ $a->[2] * $b->[3] - $a->[3] * $b->[2];

	$q->[2] = $a->[0] * $b->[2]
		+ $b->[0] * $a->[2]
		+ $a->[3] * $b->[1] - $a->[1] * $b->[3];

	$q->[3] = $a->[0] * $b->[3]
		+ $b->[0] * $a->[3]
		+ $a->[1] * $b->[2] - $a->[2] * $b->[1];
	return $q;
}

=item B<dot>

Returns the dot product of two quaternions.

 my $q1=Math::Quaternion->new(1,2,3,4);
 my $q2=Math::Quaternion->new(2,4,5,6);
 my $q3 = Math::Quaternion::dot($q1,$q2);

=cut

sub dot {
	my ($q1,$q2) = @_;
	my ($a0,$a1,$a2,$a3) = @$q1;
	my ($b0,$b1,$b2,$b3) = @$q2;
	return $a0*$b0 + $a1*$b1 + $a2*$b2 + $a3*$b3 ;
}

=item B<plus>

Performs a quaternion addition of its two arguments.

 my $q1 = Math::Quaternion->new(1,2,3,4);
 my $q2 = Math::Quaternion->new(5,6,7,8);
 my $q3 = Math::Quaternion::plus($q1,$q2);         # (6 8 10 12)

=cut


sub plus {
	my ($a,$b,$reversed)=@_;
	my $q = Math::Quaternion->new(
		$a->[0] + $b->[0],
		$a->[1] + $b->[1],
		$a->[2] + $b->[2],
		$a->[3] + $b->[3],
	);

	return $q;

}

=item B<minus>

Performs a quaternion subtraction of its two arguments.

 my $q1 = Math::Quaternion->new(1,2,3,4);
 my $q2 = Math::Quaternion->new(5,6,7,8);
 my $q3 = Math::Quaternion::minus($q1,$q2);         # (-4 -4 -4 -4)

=cut

sub minus {
	my ($a,$b,$reversed)=@_;
	($a,$b) = ($b,$a) if $reversed;
	my $q = Math::Quaternion->new(
		$a->[0] - $b->[0],
		$a->[1] - $b->[1],
		$a->[2] - $b->[2],
		$a->[3] - $b->[3],
	);

	return $q;

}

=item B<power>

Raise a quaternion to a scalar or quaternion power.

 my $q1 = Math::Quaternion->new(1,2,3,4);
 my $q2 = Math::Quaternion::power($q1,4);     # ( 668 -224 -336 -448 )
 my $q3 = $q1->power(4);                # ( 668 -224 -336 -448 )
 my $q4 = $q1**(-1);			 # Same as $q1->inverse

 use Math::Trig;
 my $q5 = exp(1)**( Math::Quaternion->new(pi,0,0) ); # approx (-1 0 0 0)

=cut

sub power {
	my ($a,$b,$reversed)=@_;
	($a,$b) = ($b,$a) if $reversed;

	if (ref $a) {
		$a = Math::Quaternion->new($a);
	}

	if (ref $b) {
		# For quaternion^quaternion, use exp and log.
		return Math::Quaternion::exp(Math::Quaternion::multiply($b,Math::Quaternion::log($a)));
	}

	# For quat raised to a scalar power, do it manually.

	my ($a0,$a1,$a2,$a3) = @$a;

	my $s = sqrt($a->squarednorm);
	my $theta = Math::Trig::acos($a0/$s);
	my $vecmod = sqrt($a1*$a1+$a2*$a2+$a3*$a3);
	my $stob = ($s**$b);
	my $coeff = $stob/$vecmod*sin($b*$theta);
	
	my $u1 = $a1*$coeff;
	my $u2 = $a2*$coeff;
	my $u3 = $a3*$coeff;


	return Math::Quaternion->new(
		$stob * cos($b*$theta), $u1,$u2,$u3
	);
	

}

=item B<negate>

Negates the given quaternion.

 my $q = Math::Quaternion->new(1,2,3,4);
 my $q1 = $q->negate;             # (-1,-2,-3,-4)

=cut

sub negate {

	my $q = shift;
	return  Math::Quaternion->new(
		-($q->[0]),
		-($q->[1]),
		-($q->[2]),
		-($q->[3]),
	);

}


=item B<squarednorm>

Returns the squared norm of its argument.

 my $q1 = Math::Quaternion->new(1,2,3,4);
 my $sn = $q1->squarednorm;             # 30

=cut

sub squarednorm {
	my $q = shift;
	return    $q->[0]*$q->[0] 
		+ $q->[1]*$q->[1] 
		+ $q->[2]*$q->[2] 
		+ $q->[3]*$q->[3];

}

=item B<scale>

Performs a scalar multiplication of its two arguments.

 my $q = Math::Quaternion->new(1,2,3,4);
 my $qq = Math::Quaternion::scale($q,2);   # ( 2 4 6 8)
 my $qqq= $q->scale(3);                    # ( 3 6 9 12 )

=cut

sub scale {
	my ($q,$s)=@_;
	return Math::Quaternion->new(
		$q->[0] * $s,
		$q->[1] * $s,
		$q->[2] * $s,
		$q->[3] * $s
	);
}

=item B<rotation>


Generates a quaternion corresponding to a rotation.

If given three arguments, interprets them as an angle and the
three components of an axis vector.

 use Math::Trig;            # Define pi.  my $theta = pi/2;
 # Angle of rotation my $rotquat =
 Math::Quaternion::rotation($theta,0,0,1);
 
 # $rotquat now represents a rotation of 90 degrees about Z axis.
 
 my ($x,$y,$z) = (1,0,0);	# Unit vector in the X direction.
 my ($xx,$yy,$zz) = $rotquat->rotate_vector($x,$y,$z);
 
 # ($xx,$yy,$zz) is now ( 0, 1, 0), to within floating-point error.


rotation() can also be passed a scalar angle and a reference to
a vector (in either order), and will generate the corresponding
rotation quaternion.

 my @axis = (0,0,1);    # Rotate about Z axis
 $theta = pi/2;
 $rotquat = Math::Quaternion::rotation($theta,\@axis);


If the arguments to rotation() are both references, they are
interpreted as two vectors, and a quaternion is returned which
rotates the first vector onto the second.

 my @startvec = (0,1,0);  # Vector pointing north
 my @endvec   = (-1,0,0); # Vector pointing west
 $rotquat = Math::Quaternion::rotation(\@startvec,\@endvec);
 
 my @newvec = $rotquat->rotate_vector(@startvec); # Same as @endvec

=cut

sub rotation {
	my ($theta,$x,$y,$z);
	if (2==@_) {
		if (ref($_[0])) {
			if (ref($_[1])) {
				# Both args references to vectors
				my ($ax,$ay,$az)=@{$_[0]};
				my ($bx,$by,$bz)=@{$_[1]};
				# Find cross product. This is a vector
				# perpendicular to both argument vectors,
				# and is therefore the axis of rotation.
				$x = $ay*$bz-$az*$by;
				$y = $az*$bx-$ax*$bz;
				$z = $ax*$by-$ay*$bx;
				# find the dot product.
				my $dotprod = $ax*$bx+$ay*$by+$az*$bz;
				my $mod1 = sqrt($ax*$ax+$ay*$ay+$az*$az);
				my $mod2 = sqrt($bx*$bx+$by*$by+$bz*$bz);
				# Find the angle of rotation.
				$theta=Math::Trig::acos($dotprod/($mod1*$mod2));
			} else {
				# 0 is a ref, 1 is not.
				$theta = $_[1];
				($x,$y,$z)=@{$_[0]};
			}
		} else {
			if (ref($_[1])) {
				# 1 is a ref, 0 is not
				$theta = $_[0];
				($x,$y,$z)=@{$_[1]};
			} else {
			 croak("Math::Quaternion::rotation() passed 2 nonref args");
			}

		}
		
	} elsif (4==@_) {
		($theta,$x,$y,$z) = @_;
	} else {
		croak("Math::Quaternion::rotation() passed wrong no of arguments");
	}

	my $modulus = sqrt($x*$x+$y*$y+$z*$z); # Make it a unit vector
	$x /= $modulus;
	$y /= $modulus;
	$z /= $modulus;

	my $st = sin(0.5 * $theta);
	my $ct = cos(0.5 * $theta);

	return Math::Quaternion->new(
		$ct, $x * $st, $y * $st, $z * $st
	);
}

=item B<rotation_angle>

Returns the angle of rotation represented by the quaternion
argument.

 my $q = Math::Quaternion::rotation(0.1,2,3,4);
 my $theta = $q->rotation_angle; # Returns 0.1 .

=cut

sub rotation_angle {
	my $q = shift;
	return 2.0 * Math::Trig::acos($q->[0]);
}

=item B<rotation_axis>

Returns the unit vector representing the axis about which
rotations will be performed, for the rotation represented by the
quaternion argument.

 my $q = Math::Quaternion::rotation(0.1,1,1,0);
 my @v = $q->rotation_axis; # Returns (0.5*sqrt(2),0.5*sqrt(2),0)

=cut

sub rotation_axis {
	my $q = shift;
	my $theta = Math::Trig::acos($q->[0]);
	my $st = sin($theta);
	if (0==$st) { return (0,0,1); } # Rotation of angle zero about Z axis
	my ($x,$y,$z) = @{$q}[1,2,3];

	return ( $x/$st, $y/$st, $z/$st );
}




=item B<rotate_vector>

When called as a method on a rotation quaternion, uses this
quaternion to perform the corresponding rotation on the vector
argument.

 use Math::Trig;                     # Define pi.
 
 my $theta = pi/2;                   # Rotate 90 degrees
 
 my $rotquat = Math::Quaternion::rotation($theta,0,0,1); # about Z axis
 
 my ($x,$y,$z) = (1,0,0);	# Unit vector in the X direction.
 my ($xx,$yy,$zz) = $rotquat->rotate_vector($x,$y,$z)
 
 # ($xx,$yy,$zz) is now ( 0, 1, 0), to within floating-point error.

=cut


sub rotate_vector {
	my ($q,$x,$y,$z) = @_;

	my $p = Math::Quaternion->new($x,$y,$z);
	my $qq = multiply($q,multiply($p,inverse($q)));
	return  @{$qq}[1,2,3];
}


=item B<matrix4x4>

Takes one argument: a rotation quaternion.
Returns a 16-element array, equal to the OpenGL
matrix which represents the corresponding rotation.

 my $rotquat = Math::Quaternion::rotation($theta,@axis); # My rotation.
 my @m = $rotquat->matrix4x4;

=cut

sub matrix4x4 {
	my $q = shift;
	my ($w,$x,$y,$z) = @{$q};

	return (
		1 - 2*$y*$y - 2*$z*$z,
		2*$x*$y + 2*$w*$z,
		2*$x*$z - 2*$w*$y,
		0,

		2*$x*$y - 2*$w*$z,
		1 - 2*$x*$x - 2*$z*$z,
		2*$y*$z + 2*$w*$x,
		0,

		2*$x*$z + 2*$w*$y,
		2*$y*$z - 2*$w*$x,
		1 - 2*$x*$x - 2*$y*$y,
		0,

		0,
		0,
		0,
		1
	);
}

=item B<matrix3x3>

Takes one argument: a rotation quaternion.
Returns a 9-element array, equal to the 3x3
matrix which represents the corresponding rotation.

 my $rotquat = Math::Quaternion::rotation($theta,@axis); # My rotation.
 my @m = $rotquat->matrix3x3;

=cut

sub matrix3x3 {
	my $q = shift;
	my ($w,$x,$y,$z) = @{$q};

	return (
		1 - 2*$y*$y - 2*$z*$z,
		2*$x*$y + 2*$w*$z,
		2*$x*$z - 2*$w*$y,

		2*$x*$y - 2*$w*$z,
		1 - 2*$x*$x - 2*$z*$z,
		2*$y*$z + 2*$w*$x,

		2*$x*$z + 2*$w*$y,
		2*$y*$z - 2*$w*$x,
		1 - 2*$x*$x - 2*$y*$y,
	);
}

=item B<matrix4x4andinverse>

Similar to matrix4x4, but returnes a list of two array
references.  The first is a reference to the rotation matrix;
the second is a reference to its inverse.  This may be useful
when rendering sprites, since you can multiply by the rotation
matrix for the viewer position, perform some translations, and
then multiply by the inverse: any resulting rectangles drawn
will always face the viewer.


 my $rotquat = Math::Quaternion::rotation($theta,@axis); # My rotation.
 my ($matref,$invref) = $rotquat->matrix4x4andinverse;

=cut


sub matrix4x4andinverse {
	my $q = shift;
	my ($w,$x,$y,$z) = @{$q};
	my (@m,@mi);

	$mi[ 0] = $m[ 0] = 1 - 2*$y*$y - 2*$z*$z;
	$mi[ 4] = $m[ 1] = 2*$x*$y + 2*$w*$z;
	$mi[ 8] = $m[ 2] = 2*$x*$z - 2*$w*$y;
	$mi[12] = $m[ 3] = 0;

	$mi[ 1] = $m[ 4] = 2*$x*$y - 2*$w*$z;
	$mi[ 5] = $m[ 5] = 1 - 2*$x*$x - 2*$z*$z;
	$mi[ 9] = $m[ 6] = 2*$y*$z + 2*$w*$x;
	$mi[13] = $m[ 7] = 0;

	$mi[ 2] = $m[ 8] = 2*$x*$z + 2*$w*$y;
	$mi[ 6] = $m[ 9] = 2*$y*$z - 2*$w*$x;
	$mi[10] = $m[10] = 1 - 2*$x*$x - 2*$y*$y;
	$mi[14] = $m[11] = 0;

	$mi[ 3] = $m[12] = 0;
	$mi[ 7] = $m[13] = 0;
	$mi[11] = $m[14] = 0;
	$mi[15] = $m[15] = 1;

	return (\@m,\@mi);

}

=item B<stringify>

Returns a string representation of the quaternion. This is used
to overload the '""' operator, so that quaternions may be
freely interpolated in strings.

 my $q = Math::Quaternion->new(1,2,3,4);
 print $q->stringify;                # "( 1 2 3 4 )"
 print "$q";                         # "( 1 2 3 4 )"


=cut

sub stringify {
	my $self = shift;
	return "( ".join(" ",@$self)." )";
}

=item B<slerp>

Takes two quaternion arguments and one scalar; performs
spherical linear interpolation between the two quaternions. The
quaternion arguments are assumed to be unit quaternions, and the
scalar is assumed to lie between 0 and 1: a scalar argument of
zero will return the first quaternion argument, and a scalar
argument of one will return the second.

 use Math::Trig;
 my @axis = (0,0,1);
 my $rq1 = Math::Quaternion::rotation(pi/2,\@axis);   # 90  degs about Z
 my $rq2 = Math::Quaternion::rotation(pi,\@axis);     # 180 degs about Z
 
 my $interp = Math::Quaternion::slerp($rq1,$rq2,0.5); # 135 degs about Z

=cut

sub slerp {
	my ($q0,$q1,$t) = @_;

	my $dotprod = dot($q0,$q1);
	if ($dotprod<0) {
		# Reverse signs so we travel the short way round
		$dotprod = -$dotprod;
		$q1 = negate($q1);
	}

	my $theta = Math::Trig::acos($dotprod);

	if (abs($theta) < 1e-5) {
		# In the limit theta->0 , spherical interpolation is
		# approximated by linear interpolation, which also
		# avoids division-by-zero problems.

		return plus(scale($q0,(1-$t)) ,scale($q1,$t));

	}

	my $st = sin($theta);
	my $ist = 1.0/$st;

	my $q = plus(
		scale($q0,($ist * sin( (1-$t)*$theta ))),
		scale($q1,($ist*sin($t*$theta)))
	);

	
	return normalize($q);

}


=item B<exp>

Exponential operator e^q. Any quaternion q can be written as x+uy,
where x is a real number, and u is a unit pure quaternion.  Then,
exp(q) == exp(x) * ( cos(y) + u sin(y) ).

 my $q = Math::Quaternion->new(1,2,3,4);
 print Math::Quaternion::exp($q);

=cut

sub exp {
	my $q = shift;

	if (isreal($q)) {
		return Math::Quaternion->new(CORE::exp($q->[0]),0,0,0);
	}

	my ($q0,$q1,$q2,$q3)=@$q;

	my $y = sqrt($q1*$q1+$q2*$q2+$q3*$q3); # Length of pure-quat part.
	my ($ux,$uy,$uz) = ($q1/$y,$q2/$y,$q3/$y); # Unit vector.

	my $ex = CORE::exp($q0);
	my $exs = $ex*sin($y);

	return Math::Quaternion->new($ex*cos($y),$exs*$ux,$exs*$uy,$exs*$uz);
}

=item B<log>

Returns the logarithm of its argument. The logarithm of a negative
real quaternion can take any value of them form (log(-q0),u*pi) for
any unit vector u. In these cases, u is chosen to be (1,0,0).

 my $q = Math::Quaternion->new(1,2,3,4);
 print Math::Quaternion::log($q);

=cut

sub log {
	my $q = shift;

	if (ref $q) {
		if ("Math::Quaternion" ne ref $q) {
			$q = Math::Quaternion->new($q);
		}
	} else {
		$q = Math::Quaternion->new($q);
	}

	if (isreal($q)) {
		if ($q->[0] > 0) {
			return Math::Quaternion->new(CORE::log($q->[0]));
		} else {
			return Math::Quaternion->new(CORE::log(-($q->[0])),pi,0,0);
		}
	}

	my ($q0,$q1,$q2,$q3)=@$q;

	my $modq = sqrt($q0*$q0 + $q1*$q1 + $q2*$q2 + $q3*$q3);

	my $x = CORE::log($modq);
	my $qquatmod = sqrt($q1*$q1+$q2*$q2+$q3*$q3); # mod of quat part
	my $y = atan2($qquatmod,$q0);
	my $c = $y/$qquatmod;

	return Math::Quaternion->new($x,$c*$q1,$c*$q2,$c*$q3);
	
}



=back

=head1 AUTHOR

Jonathan Chin, E<lt>jon-quaternion.pm@earth.liE<gt>

=head1 ACKNOWLEDGEMENTS

Thanks to Rene Uittenbogaard for useful suggestions.

=head1 SEE ALSO

=over 4

=item L<http://mathworld.wolfram.com/Quaternion.html> 

=item L<http://sjbaker.org/steve/omniv/eulers_are_evil.html> 

=item Acts 12:4 

=back

=head1 COPYRIGHT AND LICENSE

Copyright 2003 by Jonathan Chin

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself. 

=cut

1;
__END__

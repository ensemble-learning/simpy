#                              -*- Mode: Perl -*- 
# Matrix.pm -- 
# ITIID           : $ITI$ $Header $__Header$
# Author          : Ulrich Pfeifer
# Created On      : Tue Oct 24 18:34:08 1995
# Last Modified By: Ulrich Pfeifer
# Last Modified On: Fri Sep 14 21:57:54 2001
# Language        : Perl
# Update Count    : 191
# Status          : Unknown, Use with caution!
#
# Copyright (C) 2001, Brian J. Watson <bjbrew@power.net>, all rights reserved.
# Copyright (C) 2001, Ulrich Pfeifer <pfeifer@wait.de>, all rights reserved.
# Copyright (C) 1995, Universität Dortmund, all rights reserved.
# Copyright (C)  2001, Matthew Brett <matthew.brett@mrc-cbu.cam.ac.uk>
#
# Permission to use this software is granted under the same
# restrictions as for Perl itself.
#
# Revision 0.3  2001/04/17 11:10:15
# Extensions from Brian Watson
#
# Revision 0.2  1996/07/10 17:48:14  pfeifer
# Fixes from Mike Beachy <beachy@chem.columbia.edu>
#
# Revision 0.1  1995/10/25  09:48:39  pfeifer
# Initial revision
#

=head1 NAME

Math::Matrix - Multiply and invert Matrices

=head1 SYNOPSIS

use Math::Matrix;

=head1 DESCRIPTION

The following methods are available:

=head2 new

Constructor arguments are a list of references to arrays of the same
length.  The arrays are copied. The method returns B<undef> in case of
error.

        $a = new Math::Matrix ([rand,rand,rand],
                               [rand,rand,rand],
                               [rand,rand,rand]);

If you call  C<new> as method, a zero filled matrix with identical deminsions is returned.

=head2 clone

You can clone a matrix by calling:

        $b = $a->clone;

=head2 size

You can determine the dimensions of a matrix by calling:

        ($m, $n) = $a->size;

=head2 concat

Concatenates two matrices of same row count. The result is a new
matrix or B<undef> in case of error.

        $b = new Math::Matrix ([rand],[rand],[rand]);
        $c = $a->concat($b);

=head2 transpose

Returns the transposed matrix. This is the matrix where colums and
rows of the argument matrix are swaped.

=head2 multiply

Multiplies two matrices where the length of the rows in the first
matrix is the same as the length of the columns in the second
matrix. Returns the product or B<undef> in case of error.

=head2 solve

Solves a equation system given by the matrix. The number of colums
must be greater than the number of rows. If variables are dependent
from each other, the second and all further of the dependent
coefficients are 0. This means the method can handle such systems. The
method returns a matrix containing the solutions in its columns or
B<undef> in case of error.

=head2 multiply_scalar

Multiplies a matrix and a scalar resulting in a matrix of the same
dimensions with each element scaled with the scalar.

  $a->multiply_scalar(2);  scale matrix by factor 2

=head2 add

Add two matrices of the same dimensions.

=head2 substract

Shorthand for C<add($other-E<gt>negative)>

=head2 equal

Decide if two matrices are equal.  Beware of rounding errors!

=head2 slice

Extract columns:

  a->slice(1,3,5);

=head2 determinant

Compute the determinant of a matrix.

=head2 dot_product

Compute the dot product of two vectors.

=head2 absolute

Compute the absolute value of a vector.

=head2 normalizing

Normalize a vector.

=head2 cross_product

Compute the cross-product of vectors.

=head2 print

Prints the matrix on STDOUT. If the method has additional parameters,
these are printed before the matrix is printed.

=head1 EXAMPLE

        use Math::Matrix;

        srand(time);
        $a = new Math::Matrix ([rand,rand,rand], 
                         [rand,rand,rand], 
                         [rand,rand,rand]);
        $x = new Math::Matrix ([rand,rand,rand]);
        $a->print("A\n");
        $E = $a->concat($x->transpose);
        $E->print("Equation system\n");
        $s = $E->solve;
        $s->print("Solutions s\n");
        $a->multiply($s)->print("A*s\n");

=head1 AUTHOR

Ulrich Pfeifer E<lt>F<pfeifer@ls6.informatik.uni-dortmund.de>E<gt>

Brian J. Watson E<lt>F<bjbrew@power.net>E<gt>

Matthew Brett E<lt>matthew.brett@mrc-cbu.cam.ac.ukE<gt>

=cut

package Math::Matrix;
use vars qw($VERSION $eps);
use strict;

$VERSION = 0.4;

use overload
       '~'  => 'transpose',
       '+'  => 'add',
       '-'  => 'subtract',
       '*'  => 'multiply',
       '""' => 'as_string';

sub version {
    return "Math::Matrix $VERSION";
}

# Implement - array copy, inheritance 

# class call - new matrix as input
# object call - creates matrix with same dimensions matrix

sub new {
    my $that = shift;
    my $class = ref($that) || $that;
    my $self = [];
    if (ref($that) && (@_ == 0)) { # object call no args -> copy matrix
	for (@$that) {
	    push(@{$self}, [map {0} @{$_}]);
	}
    } else { # class call / object call -> matrix as input
	my $len = scalar(@{$_[0]});
	for (@_) {
	    return undef if scalar(@{$_}) != $len;
	    push(@{$self}, [@{$_}]);
	}
    }
    bless $self, $class;
}

sub clone {
    my $that = shift;
    my $self = [];

    for (@$that) {
        push(@{$self}, [@{$_}]);
    }
    bless $self, ref($that)||$that;
}

sub size {
    my $self = shift;
    my $m = @{$self};
    my $n = @{$self->[0]};
    ($m, $n);
}

sub concat {
    my $self   = shift;
    my $other  = shift;
    my $result =  $self->clone();

    return undef if scalar(@{$self}) != scalar(@{$other});
    for my $i (0 .. $#{$self}) {	
	push @{$result->[$i]}, @{$other->[$i]};
    }
    $result;
}

sub transpose {
    my $self = shift;
    my $class = ref($self);
    my @result;
    my $m;

    for my $col (@{$self->[0]}) {
        push @result, [];
    }
    for my $row (@{$self}) {
        $m=0;
        for my $col (@{$row}) {
            push(@{$result[$m++]}, $col);
        }
    }
    $class->new(@result);
}

sub vekpro {
    my($a, $b) = @_;
    my $result=0;

    for my $i (0 .. $#{$a}) {
        $result += $a->[$i] * $b->[$i];
    }
    $result;
}

sub multiply {
    my $self  = shift;
    my $class = ref($self);
    my $other = shift->transpose;
    my @result;
    my $m;

    return undef if $#{$self->[0]} != $#{$other->[0]};
    for my $row (@{$self}) {
        my $rescol = [];
	for my $col (@{$other}) {
            push(@{$rescol}, vekpro($row,$col));
        }
        push(@result, $rescol);
    }
    $class->new(@result);
}

$eps = 0.00001;

sub solve {
    my $self  = shift;
    my $class = ref($self);

    my $m    = $self->clone();
    my $mr   = $#{$m};
    my $mc   = $#{$m->[0]};
    my $f;
    my $try;

    return undef if $mc <= $mr;
    ROW: for(my $i = 0; $i <= $mr; $i++) {
	$try=$i;
	# make diagonal element nonzero if possible
	while (abs($m->[$i]->[$i]) < $eps) {
	    last ROW if $try++ > $mr;
	    my $row = splice(@{$m},$i,1);
	    push(@{$m}, $row);
	}

	# normalize row
	$f = $m->[$i]->[$i];
	for(my $k = 0; $k <= $mc; $k++) {
            $m->[$i]->[$k] /= $f;
	}
	# subtract multiple of designated row from other rows
        for(my $j = 0; $j <= $mr; $j++) {
	    next if $i == $j;
            $f = $m->[$j]->[$i];
            for(my $k = 0; $k <= $mc; $k++) {
                $m->[$j]->[$k] -= $m->[$i]->[$k] * $f;
            }
        }
    }
# Answer is in augmented column
    transpose $class->new(@{$m->transpose}[$mr+1 .. $mc]);
}

sub print {
    my $self = shift;

    print @_ if scalar(@_);
    print $self->as_string;
}

sub as_string {
    my $self = shift;
    my $out = "";
    for my $row (@{$self}) {
        for my $col (@{$row}) {
            $out = $out . sprintf "%10.5f ", $col;
        }
        $out = $out . sprintf "\n";
    }
    $out;
}

sub new_identity {
  my $type = shift;
  my $class = ref($type) || $type;
  my $self = [];
  my $size = shift;

  for my $i (1..$size) {
    my $row = [];
    for my $j (1..$size) {
      push @$row, $i==$j ? 1 : 0;
    }
    push @$self, $row;
  }
  bless $self, $class;
}

sub eye {
    &new_identity(@_);
}

sub multiply_scalar {
  my $self = shift;
  my $factor = shift;
  my $result = $self->new();

  my $last = $#{$self->[0]};
  for my $i (0 .. $#{$self}) {
    for my $j (0 .. $last) {
      $result->[$i][$j] = $factor * $self->[$i][$j];
    }
  }
  $result;
}

sub negative {
  shift->multiply_scalar(-1);
}

sub subtract {
    my $self = shift;
    my $other = shift;
    $self->add($other->negative);
}

sub equal {
  my $A = shift;
  my $B = shift;
  my $ok = 1;

  my $last = $#{$A->[0]};
  for my $i (0 .. $#{$A}) {
    for my $j (0 .. $last) {
      $A->[$i][$j] == $B->[$i][$j] or $ok=0;
    }
  }
  $ok;
}

sub add {
  my $self = shift;
  my $other = shift;
  my $result = $self->new();

  return undef
    if $#{$self} != $#{$other};

  my $last= $#{$self->[0]};
  return undef
    if $last != $#{$other->[0]};
  for my $i (0 .. $#{$self}) {
    for my $j (0 .. $last) {
      $result->[$i][$j] = $self->[$i][$j] + $other->[$i][$j];
    }
  }
  $result;
}

sub slice {
  my $self = shift;
  my $class = ref($self);
  my $result = $class->new([]);

  foreach my $j (@_) {
    for my $i (0..$#{$self}) {
      push @{$result->[$i]}, $self->[$i][$j];
    }
  }
  $result;
}

sub determinant {
  my $self = shift;
  my $class = ref($self);
  my $last= $#{$self->[0]};

  return undef
    unless $last == $#{$self};

  if ($last == 0) {
    return $self->[0][0];
  } else {
    my $result = 0;
    foreach my $col (0..$last) {
      my $matrix = $self->slice(0..$col-1,$col+1..$last);
      $matrix = $class->new(@$matrix[1..$last]);
      my $term += $matrix->determinant();
      $term *= $self->[0][$col];
      $term *= $col % 2 ? -1 : 1;
      $result += $term;
    }
    return $result;
  }
}

#
# For vectors only
#

sub dot_product {
  my $vector1 = shift;
  my $vector2 = shift;

  $vector1 = $vector1->transpose()
    unless @$vector1 == 1;
  return undef
    unless @$vector1 == 1;

  $vector2 = $vector2->transpose()
    unless @{$vector2->[0]} == 1;
  return undef
    unless @{$vector2->[0]} == 1;

  return $vector1->multiply($vector2)->[0][0];
}

sub absolute {
  my $vector = shift;
  sqrt $vector->dot_product($vector);
}

sub normalize {
  my $vector = shift;
  my $length = $vector->absolute();
  return undef
    unless $length;
  $vector->multiply_scalar(1 / $length);
}

sub cross_product {
  my $vectors = shift;
  my $class = ref($vectors);

  my $dimensions = @{$vectors->[0]};
  return undef
    unless $dimensions == @$vectors + 1;

  my @axis;
  foreach my $column (0..$dimensions-1) {
    my $tmp = $vectors->slice(0..$column-1,
                              $column+1..$dimensions-1);
    my $scalar = $tmp->determinant;
    $scalar *= ($column % 2) ? -1 : 1;
    push @axis, $scalar;
  }
  my $axis = $class->new(\@axis);
  $axis = $axis->multiply_scalar(($dimensions % 2) ? 1 : -1);
}

1;

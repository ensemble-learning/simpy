package Packages::Molecule;

use strict;
require Exporter;
use Carp;

our (@ISA, @EXPORT, $VERSION, @EXPORT_OK);
local $Carp::Carplevel = 1;
my ($cpack, $cfile) = caller();

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
$VERSION = "1.00";

#Initilization of Class
sub create {
    my $invocant = shift;
    my $class = ref($invocant) || $invocant;
    my $info = {
	   id       => -1,
	   atom     => [],
	   atomlist => {},
	   com      => {
			x    => 0,
			y    => 0,
			z    => 0,
			},
	   count    => 0,
	   mass     => 0,
    };
    my $self = sub {
        my $field = lc shift;

        # Access checks
        croak "No valid field '$field' in object"
           unless exists($info->{$field});
        if (@_) {
	    my $val = shift;
	    return $info->{atom}->[$val -1] if ($field eq "atom");
	    $info->{$field} = $val;
	}
        return $info->{$field};
    };

    bless($self, $class);
    return $self;
}

for my $field (qw(id atom com count atomlist mass)) {
    no strict "refs";
    *$field = sub {
        my $self = shift;
        return $self->(uc $field, @_);
     };
}

my $getRotations = sub {
    my $self = shift;
    my $axes = shift;
    my @vals = @_;
    my ($rotations, $rotVec, $i, $j, $rec);

    while ($axes =~ /(x|y|z)/g) {
        $i = $1;
        next if (exists($rotations->{$i}));
        $j = 0;
        while ($j <= $#vals) {
            if ($vals[$j] =~ /(\-?\d+\.?\d*)/) {
                splice @vals, $j, 1;
		$rotVec->{$i} = $1;
                last;
            } else {
                $j++;
	    }
        }
    }
    for $i ("x", "y", "z") {
	next if ($rotVec->{$i});
	$rotVec->{$i} = 0;
    }
    return $rotVec;
};

my $getRotationMatrix = sub {
    my $self = shift;
    my $rotationVector = shift;
    my ($dim, $count, $rotationMatrix, @cosA, @sinA);

    $count = 0;
    foreach $dim ("x", "y", "z") {
	$sinA[$count] = sin($rotationVector->{$dim});
	$cosA[$count] = cos($rotationVector->{$dim});
	$count++;
    }

   $rotationMatrix = {
                         "x" => [$cosA[1]*$cosA[2], $sinA[0]*$sinA[1]*$cosA[2]+$cosA[0]*$sinA[2], -$cosA[0]*$sinA[1]*$cosA[2]+$sinA[0]*$sinA[2]],
                         "y" => [-$cosA[1]*$sinA[2], -$sinA[0]*$sinA[1]*$sinA[2]+$cosA[0]*$cosA[2], $cosA[0]*$sinA[1]*$sinA[2]+$sinA[0]*$cosA[2]],
                         "z" => [$sinA[1], -$sinA[0]*$cosA[1], $cosA[0]*$cosA[1]],

		     };
    return $rotationMatrix;
};

sub print {
    my $self = shift;
    my ($i, $currAtom);

    print "Molecule " . $self->id . " Info\n====================\n";
    for $i qw(id count) {
	printf "%-10s  " . $self->$i . "\n", uc $i;
    }

    for $i (1 ... $self->count) {
	$currAtom = $self->atom->[$i - 1];
	$currAtom->print;
    }
}

sub addAtom {
    my $self = shift;
    my $currAtom = shift;
    my $atomCount = $self->count;
    my $mass = $self->mass;
    my $index;

    if (! exists($self->atomlist->{ $currAtom->code })) {
	push @{ $self->atom }, $currAtom;
	$index = $#{ $self->atom };
	$atomCount++;
	$self->atomlist->{ $currAtom->code } = $index;
	for $index ("x", "y", "z") {
	    $self->com->{$index} += ($currAtom->$index * $currAtom->mass);
	}
	$mass += $currAtom->mass;
    }
    $self->("mass", $mass);
    $self->("count", $atomCount)
}

sub removeAtom {
    my $self = shift;
    my $currAtom = shift;
    my $atomCount = $self->count;
    my ($atomIndex, $index, $mass);

    if ($self->atomlist->{ $currAtom->code }) {
	$atomIndex = $self->atomlist->{ $currAtom->code };
	delete($self->atomlist->{ $currAtom->code });
	splice(@{ $self->atom }, $atomIndex, 1);
	$atomCount--;
    }
    $self->("count", $atomCount);
    for $index ("x", "y", "z") {
        $self->com->{$index} -= ($currAtom->$index * $currAtom->mass);
    }
    $mass = $self->mass - $currAtom->mass;
    $self->("mass", $mass);
}

sub CoM {
    my $self = shift;
    my ($i, $retStr);

    return 0 if ($self->mass == 0);
    for $i ("x", "y", "z") {
	$retStr->{$i} = $self->com->{$i}/$self->mass;
    }

    return $retStr;
}

sub rotate {
    my $self = shift;
    my $axes = lc shift;
    return if (! @_);
    my $angles = join " ", @_;
    my ($rotMatrix, $rotVec, $currPos, $newPos, $atom, $dim);

    return if ($axes !~ /^(x|y|z)/);
    return if ($angles !~ /^\s*(\-?\d+\.?\d*)/);
    $rotVec = $self->$getRotations($axes, @_);
    $rotMatrix = $self->$getRotationMatrix($rotVec);
    $self->com->{x} = $self->com->{y} = $self->com->{z} = 0;
    for $atom (@{ $self->atom }) {
	for $dim ("x", "y", "z") {
	    $currPos = $atom->$dim;
	    $newPos = ($atom->x * $rotMatrix->{$dim}[0] + $atom->y * $rotMatrix->{$dim}[1] + $atom->z * $rotMatrix->{$dim}[2]);
	    $atom->move($dim, $newPos);
	    $self->com->{$dim} += $newPos;
	}
    }
}

sub trans {
    my $self = shift;
    my $dimension = lc shift;
    my $val = shift;
    my $atom;

    return if ($dimension !~ /^(x|y|z)/);
    $dimension = $1;
    $val = "+ ${val}" if ($val !~ /^(\-|\+)\s+/);
    return if ($val !~ /^(\+|\-)\s+(\-?\d+\.?\d*)/);
    for $atom (@{ $self->atom }) {
	$atom->move($dimension, $val);
    }
    for ("x", "y", "z") {
	$self->com->{$_} += $2;
    }
}

1;

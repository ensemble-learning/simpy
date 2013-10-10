package Packages::Atom;

use strict;
use Carp;
use Packages::Bond;
require Exporter;

our (@ISA, @EXPORT, $VERSION, @EXPORT_OK);
local $Carp::Carplevel = 1;
my ($cpack, $cfile) = caller();
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
$VERSION = "1.00";

# Initialize the atom object
sub new {
    my $invocant = shift;
    my $class = ref($invocant) || $invocant;
    my $i;
    my $data = {
	code      => -1,
	index     => -1,
	label     => undef,
	atmname   => undef,
	resname   => undef,
	resid     => undef,
	chain     => undef,
	x         => undef,
	y	  => undef,
	z	  => undef,
	fftype    => undef,
	numbonds  => undef,
	lonepairs => undef,
	charge    => undef,
	bondlist  => Packages::Bond->new(),
	molid     => undef,
	occupancy => undef,
	resonance => undef,
	radii     => undef,
	energy    => 0,
	mass      => 0,
	element   => undef,
	};
    my $self = sub {
	my $field = lc shift;

        # Access checks
	croak "No valid field '$field' in object"
	   unless (exists($data->{$field}) or ($field =~ /print|fields/));
	#croak "Invalid access to Bond class"
	   #unless $field ne "BONDS";
	if ($field eq "print") {
	    print "\nAtom " . $data->{code} . " Info\n===============\n";
	    for $i (keys %{ $data }) {
		printf "%-10s " . $data->{$i} . "\n", uc($i) if (defined($data->{$i}) and $i ne "bondlist" and $data->{$i});
	    }
	    $data->{bondlist}->print;
	    return ();
	} elsif ($field eq "fields") {
	    my $fieldlist = "";
	    for (keys %{ $data }) {
		$fieldlist .= "$_ ";
	    }
	    return $fieldlist;
	}
	if (@_) { $data->{$field} = shift }
	return $data->{$field};
    };
    bless($self, $class);
    return $self;
}

for my $field (qw(code index label atmname resname resid chain x y z fftype numbonds lonepairs charge bondlist molid occupancy resonance radii energy print fields mass element)) {
    no strict "refs";
    *$field = sub {
	my $self = shift;
	return $self->(lc $field, @_);
     };
}

sub move {
    my $self = shift;
    my $field = shift;
    my $val = shift;
    my $curr = $self->$field;

    $self->($field, $val) if ($val =~ /^\-?\d+\.?\d*/);
    $self->($field, ($curr - $1)) if ($val =~ /^\-\s+(\-?\d+\.?\d*)/);
    $self->($field, ($curr + $1)) if ($val =~ /^\+\s+(\-?\d+\.?\d*)/);
}

sub bonds {
    my $self = shift;
    my $val = shift;
    my $currBond = $self->bondlist->index($val);

    return $currBond->data;
}

1;

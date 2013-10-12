package Packages::ListItem;

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
sub spawn {
    my $invocant = shift;
    my $class = ref($invocant) || $invocant;
    my $info = {
	   data  => {},
	   next  => undef,
	   prev  => undef,
    };
    my $self = sub {
        my $field = lc shift;

        # Access checks
        croak "No valid field '$field' in object"
           unless exists($info->{$field});
        if (@_) { $info->{$field} = shift }
        return $info->{$field};
    };

    bless($self, $class);
    return $self;
}

for my $field (qw(data next prev)) {
    no strict "refs";
    *$field = sub {
        my $self = shift;
        return $self->(uc $field, @_);
     };
}

sub DESTROY {
    my $self = shift;
    $self->("data", undef);
    $self->("next", undef);
    $self->("prev", undef);
}

1;

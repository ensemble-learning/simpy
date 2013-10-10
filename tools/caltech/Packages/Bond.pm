package Packages::Bond;
# Creates a linked list with additional fields for bonds
use strict;
use Packages::BondItem;
use Carp;
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
    my $data = {
	    count => 0,
	    first => Packages::BondItem->spawn(),
	    last  => Packages::BondItem->spawn(),
	};
    my $self = sub {
	my $field = lc shift;

        # Access checks
	croak "No valid field '$field' in object"
	   unless exists($data->{$field});
	if (@_) { $data->{$field} = shift }
	return $data->{$field};
    };
    bless($self, $class);
    return $self;
}

for my $field (qw(count first last)) {
    no strict "refs";
    *$field = sub {
	my $self = shift;
	return $self->(uc $field, @_);
     };
}

sub clone {
    my $self = shift;
    my $newBondlist = shift;
    my ($bondAtom, $count, $i, $prev, $curr);

    $count = $self->count;
    return () if (! $count);
    $newBondlist->first->("data", $self->first->data);
    $newBondlist->last->("data", $self->last->data);
    $newBondlist->first->("next", Packages::BondItem->spawn());

    $curr = $newBondlist->first->next;
    $prev = $newBondlist->first;
    
    $bondAtom = $self->first;
    for $i (2 .. $count) {
	$curr->("prev", $prev);
	$prev->("next", $curr);
	$bondAtom = $bondAtom->next;
	$curr->("data", $bondAtom->data);
	$curr->("next", Packages::BondItem->spawn());
	$prev = $curr;
	$curr = $curr->next;
    }
    $newBondlist->("last", $prev);
    $newBondlist->first->("prev", $newBondlist->last);
    $newBondlist->last->("next", $newBondlist->first);
}

sub print {
    my $self = shift;
    my ($atom);
    return "" if ($self->count == 0);
    
    print "Bonds\n--------\n";
    $atom = $self->first;
    print "Atom " . $atom->data->index;
    $atom = $atom->next;
    while ($atom ne $self->first) {
	print " Atom " . $atom->data->index;
	$atom = $atom->next;
    }
    print "\n";
}

sub store {
    my $self = shift;
    my $partner = shift;
    my $field = shift;
    my $val = shift;
    my $curr = $self->findParent($partner);

    $curr->("$field", "$val");
}

sub create {
    my $self = shift;
    my $partner = shift;
    my ($currCount, $rec, $prev, $next);
    croak "Attempted to create link to nothing!\n"
	unless (defined($partner));
    $currCount = $self->count;
    $rec = Packages::BondItem->spawn();
    $rec->("data", $partner);

    if ($currCount == 0) { 
    # this is the first entry in to the linked Bond so update first and last pointers
	$self->("first", $rec);
	$self->("last", $rec);
	$prev = $self->first;
	$next = $self->last;
	$self->first->("next", $next);
	$self->first->("prev", $next);
	$self->last->("next", $prev);
	$self->last->("prev", $prev);
    } else { # create a new item as the last item
	return () if ($self->find($partner)); # do not create link to already linked items
	$prev = $self->last;
	$next = $self->first;
	$rec->("prev", $prev);
	$rec->("next", $next);
	$self->first->("prev", $rec);
	$self->last->("next", $rec);
	$self->("last", $rec);
    }
    $currCount++;
    $self->("count", $currCount);
}

sub remove {
    my $self = shift;
    my $partner = shift;
    my ($currCount, $currlink, $islink, $prev, $next);

    $islink = 0;
    if ($self->first->data->code eq $partner->code) { # removing the head of the Bond
	if ($self->count == 1) { # removing the only link so set everything to null
	    $self->("first", undef);
	    $self->("last", undef);
	} else {
	    $currlink = $self->first->next;
	    $currlink->("prev", $self->last);
	    $currlink->("next", $self->first->next->next);
	    $self->last->("next", $currlink);
	    $self->first->("next", undef);
	    $self->first->("prev", undef);
	    $self->("first", undef);
	    $self->("first", $currlink);
	}
	$islink = 1;
    } elsif ($self->last->data->code eq $partner->code) { # removing the tail of the Bond
	$currlink = $self->last->prev;
	$currlink->("next", $self->first);
	$currlink->("prev", $self->last->prev->prev);
	$self->first->("prev", $currlink);
	$self->last->("next", undef);
	$self->last->("prev", undef);
	$self->("last", undef);
	$self->("last", $currlink);
	$islink = 1;
    } else { # search to see if partner is in the link Bond of curr atom
        $currlink = $self->find($partner);
        return () if (! $currlink);
        $prev = $currlink->prev;
        $next = $currlink->next;
        $prev->("next", $next);
        $next->("prev", $prev);
        $currlink->("next", undef);
        $currlink->("prev", undef);
        $currlink->DESTROY();
	$currlink = undef;
        $islink = 1;
    }

    return () if (! $islink);
    $currCount = $self->count;
    $currCount--;
    $self->("count", $currCount);
}

sub index {
    my $self = shift;
    my $sIndex = shift;
    my ($currItem, $count, $totItems);

    $totItems = $self->count;
    return $self->first if ($sIndex <= 1);
    return $self->last if ($sIndex == $self->count);
    return () if ($sIndex > $totItems);
    $count = 2;
    $currItem = $self->first->next;
    while ($currItem ne $self->first) {
	return $currItem if ($count == $sIndex);
	$currItem = $currItem->next;
    }
}

sub findIndex {
    my $self = shift;
    my $aIndex = shift;
    my ($currItem);

    return $self->first if ($self->first->data->index == $aIndex);

    $currItem = $self->first->next;
    while ($currItem ne $self->first) {
        return $currItem if ($currItem->data->index == $aIndex);
        $currItem = $currItem->next;
    }
}

sub find {
    my $self = shift;
    my $partner = shift;
    my ($currItem, $index);

    $index = $partner->index;

    return $self->first if ($self->first->data->index == $index);

    $currItem = $self->first->next;
    while ($currItem ne $self->first) {
        return $currItem if ($currItem->data->index == $index);
        $currItem = $currItem->next;
    }
}

sub findParent {
    my $self = shift;
    my $partner = shift;
    my ($currItem, $index);

    $index = $partner->data->index;

    return $self->first if ($self->first->data->index == $index);

    $currItem = $self->first->next;
    while ($currItem ne $self->first) {
        return $currItem if ($currItem->data->index == $index);
        $currItem = $currItem->next;
    }
}

1;

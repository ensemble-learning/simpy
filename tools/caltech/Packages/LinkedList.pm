package Packages::LinkedList;

use strict;
use Packages::ListItem;
use Carp;
require Exporter;

our (@ISA, @EXPORT, $VERSION, @EXPORT_OK);
local $Carp::Carplevel = 1;
my ($cpack, $cfile) = caller();
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
$VERSION = "1.00";

# Initialize the link list object
sub new {
    my $invocant = shift;
    my $class = ref($invocant) || $invocant;
    my $data = {
	    count => 0,
	    first => Packages::ListItem->spawn(),
	    last  => Packages::ListItem->spawn(),
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

sub insert {
    my $self = shift;
    my $partner = shift;
    my $insertPt = shift;
    my ($currCount, $rec, $prev, $next);
    croak "Attempted to create link to nothing!\n"
	unless (defined($partner));
    $currCount = $self->count;
    $rec = Packages::ListItem->spawn();
    $rec->("data", $partner);

    if ($currCount == 0) { 
    # this is the first entry in to the linked list so update first and last pointers
	$self->("first", $rec);
	$self->("last", $rec);
	$prev = $self->first;
	$next = $self->last;
	$self->first->("next", $next);
	$self->first->("prev", $next);
	$self->last->("next", $prev);
	$self->last->("prev", $prev);
    } elsif (! defined($insertPt)) { # create a new item as the last item
	#return () if ($self->find($partner)); # do not create link to already linked items
	$prev = $self->last;
	$next = $self->first;
	$rec->("prev", $prev);
	$rec->("next", $next);
	$self->first->("prev", $rec);
	$self->last->("next", $rec);
	$self->("last", $rec);
    } else { # insert into the middle of the list
        $prev = $insertPt->prev;
        $next = $insertPt;
        $rec->("prev", $prev);
        $rec->("next", $next);
	$prev->("next", $rec);
	$next->("prev", $rec);
    }
    $currCount++;
    $self->("count", $currCount);
}

sub delete {
    my $self = shift;
    my $partner = shift;
    my ($currCount, $currlink, $islink, $prev, $next);

    $islink = 0;
    if ($self->first->data->code eq $partner->code) { # removing the head of the list
	if ($self->count == 1) { # removing the only link so set everything to null
	    $self->("first", undef);
	    $self->("last", undef);
	} else {
	    $currlink = $self->first->next;
	    $currlink->("prev", $self->first->prev);
	    $self->last->("next", $currlink);
	    $self->first->("next", undef);
	    $self->first->("prev", undef);
	    $self->("first", undef);
	    $self->("first", $currlink);
	}
	$islink = 1;
    } elsif ($self->last->data->code eq $partner->code) { # removing the tail of the list
	$currlink = $self->last->prev;
	$currlink->("next", $self->last->next);
	$self->first->("prev", $currlink);
	$self->last->("next", undef);
	$self->last->("prev", undef);
	$self->("last", undef);
	$self->("last", $currlink);
	$islink = 1;
    } else { # search to find curritem in list
	$currlink = $self->find($partner);
	return () if (! $currlink);
	$prev = $currlink->prev;
	$next = $currlink->next;
	$prev->("next", $next);
	$next->("prev", $prev);
	$currlink->("next", undef);
	$currlink->("prev", undef);
	$currlink->DESTROY();
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
    my ($currItem);

    return $self->first if ($self->first->data->index == $partner->index);

    $currItem = $self->first->next;
    while ($currItem ne $self->first) {
        return $currItem if ($currItem->data->index == $partner->index);
        $currItem = $currItem->next;
    }
}

1;

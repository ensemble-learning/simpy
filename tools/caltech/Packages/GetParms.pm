package Packages::GetParms;

use strict;
    
our (@ISA, @EXPORT, $VERSION, @EXPORT_OK);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(IsValidCrossovers GetHashKeys);
$VERSION = "1.0";

# Private Methods
my $hKeys = sub {
#    This sub will decend into the self hash and extract all the keys
    my $self = shift;
    my ($returnval) = GetHashKeys($self);
    $returnval =~ s/\$self\-\>//g;
    return $returnval;
};

my $GetParms = sub {
    my $self = shift;
    my $returnval = 1;
    my $parmfile = $self->{"Files"}->{"parameter"};
    my ($in_data, $c_list);

    $self->{isLayout} = 0;
    $self->{CreatePDB} = 1;
    if (open PARMFILE, $parmfile) {
	while (<PARMFILE>) {
	    chomp;
	    $in_data = $_;
	    if ($in_data =~ /Mol:\s+(\d+).(\d+)/) {
		$self->{"Molecule"}->{"name"} = $1 . $2;
		$self->{"Molecule"}->{"periodicity"} = $1 + $2;
		$self->{"Molecule"}->{"major_groove"} = $1;
		$self->{"Molecule"}->{"minor_groove"} = $2;
	    } elsif ($in_data =~ /Total bases:\s+(\d+)/) {
		$self->{"Molecule"}->{"total_bases"} = $1;
	    } elsif ($in_data =~ /Bases at end:\s+(\d+)/) {
		$self->{"Molecule"}->{"bases_at_end"} = $1;
	    } elsif ($in_data =~ /Crossovers:\s+(.+)/) {
		$c_list = $1;
		@{ $self->{"Molecule"}->{"crossovers"} } = split(/ /, $c_list);
	    } elsif ($in_data =~ /Topology file:\s+(.+)\s*$/) {
		$self->{"Files"}->{"topology"} = $1;
	    } elsif ($in_data =~ /Trajectory file:\s+(.+)\s*$/) {
		$self->{"Files"}->{"trajectory"} = $1;
	    } elsif ($in_data =~ /3 prime in:\s+([1|0])/) {
		$self->{"Molecule"}->{"is3PrimeIn"} = $1;
	    } elsif ($in_data =~ /Cuttoff:\s+(\d+\.\d+)/) {
		$self->{"Amber_Options"}->{"cutoff"} = $1;
	    } elsif ($in_data =~ /Run Ptraj:\s+([1|0])/) {
		$self->{"Amber_Options"}->{"run_ptraj"} = $1;
	    } elsif ($in_data =~ /Run Anal:\s+([1|0])/) {
		$self->{"Amber_Options"}->{"run_anal"} = $1;
	    } elsif ($in_data =~ /Host:\s+(\w+)$/) {
		$self->{"Movie"}->{"host"} = $1;
	    } elsif ($in_data =~ /Movie Dimensions:\s+(\d+)x(\d+)/) {
		$self->{"Movie"}->{"size"} = $1 . "," . $2;
	    } elsif ($in_data =~ /Movie Length:\s+(\d+)/) {
		$self->{"Movie"}->{"length"} = $1;
	    } elsif ($in_data =~/Set Layout: ([1|0])/) {
		$self->{isLayout} = $1;
	    } elsif ($in_data =~ /Helix (\d+) Strand (\d+): (.+)/) {
		$self->{isLayout} = 1;
		$self->{Layout}{$1}{$2} = $3;
	    } elsif ($in_data =~ /Create PDB: ([1|0])/) {
		$self->{CreatePDB} = $1;
	    } elsif ($in_data =~ /Reference: (\w+)/) {
		if (-e $1) {
		    $self->{has_reference} = 1;
		    $self->{reference} = $1;
		}
	    } elsif ($in_data =~ /Trajectory type: (\w+)/) {
		$self->{Files}{trjType} = 2 if (lc($1) eq "lammps");
	    } 
	}
	close PARMFILE;
    } else {
	$returnval = 1;
    }
    return $returnval;
};

my $ValidateParms = sub {
    my $self = shift;

    my $returnval = 1;

    if ($#{ $self->{Molecule}{crossovers} } == -1 ||
	$self->{Molecule}{bases_at_end} == 0 ||
	$self->{Molecule}{total_bases} == 0 ||
	$self->{Molecule}{name} eq "" ||
	$self->{Molecule}{is3PrimeIn} == -1 ||
	$self->{Molecule}{periodicity} == 0) {
	print "PARAMS: #crossovers: " . $#{ $self->{"Molecule"}->{"crossovers"} } . "\n" .
	"BASES AT END: " .  $self->{"Molecule"}->{"bases_at_end"} . "\n" .
	"TOTAL BASES: " . $self->{"Molecule"}->{"total_bases"} . "\n" .
	"MOLECULE NAME: " . $self->{"Molecule"}->{"name"} . "\n" .
	"3 Prime End In: " . $self->{"Molecule"}->{"is3PrimeIn"} . "\n" .
	"PERIODICITY: " . $self->{"Molecule"}->{"periodicity"} . "\n";
	$returnval = 0;
#    } elsif ( ($self->{"Molecule"}->{"total_bases"} % 4) > 0 ) {
#	print "Invalid total number of bases: Must be a multiple of 4\n";
#	$returnval = 0;
    } elsif (! IsValidCrossovers($self->{"Molecule"}->{"total_bases"}, 
				 $self->{"Molecule"}->{"crossovers"} ) ) {
	print "Invalid crossover specification. ", 
	"Number must be less than strand length\n";
	$returnval = 0;
    } elsif ((! $self->{"Files"}->{"topology"}) or (! -e $self->{"Files"}->{"topology"})) {
	print "Cannot locate topology file " . $self->{"Files"}->{"topology"} . "\n";
	$returnval = 0;
    } elsif ((! $self->{"Files"}->{"trajectory"}) or (! -e $self->{"Files"}->{"trajectory"})) {
	print "Cannot locate trajectory file " . $self->{"Files"}->{"trajectory"} . "\n";
	$returnval = 0;
    }

    if ($self->{"Movie"}->{"host"} &&
	$self->{"Movie"}->{"size"} && $self->{"Movie"}->{"length"}) {
	$self->{"Movie"}->{"play"} = 1;
	$self->{"Movie"}->{"user_name"} = $ENV{USER};
    }

    return $returnval;
};

# Public Methods

sub IsValidParams {
    my $self = shift;
    $self->{"Files"}->{"parameter"} = $_[0];

    my $returnval;

    if (! -e $self->{"Files"}->{"parameter"}) {
	print "Cannot locate " . $self->{"Files"}->{"parameter"} . ": $!\n";
	$returnval = 0;
    }

    if ($self->$GetParms) {
	if ($self->$ValidateParms) {
	    return 1;
	} else {
	    return 0;
	}
    }else {
	return 0;
    }
}

sub IsValidCrossovers {
    my ($tot_bases, $crossovers) = @_;
    my ($returnval) = 1;

    if ($#{$crossovers} > -1) {
	for (@$crossovers) {
	    if ($_ > ($tot_bases /4)) {
		$returnval = 0;
		last;
	    }
	}
    }else {
	print "ERROR: No crossover specified!\n";
    }

    return $returnval;
}

sub GetHashKeys {
    my $self = shift;
    my ($curr_hash, $curr_path) = @_;

    my ($returnval, $h_keys, $h_vals);

    if (! $curr_hash) {
	$curr_hash = $self;
	$curr_path = "";
    }

    while ( ($h_keys, $h_vals) = each %{ $curr_hash } ) {
	while (keys %{ $curr_hash->$h_keys }) {
	    $curr_path .= "$h_keys->";
	    $returnval .= GetHashKeys(\$curr_hash->$h_keys, $curr_path);
	}
	$curr_path =~ s/\-\>$//g;
	return "$curr_path ";
    }
}

sub init {
    my $self = shift;
    my $data = {
	Files => {
	    topology => "",
	    trajectory => "",
	    trjType => 1,
	},
	Molecule => {
	    name => "",
	    total_bases => 0,
	    bases_at_end => 0,
	    is3PrimeIn => -1,
	    major_groove => 0,
	    minor_groove => 0,
	    periodicity => 0,
	    crossovers => [],
	},
	Movie => {
	    play => 0,
	    user_name => "",
	    host_name => "",
	    size => 0,
	    length => 0,
	},
	Amber_Options => {
	    cutoff => -1,
	    run_anal => 1,
	    run_ptraj => 1,
	},
    };
    return $data;
};

sub new {
    my $invocant = shift;
    my $class = ref($invocant) || $invocant;
    my $self = init();
    $self = bless($self, $class);
    return $self;
}

for my $field (GetHashKeys()) {
    my $slot = __PACKAGE__ . "::$field";
    no strict "refs";
    *$field = sub {
	my $self = shift;
	$self->{$slot} = shift if @_;
	return $self->{$slot};
    };
}

1;

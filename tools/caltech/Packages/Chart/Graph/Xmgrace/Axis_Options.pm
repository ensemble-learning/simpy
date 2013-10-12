#
# Chart::Graph::Xmgrace::Axis_Options.pm contains all the fundamental options
# for an axis.
#
# Options included:
#                  1) Chart::Graph::Xmgrace::Bar_Options
#                  2) Chart::Graph::Xmgrace::Label_Options
#                  3) Chart::Graph::Xmgrace::Tickmark
#                  4) Chart::Graph::Xmgrace::Tick_Options
#                  5) Chart::Graph::Xmgrace::Ticklabel_Options
#                 
## This software product is developed by Esmond Lee and David Moore,
## and copyrighted(C) 1998 by the University of California, San Diego
## (UCSD), with all rights reserved. UCSD administers the CAIDA grant,
## NCR-9711092, under which part of this code was developed.
##
## There is no charge for this software. You can redistribute it and/or
## modify it under the terms of the GNU General Public License, v. 2 dated
## June 1991 which is incorporated by reference herein. This software is
## distributed WITHOUT ANY WARRANTY, IMPLIED OR EXPRESS, OF MERCHANTABILITY
## OR FITNESS FOR A PARTICULAR PURPOSE or that the use of it will not
## infringe on any third party's intellectual property rights.
##
## You should have received a copy of the GNU GPL along with this program.
##
##
## IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY
## PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL
## DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS
## SOFTWARE, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF
## THE POSSIBILITY OF SUCH DAMAGE.
##
## THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE
## UNIVERSITY OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE,
## SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS. THE UNIVERSITY
## OF CALIFORNIA MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES
## OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED
## TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
## PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE
## ANY PATENT, TRADEMARK OR OTHER RIGHTS.
##
##
## Contact: graph-dev@caida.org
##
##

$VERSION = 3.2;

use Chart::Graph::Xmgrace::Axis_Option;		
use Carp;

package  Chart::Graph::Xmgrace::Bar_Options;
@ISA = qw(Chart::Graph::Xmgrace::Axis_Option);		# derived from Axis_Option

sub _init {
    my $self = shift;
    $self->{name} = "bar";
    $self->{print_order} = ["status","color","linestyle","linewidth"];
    $self->{length} = 4;
    $self->{options} = {
			"status" => "on",
			"color" => "1",
			"linestyle" => "1",
			"linewidth" => "1.0",
		       };
}

package Chart::Graph::Xmgrace::Label_Options;
@ISA = qw(Chart::Graph::Xmgrace::Axis_Option);

sub _init {
    my $self = shift;
    $self->{name} = "label";
    $self->{print_order} = ["label","layout","place","char size",
			    "font", "color", "place"];
    $self->{length} = 4;
    $self->{options} = {
			"label" => "",
			"layout" => "para",
			"place" => "auto",
			"char size" => "1.000000",
			"font" => "0",
			"color" => "1",
			"place" => "auto",
		       };
}

sub char_size ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"char size"} = $val;
    return 1;
}

sub print ($$$ ) {
    my $self = shift;
    my ($handle, $name) = @_;
    my $string = "";
    my $substr = "";
  
    foreach $option (@{ $self->{"print_order"} }) {
	my $option_ref = $self->{"options"};
	    
	if ($option eq "label") {
	    $string = "$name  $option \"$option_ref->{$option}\"\n";
	} else {
#	    if ($option eq "place") {
#		if (ref($option_ref->{$option} eq ARRAY) {
#		    print $handle "$name  $self->{name} $option spec\n"; 
#		    $substr = join (", ", (@{ $option_ref->{$option} })); 
#		    $string = "$name  $self->{name} $option " 
#		}
#	    }
	    
	    # print function handles both scalars and lists
	    if (ref($option_ref->{$option}) eq ARRAY) {
		if ($option eq "place") {
		    $self->_printline($handle,"$name  $self->{name} $option spec\n", $self->{"length"}); 
		}

		$substr = join (", ", (@{ $option_ref->{$option} })); 
		$string = "$name  $self->{name} $option $substr\n";
	    } else {
		$string = "$name  $self->{name} $option $option_ref->{$option}\n";
	    }
	}
	$self->_printline($handle, $string, $self->{"length"}); 
    }
}

package Chart::Graph::Xmgrace::Tickmark;
@ISA = qw(Chart::Graph::Xmgrace::Axis_Option);

sub _init {
    my $self = shift;
    my $type_of_mark = shift;	# major or minor
    $self->{name} = "$type_of_mark";
    $self->{print_order} = ["size","color","linewidth","linestyle",
			    "grid"];
    $self->{length} = 4;
    $self->{options} = {
			"size" => "1",
			"color" => "1",
			"linewidth" => "1.0",
			"linestyle" => "1.000000",
			"grid" => "off",
		       };
}

sub print ($$$$ ) {
    my $self = shift;
    my ($handle, $name, $subname) = @_;
    my $string = "";
    my $substr = "";
  
    foreach $option (@{ $self->{"print_order"} }) {
	my $option_ref = $self->{"options"};
	
	# print function handles both scalars and lists
	if (!ref($option_ref->{$option})) {
	    $string = "$name  $subname $self->{name} $option $option_ref->{$option}\n";
	    $self->_printline($handle, $string, $self->{"length"});       
	} elsif (ref($option_ref->{$option}) eq ARRAY) {
	    $substr = join (", ", (@{ $option_ref->{$option} })); 
	    $string = "$name  $subname $self->{name} $option $substr\n";
	    $self->_printline($handle, $string, $self->{"length"}); 
	} else { # a blessed object, uses Axis Option printline
	    print "shouldn't be in here\n";
	}
    }
}

package Chart::Graph::Xmgrace::Tick_Options;
@ISA = qw(Chart::Graph::Xmgrace::Axis_Option);

sub _init {
    my $self = shift;
    $self->{name} = "tick";
    $self->{print_order} = ["status","default","place rounded","major",
			    "minor","in_out_status","place","type"];
    $self->{length} = 4;
    $self->{options} = {
			"status" => "on",
			"default" => "6",
			"place rounded" => "true",
			"major" => new Chart::Graph::Xmgrace::Tickmark("major"),
			"minor" => new Chart::Graph::Xmgrace::Tickmark("minor"),
			"in_out_status" => "in",
			"place" => "both",
			"type" => "auto",
		       };
}

sub place_rounded ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"place rounded"} = $val;
    return 1;
}

sub print ($$$$ ) {
    my $self = shift;
    my ($handle, $name) = @_;
    my $string = "";
    my $substr = "";
  
    foreach $option (@{ $self->{"print_order"} }) {
	my $option_ref = $self->{"options"};	
	# special case for "status" and "in_out_status"
	if ($option eq "status" or $option eq "in_out_status") {   
	    $string = "$name  $self->{name} $option_ref->{$option}\n";
	    $self->_printline($handle, $string, $self->{"length"});       
	} else {

	# print function handles both scalars and lists
	    if (!ref($option_ref->{$option})) {
		$string = "$name  $self->{name} $option $option_ref->{$option}\n";
		$self->_printline($handle, $string, $self->{"length"});       
	    } elsif (ref($option_ref->{$option}) eq ARRAY) {
		$substr = join (", ", (@{ $option_ref->{$option} })); 
		$string = "$name  $self->{name} $option $substr\n";
		$self->_printline($handle, $string, $self->{"length"}); 
	    } else { # a blessed object, uses Tickmark print
		$option_ref->{$option}->print($handle, $name, $self->{name});
	    }
	}
    }
}



package Chart::Graph::Xmgrace::Ticklabel_Options;
@ISA = qw(Chart::Graph::Xmgrace::Axis_Option);

sub _init {
    my $self = shift;
    $self->{name} = "ticklabel";
    $self->{ticklabels} = undef; # user can put in a ref of array of arrays of ticklabels
    # i.e. [ ["small\\nfoo", 10], ["medium\\nfoo", 20], ["large\\nfoo", 30] ]
    $self->{print_order} = ["status", "prec", "format", "append", "prepend", 
			    "angle", "skip", "stagger", "place", "offset",
			    "sign", "start type", "start", "stop type", "stop",
			    "char size", "font", "color", "type"];
    $self->{length} = 4;
    $self->{options} = {
			"status" => "on",
			"prec" => "5",
			"format" => "general",
			"append" => "",
			"prepend" => "",
			"angle" => "0",
			"skip" => "0",
			"stagger" => "0",
			"place" => "normal",
			"offset" => ["0.000000", "0.000000"],
			"sign" => "normal",
			"start type" => "auto",
			"start" => "0.000000",
			"stop type" => "auto",
			"stop" => "0.000000",
			"char size" => "1.000000",
			"font" => "0",
			"color" => "1",
			"type" => "auto",
		       };
}

sub ticklabels ($$) {
    $self = shift;
    $val = shift;
    $self->{ticklabels} = $val;
    return 1;
}

sub char_size ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"char size"} = $val;
    return 1;
}

sub start_type ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"start type"} = $val;
    return 1;
}

sub stop_type ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"stop type"} = $val;
    return 1;
}

sub print ($$$ ) {
    my $self = shift;
    my ($handle, $name) = @_;
    my $string = "";
    my $substr = "";
  
    foreach $option (@{ $self->{"print_order"} }) {
	my $option_ref = $self->{"options"};

	if ($option eq "status" or $option eq "in_out_status") {   
	  
	    # we first check the status of the option, whether it's on/off
	    # if it's off, we don't print it out
	    if ($option_ref->{"status"} eq "off") {
		$string = "$name  $self->{name} $option_ref->{$option}\n";
		$self->_printline($handle, $string, $self->{"length"});
		last;
	    }
	    $string = "$name  $self->{name} $option_ref->{$option}\n";
	    $self->_printline($handle, $string, $self->{"length"});

	} elsif ($option eq "type" and $self->{ticklabels}) {

	    $string = "$name  $self->{name} $option $option_ref->{$option}\n";
	    $self->_printline($handle, $string, $self->{"length"});
	    my @ticks = @{$self->{"ticklabels"}};
	    my $index;
	    my $count = $#ticks+1;
	    $string = "$name  tick spec $count\n";
	    $self->_printline($handle, $string, $self->{"length"});
	    # print out the tick labels
	    foreach $index (0..$#ticks) {
		my ($label, $value) = @{$ticks[$index]};
		my $ticklabel = "$name  tick major $index,\t\t$value\n";
		$self->_printline($handle, $ticklabel, $self->{"length"});
		$ticklabel = "$name  $self->{name} $index, \"$label\"\n";
		$self->_printline($handle, $ticklabel, $self->{"length"});
	    }
	} else {
	    if ($option eq "append" or $option eq "prepend") {
		$string = "$name  $self->{name} $option \"$option_ref->{$option}\"\n";
	    } else {
		# print function handles both scalars and lists
		if (ref($option_ref->{$option}) eq ARRAY) {
		    $substr = join (", ", (@{ $option_ref->{$option} })); 
		    $string = "$name  $self->{name} $option $substr\n";
		} else {
		    $string = "$name  $self->{name} $option $option_ref->{$option}\n";
		}
	    }
	    $self->_printline($handle, $string, $self->{"length"}); 
	}
    }
}

1;



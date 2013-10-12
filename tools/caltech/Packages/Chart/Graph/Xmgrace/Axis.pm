#
# Axis.pm is a class that contains the fundamental options for 
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
package Chart::Graph::Xmgrace::Axis;
@ISA = qw(Chart::Graph::Xmgrace::Base_Option);
use Carp;
# for Bar_Options, Label_Options, Tick_Options, Ticklabel_Options
use Chart::Graph::Xmgrace::Axis_Options;

sub _init {
    my $self = shift;
    my $type_of_axis = shift;	#xaxis, yaxis, altxaxis, altyaxis  
    $self->{name} = $type_of_axis;
    $self->{print_order} = ["status", "type zero", "offset", "bar options", 
			    "label options", "tick options","ticklabel options"];
    $self->{length} = 4;
    $self->{options} = {
			"status" => "on",
			"type zero" => "false",
			"offset" => ["0.000000","0.000000"],
			"bar options" => new Chart::Graph::Xmgrace::Bar_Options,
			"label options" => new Chart::Graph::Xmgrace::Label_Options,
			"tick options" => new Chart::Graph::Xmgrace::Tick_Options,
			"ticklabel options" => new Chart::Graph::Xmgrace::Ticklabel_Options,
		       };
}

sub type_zero ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"type zero"} = $val;
    return 1;
}

sub bar_options ($) {
    $self = shift;
    return $self->{options}->{"bar options"};
}

sub label_options ($) {
    $self = shift;
    return $self->{options}->{"label options"};
}

sub tick_options ($) {
    $self = shift;
    return $self->{options}->{"tick options"};
}

sub ticklabel_options ($) {
    $self = shift;
    return $self->{options}->{"ticklabel options"};
}

sub print($$ ) {
    my $self = shift;
    my $handle = shift;
    my $string = "";
    my $substr = "";

    foreach $option (@{ $self->{"print_order"} }) {
	my $option_ref = $self->{"options"};

	if ($option eq "status" or $option eq "in_out_status") {   
	    
	    # we first check the status of the axes, whether it's on/off
	    # if it's off, we don't print it out
	    if ($option_ref->{"status"} eq "off") {
		$string = "$self->{name}  $option_ref->{$option}\n";
		$self->_printline($handle, $string, $self->{"length"});
		last;
	    }

	    $string = "$self->{name}  $option_ref->{$option}\n";
	    $self->_printline($handle, $string, $self->{"length"});       
	    
	} else {
	    # print function handles both scalars and lists
	    if (!ref($option_ref->{$option})) {
		$string = "$self->{name}  $option $option_ref->{$option}\n";
		$self->_printline($handle, $string, $self->{"length"});       
	    } elsif (ref($option_ref->{$option}) eq ARRAY) {
		$substr = join (", ", (@{ $option_ref->{$option} })); 
		$string = "$self->{name}  $option $substr\n";
		$self->_printline($handle, $string, $self->{"length"}); 
	    } else { # a blessed object, uses Axis Option print
		$option_ref->{$option}->print($handle, $self->{name});
	    } 
	}
    }
}

1;

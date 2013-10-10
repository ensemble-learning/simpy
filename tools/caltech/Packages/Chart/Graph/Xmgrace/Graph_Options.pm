#
# Graph_Options.pm contains all the different options (objects) for a particular graph
#
# options included:
#                  1) Chart::Graph::Xmgrace::Graph_Global_Options
#                  2) Chart::Graph::Xmgrace::Global_Options
#                  3) Chart::Graph::Xmgrace::Stack_Options
#                  4) Chart::Graph::Xmgrace::View_Options
#                  5) Chart::Graph::Xmgrace::Title_Options
#                  6) Chart::Graph::Xmgrace::Subtitle_Options
#                  7) Chart::Graph::Xmgrace::Legend_Options
#                  8) Chart::Graph::Xmgrace::Frame_Options
#                  9) Chart::Graph::Xmgrace::Extra_Options
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

use Chart::Graph::Xmgrace::Base_Option;
use Carp;
package Chart::Graph::Xmgrace::Graph_Global_Options;
@ISA = qw(Chart::Graph::Xmgrace::Base_Option);

sub _init {
    my $self = shift;
    my $graph_number = shift;

    $self->{name} = $graph_number;
    $self->{print_order} = ["stacked","hidden","type","bar hgap","fixedpoint",
			    "fixedpoint type","fixedpoint xy","fixedpoint format",
			    "fixedpoint prec","with"];
    $self->{length} = 0;
    $self->{options} = {
			"stacked" => "false",
			"hidden" => "false",
			"type" => "XY",
			"bar hgap" => '0.000000',
			"fixedpoint" => "off",
			"fixedpoint type" => "0",
			"fixedpoint xy" => ["0.000000","0.000000"],
			"fixedpoint format" => "general general",
			"fixedpoint prec" => ["6","6"],
			"with" => $graph_number,
		       };
}

sub print($$ ) {
    my $self = shift;
    my $handle = shift;
    my $string = "";
    my $substr = "";		# for making arrays into a string

    foreach $option (@{ $self->{"print_order"} }) {
	my $option_ref = $self->{"options"};

	# print function handles both scalars and lists
	if ($self->{name}) {
	    if (ref($option_ref->{$option}) eq ARRAY) {
		$substr = join (", ", (@{ $option_ref->{$option} })); 
		$string = "$self->{name} $option $substr\n";
	    } elsif ($option eq "with") {
		$string = "$option $option_ref->{$option}\n";
	    } else {
		$string = "$self->{name} $option $option_ref->{$option}\n";
	    }
	} else {		# global options don't have a name field
	    if (ref($option_ref->{$option}) eq ARRAY) {
		$substr = join (", ", (@{ $option_ref->{$option} })); 
		$string = "$option $substr\n";
	    } else {
		$string = "$option $option_ref->{$option}\n";
	    }	    
	}
	$self->_printline($handle, $string, $self->{"length"}); 
    }
}


package Chart::Graph::Xmgrace::Global_Options;
@ISA = qw(Chart::Graph::Xmgrace::Base_Option);

sub _init {
    my $self = shift;
    $self->{print_order} = ["version","page size","page scroll","page inout",
			    "link page","reference date","date wrap",
			    "date wrap year"];
    $self->{length} = 0;
    $self->{options} = {
			"version" => "50005",
			"page size" => ["640", "480"],
			"page scroll" => "5\%",
			"page inout" => "5\%",
			"link page" => "off",
			"reference date" => "0",
			"date wrap" => "off",
			"date wrap year" => "1950",
		       };
}

sub page_size ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"page size"} = $val;
    return 1;
}

sub page_scroll ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"page size"} = $val;
    return 1;
}

sub page_inout ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"page size"} = $val;
    return 1;
}

sub link_page ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"link page"} = $val;
    return 1;
}

sub reference_date ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"reference date"} = $val;
    return 1;
}

sub date_wrap ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"date wrap"} = $val;
    return 1;
}

sub date_wrap_year ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"date wrap year"} = $val;
    return 1;
}

package Chart::Graph::Xmgrace::World_Options;
@ISA = qw(Chart::Graph::Xmgrace::Base_Option);

my %def_world_options = (
			);

sub _init {
    my $self = shift;
    $self->{name} = "world";
    $self->{print_order} = ["xmin","xmax","ymin","ymax"];
    $self->{length} = 4;
    $self->{options} = {
			"xmin" => "0",
			"xmax" => "1",
			"ymin" => "0",
			"ymax" => "1",
		       };
}

package Chart::Graph::Xmgrace::Stack_Options;
@ISA = qw(Chart::Graph::Xmgrace::Base_Option);

my %def_stack_options = (
			 
			);

sub _init {
    my $self = shift;  
    $self->{name} = "stack";
    $self->{print_order} = ["world"];
    $self->{length} = 4;
    $self->{options} = {
			"world" => ["0","0","0","0"],
		       };
}

package Chart::Graph::Xmgrace::View_Options;
@ISA = qw(Chart::Graph::Xmgrace::Base_Option);

sub _init {
    my $self = shift;  
    $self->{name} = "view";
    $self->{print_order} = ["xmin","xmax","ymin","ymax"];
    $self->{length} = 4;
    $self->{options} = {
			"xmin" => "0.150000",
			"xmax" => "1.150000",			
			"ymin" => "0.150000",
			"ymax" => "0.850000",
		       };
}

package Chart::Graph::Xmgrace::Title_Options;
@ISA = qw(Chart::Graph::Xmgrace::Base_Option);

sub _init {
    my $self = shift;  
    $self->{name} = "title";
    $self->{print_order} = ["title","font","size","color"];
    $self->{length} = 4;
    $self->{options} = {
			"title" => '""',
			"font" => "0",
			"size" => "1.250000",
			"color" => "1",					
		       };
}

sub print($$ ) {
    my $self = shift;
    my $handle = shift;
    my $string = "";
    my $substr = "";		# for making arrays into a string

    foreach $option (@{ $self->{"print_order"} }) {
	my $option_ref = $self->{"options"};
	
	if ($option eq "title") {
	    $string = "$self->{name} \"$option_ref->{$option}\"\n";
	} else {
	    # print function handles both scalars and lists
	    if (ref($option_ref->{$option}) eq ARRAY) {
		$substr = join (", ", (@{ $option_ref->{$option} })); 
		$string = "$self->{name} $option $substr\n";
	    } else {
		$string = "$self->{name} $option $option_ref->{$option}\n";
	    }
	}	
    $self->_printline($handle, $string, $self->{"length"}); 
    }
}

package Chart::Graph::Xmgrace::Subtitle_Options;
@ISA = qw(Chart::Graph::Xmgrace::Title_Options);

sub _init {
    my $self = shift;  
    $self->{name} = "subtitle";
    $self->{print_order} = ["title","font","size","color"];
    $self->{length} = 4;
    $self->{options} = {
			"title" => '""',
			"font" => "0",
			"size" => "1.000000",
			"color" => "1",		
		       };
}

package Chart::Graph::Xmgrace::Legend_Options;
@ISA = qw(Chart::Graph::Xmgrace::Base_Option);

sub _init {
    my $self = shift;  
    $self->{name} = "legend";
    $self->{print_order} = ["status","loctype","x1","y1","box color",
			    "box pattern","box linewidth","box linestyle",
			    "box fill color","box fill pattern","font",
			    "char size","color","length","vgap","hgap",
			    "invert"];
    $self->{length} = 4;
    $self->{options} = {
			"status" => "on",
			"loctype" => "view",
			"x1" => "0.85",
			"y1" => "0.8",
			"box color" => "1",
			"box pattern" => "1",
			"box linewidth" => "1.0",
			"box linestyle" => "1",
			"box fill color" => "0",
			"box fill pattern" => "1",
			"font" => "0",
			"char size" => "1.000000",
			"color" => "1",
			"length" => "4",
			"vgap" => "1",
			"hgap" => "1",
			"invert" => "false",
		       };
}

sub box_color ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"box color"} = $val;
    return 1;
}

sub box_linewidth ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"box linewidth"} = $val;
    return 1;
}

sub box_linestyle ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"box linestyle"} = $val;
    return 1;
}

sub box_fill_color ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"box fill color"} = $val;
    return 1;
}

sub box_fill_pattern ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"box fill pattern"} = $val;
    return 1;
}

sub char_size ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"char size"} = $val;
    return 1;
}

package Chart::Graph::Xmgrace::Frame_Options;
@ISA = qw(Chart::Graph::Xmgrace::Base_Option);

sub _init {
    my $self = shift;  
    $self->{name} = "frame";
    $self->{print_order} = ["type","linestyle","linewidth","color",
			    "pattern","background color","background pattern"];
    $self->{length} = 4;
    $self->{options} = {
			"type" => "0",
			"linestyle" => "1",
			"linewidth" => "1.0",
			"color" => "1",		
			"pattern" => "1",
			"background color" => "0",
			"background pattern" => "0",
		       };
}

sub background_color ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"background color"} = $val;
    return 1;
}

sub background_pattern ($$) {
    $self = shift;
    $val = shift;
    $self->{options}->{"background pattern"} = $val;
    return 1;
}


package Chart::Graph::Xmgrace::Extra_Options;
@ISA = qw(Chart::Graph::Xmgrace::Base_Option);

sub _init {
    my $self = shift;  
    $self->{name} = "extra options";
    $self->{length} = 4;
    $self->{options} = {"extras" => undef,};	# extra options are \n delimited
}


# just dumps, verbatim, whatever the user inputs
sub print($$ ) {
    my $self = shift;
    my $handle = shift;
    my $string = "";
    my @xtra_opts;
    my $pre_sub = $self->{options}->{extras};
    if ($pre_sub) {
	@xtra_opts = split(/;\n*\s*/, $pre_sub);
    
	foreach $option (@xqtra_opts) {
	    $self->_printline($handle, "$option\n", $self->{"length"});
	}
    }
}

1;

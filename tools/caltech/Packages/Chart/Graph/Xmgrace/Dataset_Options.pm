#
# Dataset_Options.pm 
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

use Chart::Graph::Xmgrace::Base_Dataset_Option;
use Carp;

package Chart::Graph::Xmgrace::Symbol_Options;
@ISA= qw(Chart::Graph::Xmgrace::Base_Dataset_Option);

sub _init {
    my $self = shift;
    my $color = shift;

    $self->{name} = "symbol";
    $self->{print_order} = ["symbol type","size","color","pattern","fill color",
			    "fill pattern","linewidth","linestyle","char",
			    "char font","skip"];
    $self->{length} = 4;
    $self->{options} = {
			"symbol type" => "0",
			"size" => "1.000000",
			"color" => $color,
			"pattern" => "1",
			"fill color" => $color,		    
			"fill pattern" => "0",
			"linewidth" => "1.0",
			"linestyle" => "1",
			"char" => "65", 
			"char font" => "0",
			"skip" => "0",
		       };
}

sub fill_color ($$) {
    my $self = shift;
    my $val = shift;
    $self->{options}->{"fill color"} = $val;    
}

sub fill_pattern ($$) {
    my $self = shift;
    my $val = shift;
    $self->{options}->{"fill pattern"} = $val;    
}

sub symbol_char ($$) {
    my $self = shift;
    my $val = shift;
    $self->{options}->{"symbol char"} = $val;    
}

sub char_font ($$) {
    my $self = shift;
    my $val = shift;
    $self->{options}->{"char font"} = $val;    
}

sub print($$$ ) {
    my $self = shift;
    my ($handle, $set) = @_;
    my $string = "";
  
    foreach $option (@{ $self->{"print_order"} }) {
	my $option_ref = $self->{"options"};

	if ($option eq "symbol type") {
	    $string = "$set $self->{name} $option_ref->{$option}\n";
	} else {
	    $string = "$set $self->{name} $option $option_ref->{$option}\n";
	}
	    
	$self->_printline($handle, $string, $self->{"length"});
	$string = "";
    }
}

package Chart::Graph::Xmgrace::Line_Options;
@ISA= qw(Chart::Graph::Xmgrace::Base_Dataset_Option);

sub _init {
    my $self = shift;
    my $color = shift;
    $self->{name} = "line";
    $self->{print_order} = ["type","linestyle","linewidth","color","pattern"];
    $self->{length} = 4;
    $self->{options} = {
			"type" => "1", # straight, left stairs, right stairs, segs, 3-segs
			"linestyle" => "1", # solid, 7 other variations of dotted lines
			"linewidth" => "1.0",
			"color" => $color,
			"pattern" => "1",
		       };
}

package Chart::Graph::Xmgrace::Baseline_Options;
@ISA= qw(Chart::Graph::Xmgrace::Base_Dataset_Option);

sub _init {
    my $self = shift;
    my $color = shift;
    $self->{name} = "baseline";
    $self->{print_order} = ["type","status"];
    $self->{length} = 4;
    $self->{options} = {
			"type" => "0",
			"status" => "off",
		       };

}

package Chart::Graph::Xmgrace::Dropline_Options;
@ISA= qw(Chart::Graph::Xmgrace::Base_Dataset_Option);

sub _init {
    my $self = shift;
    my $color = shift;
    $self->{name} = "dropline";
    $self->{print_order} = ["status"];
    $self->{length} = 4;
    $self->{options} = {
			"status" => "off",
		       };

}

package Chart::Graph::Xmgrace::Fill_Options;
@ISA= qw(Chart::Graph::Xmgrace::Base_Dataset_Option);

sub _init {
    my $self = shift;
    my $color = shift;
    $self->{name} = "fill";
    $self->{print_order} = ["type","rule","color","pattern"];
    $self->{length} = 4;
    $self->{options} = {
			"type" => "0",
			"rule" => "0",
			"color" => $color,
			"pattern" => "1",
		       };
}
    

package Chart::Graph::Xmgrace::Avalue_Options;
@ISA= qw(Chart::Graph::Xmgrace::Base_Dataset_Option);

sub _init {
    my $self = shift;
    my $color = shift;
    $self->{name} = "avalue";
    $self->{print_order} = ["status","type","char size","font","color",
			    "rot","format","prec","prepend","append","offset"];
    $self->{length} = 4;
    $self->{options} = {
			"status" => "off",
			"type" => "2",
			"char size" => "1.000000",
			"font" => "0",
			"color" => $color, # used to be "auto"
			"rot" => "0",
			"format" => "general",
			"prec" => "3",
			"prepend" => "",
			"append" => "",
			"offset" => ["0.000000", "0.000000"],
		       };

}

sub char_size ($$) {
    my $self = shift;
    my $val = shift;
    $self->{options}->{"char size"} = $val;    
}


package Chart::Graph::Xmgrace::Errorbar_Options;
@ISA= qw(Chart::Graph::Xmgrace::Base_Dataset_Option);

sub _init {
    my $self = shift;
    my $color = shift;
    $self->{name} = "errorbar";
    $self->{print_order} = ["status","place","color","pattern","size","linewidth",
			    "linestyle","riser linewidth","riser linestyle",
			    "riser clip", "riser clip length"];
    $self->{length} = 4;
    $self->{options} = {
			"status" => "on",
			"place" => "normal",
			"color" => $color,
			"pattern" => "1",
			"size" => "1.000000",
			"linewidth" => "1.0",
			"linestyle" => "1",
			"riser linewidth" => "1.0",
			"riser linestyle" => "1",
			"riser clip" => "off",
			"riser clip length" => "0.100000",
		       };
}

sub riser_linewidth ($$) {
    my $self = shift;
    my $val = shift;
    $self->{options}->{"riser linewidth"} = $val;    
}

sub riser_linestyle ($$) {
    my $self = shift;
    my $val = shift;
    $self->{options}->{"riser linestyle"} = $val;    
}

sub riser_clip_status ($$) {
    my $self = shift;
    my $val = shift;
    $self->{options}->{"riser clip status"} = $val;    
}

sub riser_clip_length ($$) {
    my $self = shift;
    my $val = shift;
    $self->{options}->{"riser clip length"} = $val;    
}

1;






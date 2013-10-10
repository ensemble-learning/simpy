#
# Axes.pm contains options for the Axes
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
# Included objects:
#                  1) Chart::Graph::Xmgrace::Axes_Options
#

$VERSION = 3.2;

package Chart::Graph::Xmgrace::Axes;
@ISA = qw(Chart::Graph::Xmgrace::Axis_Option);
use Chart::Graph::Xmgrace::Axis;
use Carp;

sub _init {
    my $self = shift;
    $self->{print_order} = ["xaxes", "yaxes", "xaxis","yaxis",
			    "altxaxis","altyaxis"];
    $self->{length} = 4;
    $self->{options} = {
			"xaxes" => new Chart::Graph::Xmgrace::Axes_Options("xaxes"),
			"yaxes" => new Chart::Graph::Xmgrace::Axes_Options("yaxes"),
			"xaxis" => new Chart::Graph::Xmgrace::Axis("xaxis"),
			"yaxis" => new Chart::Graph::Xmgrace::Axis("yaxis"),
			"altxaxis" => new Chart::Graph::Xmgrace::Axis("altxaxis"),
			"altyaxis" => new Chart::Graph::Xmgrace::Axis("altyaxis"),
		       };
    
    # by default, the alternate axes are off
    $self->altxaxis->status("off");
    $self->altxaxis->label_options->place(["0.000000","0.105000"]);
    $self->altyaxis->status("off");
    $self->altyaxis->label_options->place(["0.000000","0.105000"]);
}

sub print ($$$$ ) {
    my $self = shift;
    my ($handle) = @_;
  
    foreach $option (@{ $self->{"print_order"} }) {
	my $option_ref = $self->{"options"};
	
	# call each objects own print function
	$option_ref->{$option}->print($handle);
    }
}

package Chart::Graph::Xmgrace::Axes_Options;
@ISA = qw(Chart::Graph::Xmgrace::Base_Option); 
use Carp;

sub _init {
    my $self = shift;
    my $type_of_axes = shift;	# xaxes, yaxes
    $self->{name} = $type_of_axes;
    $self->{print_order} = ["scale","invert"];
    $self->{length} = 4;
    $self->{options} = {
			"scale" => "Normal", # Normal, Linear, Reciprocal
			"invert" => "off",
		       };
}

1;

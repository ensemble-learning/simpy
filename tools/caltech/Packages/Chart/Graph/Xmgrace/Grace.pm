#
# Grace.pm is the grace object that holds ALL options
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

package Chart::Graph::Xmgrace::Grace;
use Chart::Graph::Xmgrace::Graph_Options;
use Chart::Graph::Xmgrace::Axes;
@ISA = qw(Chart::Graph::Xmgrace::Base_Option);

sub _init {
    my $self = shift;
    my $graph_number = shift;
    $self->{print_order} = ["global options","graph global options","world options",
			    "stack options","view options","title options",
			    "subtitle options","axes options","legend options",
			    "frame options", "extra options"];
    $self->{output_type} = undef;
    $self->{output_file} = undef;
    $self->{grace_output_file} = undef;
    $self->{auto_color} = undef;
    $self->{length} = 0;
    $self->{options} = {
			"global options" => new Chart::Graph::Xmgrace::Global_Options,
			"graph global options" => new Chart::Graph::Xmgrace::Graph_Global_Options($graph_number),
			"world options" => new Chart::Graph::Xmgrace::World_Options,
			"stack options" => new Chart::Graph::Xmgrace::Stack_Options,
			"view options" => new Chart::Graph::Xmgrace::View_Options,
			"title options" => new Chart::Graph::Xmgrace::Title_Options,
			"subtitle options" => new Chart::Graph::Xmgrace::Subtitle_Options,
			"axes options" => new Chart::Graph::Xmgrace::Axes,
			"legend options" => new Chart::Graph::Xmgrace::Legend_Options,
			"frame options" => new Chart::Graph::Xmgrace::Frame_Options,
			"extra options" => new Chart::Graph::Xmgrace::Extra_Options,
		       };
}      

sub print ($$ ) {
    my $self = shift; 
    my ($handle) = @_;
    #my @cum_data_set_objects = @{$cum_data_set_objects}; 

    print $handle "# Grace project file\n# \n";
    foreach $option (@{ $self->{"print_order"} }) {
        my $option_ref = $self->{"options"};

        # call each objects own print function
        $option_ref->{$option}->print($handle);
    }
}          

sub graph_global_options ($) {
     my $self = shift;
     return $self->{options}->{"graph global options"};
}

sub global_options ($) {
     my $self = shift;
     return $self->{options}->{"global options"};
}

sub world_options ($) {
    my $self = shift;
    return $self->{options}->{"world options"};
}

sub stack_options ($) {
    my $self = shift;
    return $self->{options}->{"stack options"};
}

sub view_options ($) {
    my $self = shift;
    return $self->{options}->{"view options"};
}

sub axes_options ($) {
    my $self = shift;
    return $self->{options}->{"axes options"};
}

sub title_options ($$) {
    my $self = shift;
    return $self->{options}->{"title options"};
}

sub subtitle_options ($) {
    my $self = shift;
    return $self->{options}->{"subtitle options"};
}

sub legend_options ($) {
    my $self = shift;
    return $self->{options}->{"legend options"};
}
sub frame_options ($) {
    my $self = shift;
    return $self->{options}->{"frame options"};
}
sub extra_options ($) {
    my $self = shift;
    return $self->{options}->{"extra options"};
}
sub output_file ($$) {
    my $self = shift;
    my $value = shift;
    $self->{output_file} = $value; 
}
sub output_type ($$) {
    my $self = shift;
    my $value= shift;
    $self->{output_type} = $value; 
}
sub grace_output_file ($$) {
    my $self = shift;
    my $value = shift;
    $self->{grace_output_file} = $value; 
}
sub autocolor ($$) {
    my $self = shift;
    my $value = shift;
    $self->{autocolor} = $value; 
}

1;

#
# Graph_Presentation_Type manages what type of graphs and what kind of presentation 
# is used. XY graph, Xy chart, BAR graph, BAR chart.
#
# contains: Chart::Graph::Xmgrace::Graph
#           Chart::Graph::Xmgrace::Chart
# 
# things to do as of 4/2000:
#           1) implement Polar graph object
#           2) implement Smith chart object
#           3) implement Fixed chart object
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

use Chart::Graph::Xmgrace::Dataset_Options;
use Chart::Graph::Xmgrace::Base_Dataset_Option;

package Chart::Graph::Xmgrace::Graph_Presentation_Type;
use Carp;

sub new {
    my $that = shift;
    my $color = shift;
    $color += 1;
    my $class = ref($that) || $that;
    my $self = {
		"XY graph" => new Chart::Graph::Xmgrace::Graph("XY", $color),
		"XY chart" => new Chart::Graph::Xmgrace::Chart("XY", $color),
		"BAR graph" => new Chart::Graph::Xmgrace::Graph("BAR",$color),
		"BAR chart" => new Chart::Graph::Xmgrace::Chart("BAR",$color),
		"Polar graph" => undef,
		"Smith chart" => undef,
		"Fixed" => undef,
               };
    bless $self, $class;
    return $self;
}         

sub XY_graph {
    my $self = shift;
    return $self->{"XY graph"};
}

sub XY_chart {
    my $self = shift;
    return $self->{"XY chart"};
}

sub AUTOLOAD {
    my $self = shift;
    my $type = ref($self) || croak "$self is not an object";
    my $name = $AUTOLOAD;
    $name =~ s/.*://; #strip fully-qualified portion
    unless (($name eq "DESTROY") or  (exists $self->{$name})) {
        croak "Can't access '$name' field in object of class $type";
    }

    if (@_) {
        return $self->{$name} = shift;
    } else {
        return $self->{$name};
    }
}

package Chart::Graph::Xmgrace::Graph;
use Carp;
@ISA = qw(Chart::Graph::Xmgrace::Base_Dataset_Option);

sub _init {
    my $self = shift;
    my ($type, $color) = @_;
    $self->{type} = $type;
    $self->{print_order} = ["hidden","type","symbol","line","baseline",
			    "dropline","fill","avalue","errorbar","comment",
			    "legend"],
    $self->{length} = 4;
    $self->{options} = {
			"hidden" => "false",
			"type" => "$type",
			"symbol" => new Chart::Graph::Xmgrace::Symbol_Options($color),
			"line" => new Chart::Graph::Xmgrace::Line_Options($color),
			"baseline" => new Chart::Graph::Xmgrace::Baseline_Options($color),
			"dropline" => new Chart::Graph::Xmgrace::Dropline_Options($color),
			"fill" => new Chart::Graph::Xmgrace::Fill_Options($color),
			"avalue" => new Chart::Graph::Xmgrace::Avalue_Options($color),
			"errorbar" => new Chart::Graph::Xmgrace::Errorbar_Options($color),
			"comment" => "",
			"legend" =>  "",
		       };

    if ($type eq "XY") {
	return 1;
    } elsif ($type eq "BAR") {
	$self->symbol->fill_pattern("1");
	$self->symbol->color("1");
	$self->line->type("0");  
	return 1;
    } else {
	carp "Warning: Invalid graph type\n";
	return 0;
    }
}

package Chart::Graph::Xmgrace::Chart;
use Carp;
@ISA = qw(Chart::Graph::Xmgrace::Base_Dataset_Option);

sub _init {
    my $self = shift;
    my ($type, $color) = @_;
    $self->{type} = $type;
    $self->{print_order} = ["hidden","type","symbol","line","baseline",
			    "dropline","fill","avalue","errorbar","comment",
			    "legend"],
    $self->{length} = 4;

    $self->{options} = {
			"hidden" => "false",
			"type" => "$type",
			"symbol" => new Chart::Graph::Xmgrace::Symbol_Options($color),
			"line" => new Chart::Graph::Xmgrace::Line_Options($color),
			"baseline" => new Chart::Graph::Xmgrace::Baseline_Options($color),
			"dropline" => new Chart::Graph::Xmgrace::Dropline_Options($color),
			"fill" => new Chart::Graph::Xmgrace::Fill_Options($color),
			"avalue" => new Chart::Graph::Xmgrace::Avalue_Options($color),
			"errorbar" => new Chart::Graph::Xmgrace::Errorbar_Options($color),
			"comment" => "",
			"legend" =>  "",
		       };

    if ($type eq "XY") {
	return 1;
    } elsif ($type eq "BAR") {
	$self->symbol->fill_pattern("1");
	$self->symbol->pattern("0");
	$self->line->type("0");  
	return 1;
    } else {
	carp "Warning: Invalid graph type\n";
	return 0;
    }
}
    
1;

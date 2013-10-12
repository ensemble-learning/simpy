#
# Dataset.pm holds essential "ingredients" to each Xmgrace dataset. 
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

package Chart::Graph::Xmgrace::Dataset;
use Carp;

sub new {
    my $that = shift;
    my $class = ref($that) || $that;
    my $self = {
		"data" => undef,  # holds a reference to data
		"options" => undef,
		"set_type" => undef, # holds set type ie, XY, BAR for the dataset
		"set_number" => undef, # holds set number
		"set_color" => undef, # holds the default set color number
		"data_format" => undef,
               };
#    $self->{set_color} = $self->{set_number} + 1;
    bless $self, $class;
    return $self;
}

sub AUTOLOAD {
    my $self = shift;
    my $type = ref($self) || croak "$self is not an object";
    my $name = $AUTOLOAD;
    $name =~ s/.*://; #strip fully-qualified portion
    unless ( ($name eq "DESTROY") or (exists $self->{$name})) {
        croak "Can't access '$name' field in object of class $type";
    }

    if (@_) {
        return $self->{$name} = shift;
    } else {                                     
        return $self->{$name};
    }
}

1;      

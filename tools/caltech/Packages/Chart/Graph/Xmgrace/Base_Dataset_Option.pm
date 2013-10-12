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

package Chart::Graph::Xmgrace::Base_Dataset_Option;
use Carp;

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {};
    bless $self, $class;
    $self->_init(@_);
    return $self;
}

# 
# 
# Subroutine: _printline() 
# 
# Description: this function will print out a line that can be  
#              read with xmgrace (with the prepended @) 
# 
# 
# 
 
sub _printline ($$$$ ) { 
    my $self = shift;
    my ($handle, $string, $length) = @_;

    if (!$string) {
	print $handle "error: no string to print!\n";
	return 0;
    }

    print $handle "@";
    print $handle ' ' x $length;
    print $handle "$string"; 

    return 1;                   # just for fun 
} 

sub print($$$ ) {
    my $self = shift;
    my ($handle, $set) = @_;
    my $string = "";
  
    foreach $option (@{ $self->{"print_order"} }) {
	my $option_ref = $self->{"options"};

	if ($option eq "status" or $option eq "in_out_status") {   
	  
	    # we first check the status of the option, whether it's on/off
	    # if it's off, we don't print it out
	    if ($option_ref->{"status"} eq "off") {
		$string = "$set $self->{name} $option_ref->{$option}\n";
		$self->_printline($handle, $string, $self->{"length"});
		last;
	    }
	    $string = "$set $self->{name} $option_ref->{$option}\n";
	    $self->_printline($handle, $string, $self->{"length"});

	} else {	

	    # we skip hidden datasets
	    if ($option eq "hidden") {
		if ($option_ref->{"hidden"} eq "true") {
		    last;
		}
	    }

	    if (!ref($option_ref->{$option})) {
		if ($option eq "comment" or $option eq "legend") {
		    $string = "$set $option \"$option_ref->{$option}\"\n";	
		} elsif (!$self->{name}) {
		    $string = "$set $option $option_ref->{$option}\n";
		} else {
		    $string = "$set $self->{name} $option $option_ref->{$option}\n";
		}
		
	    } elsif (ref($option_ref->{$option}) eq ARRAY) {
		$substr = join (", ", (@{ $option_ref->{$option} })); 
		$string = "$set $self->{name} $option $substr\n";
	    } else { # a blessed object, use own print function (which in this case is the same)
		$option_ref->{$option}->print($handle, $set);
		next;
	    }

	    $self->_printline($handle, $string, $self->{"length"});
	    $string = "";
	}
    }
}

sub AUTOLOAD {
    my $self = shift;
    my $type = ref($self) || croak "$self is not an object";
    my $name = $AUTOLOAD;
    $name =~ s/.*://; #strip fully-qualified portion
    unless (($name eq "DESTROY") or (exists $self->{options}->{$name})) {
	croak "Can't access '$name' field in object of class $type";
    }

    if (@_) {
	return $self->{options}->{$name} = shift;
    } else {
	return $self->{options}->{$name};
    }
}

1;

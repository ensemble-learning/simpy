## XrtUtils.pm is a sub-module of Graph.pm. It has all the subroutines 
## needed for the Xrt3d part of the package.
##
## $Id: XrtUtils.pm,v 1.13 2006/06/07 21:09:33 emile Exp $ $Name:  $
##
## This software product is developed by Michael Young and David Moore,
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
package Chart::Graph::XrtUtils;
use Exporter ();

@ISA = qw(Exporter);
@EXPORT = qw();
%EXPORT_TAGS = (UTILS => [qw(&_set_xrtpaths &_set_ldpath &_print_matrix 
			     &_print_array &_verify_ticks &_exec_xrt3d &_exec_xrt2d
			     &_exec_netpbm &_exec_xvfb &_try_port &_convert_raster
			     &_childpid_dead &_transfer_file)],
	       );

Exporter::export_ok_tags('UTILS');

use Carp;
use POSIX ":sys_wait_h";	# for waitpid()
use Chart::Graph::Utils qw(:UTILS); 

$cvs_Id = '$Id: XrtUtils.pm,v 1.13 2006/06/07 21:09:33 emile Exp $';
$cvs_Author = '$Author: emile $';
$cvs_Name = '$Name:  $';
$cvs_Revision = '$Revision: 1.13 $';

$VERSION = 3.2;

use strict;
use File::Basename;

use vars qw($xrt2d $xrt3d);

#
# Subroutine: set_converterpaths()
#
# Description: set paths for converter programs in particular.  This
#              subroutine can take one or two arguments and tests if
#              the required converter programs are indeed available
#              for the choosen method to convert a file from one
#              graphics format to another.
#
sub _set_converterpaths {
	my @converters = @_;

	# Loop through the list of converter seeing which are available
	foreach my $converter (@converters) {
	  if (not -x $$converter) {
		if (not $$converter = _get_path($$converter)) {
		  return(0);
		} 
	  }
	}
	return(1);
  }

#
# Subroutine: _convert_raster($plot_file, $output_file)
#
# Description: A subroutine to over see the conversion process from
#              one raster graphic format to another.  It will try
#              ImageMagick convert first and if that fails try Netpbm
#              utilities if they are available in that format.
#
sub _convert_raster {
    my $FORMAT = shift;
    my $plot_file = shift;
    my $output_file = shift;

  # First try ImageMagick as it is more robust and simpler
  if (_set_converterpaths(\$convert)) {
	if (_exec_convert($convert, $FORMAT, $plot_file, $output_file)) {
	  return(1);
	} else {
	  carp "Attempt to use ImageMagick failed, will try Netpbm."
	} 
  }	else {
	carp "No ImageMagick found, will try Netpbm."
  }

  if ($FORMAT eq 'GIF') {
	_try_netpbm_combo($xwdtopnm, $ppmtogif, $plot_file, $output_file);
  }
  elsif ($FORMAT eq 'JPG') {
	_try_netpbm_combo($xwdtopnm, $ppmtojpg, $plot_file, $output_file);
  }
  elsif ($FORMAT eq 'PNG') {
	_try_netpbm_combo($xwdtopnm, $pnmtopng, $plot_file, $output_file);
  } else {
	carp "Untrapped raster image format - XrtUtils.pm internal error";
	return(0);
  }
}

#
# Subroutine _try_netpbm_combo($xwdtopbm, pbmtotarget, $xwd_file, $target_file)
#
#
# Description: Contains the logic for testing if a combination of
#              netpbm programs can be accessed and executed to perform
#              the desired conversion.  If not, it produces the
#              appropriate error messages.  Basically, it saves a
#              batch of conditional statements that would otherwise be
#              repeated.
#
sub _try_netpbm_combo {
    my ($xwdtopbm, $pbmtotarget, $xwd_file, $target_file) = @_;
	
	if (_set_converterpaths(\$xwdtopbm, \$pbmtotarget)) {
	  if (_exec_netpbm($xwdtopbm, $pbmtotarget, $xwd_file, $target_file)) {
		return(1);
	  } else {
		carp "Failure to execute any suitable image " . 
		  "converters for create file: $target_file";
		return(0);
	  } 
	} else {
	carp "Unable to find any suitable image converters to " . 
	  "create file: $target_file";
	return(0);
  }
}

# 
# Subroutine: set_xrtpaths()
# 
# Description: set paths for external programs required by xrt()
#              if they are not defined already
#
sub _set_xrtpaths {

    my $xrtver = shift;



    if (defined($xrtver)) {
	  if ($xrtver eq "xrt2d") {
	    if (not $Chart::Graph::xrt2d = _get_path("xrt2d")) {
		  return 0;
	    }
	  }

	  if ($xrtver eq "xrt3d") {
	    if (not $Chart::Graph::xrt3d = _get_path("xrt3d")) {
		  return 0;
	    }
	  }
    }

    if (not defined($xwdtopnm)) {
	if (!($xwdtopnm = _get_path("xwdtopnm"))) {
	    return 0;
	}
    }

    if (not defined($xvfb)) {
	if (not $xvfb = _get_path("Xvfb")) {
	    return 0;
	}
    }

    # make sure /usr/dt/lib is in the library path
    _set_ldpath("/usr/dt/lib");

    return 1;
}

#
# Subroutine: set_ldpath()
#
# Description: Xvfb has trouble finding libMrm, so we have to add
#              /usr/dt/lib to LD_LIBRARY_PATH
#

sub _set_ldpath {
    my ($libpath) = @_;
    
    if (not defined($ENV{LD_LIBRARY_PATH})) {
	$ENV{LD_LIBRARY_PATH} = "$libpath";
	return 1;
    }

    my @ldpath = split (/:/, $ENV{LD_LIBRARY_PATH});

    # make sure library path isn't already defined
    foreach my $i(@ldpath){
	if ($i eq $libpath) {
	    return 1;
	}
    }

    # add library path to LD_LIBRARY_PATH
    $ENV{LD_LIBRARY_PATH} = "$libpath:$ENV{LD_LIBRARY_PATH}";
    return 1;
}

# 
# Subroutine: print_matrix() 
# 
# Description: print out all the elements 
#              in a X by Y  matrix, row by row
#

sub _print_matrix {
    my ($handle, @matrix) = @_;
    
    foreach my $row (@matrix){
	foreach my $i (@{$row}){
	    print $handle "$i\t";
	}
	print $handle "\n";
    }
    return 1;
}


#
# Subroutine:  _transfer_file($handle, $data_filename)
#
# Description: open file $data_filename.  Read the contents
# and write it into the command file tab delimited.  Don't
# assume data was tab delimited to be safe.
# 
sub  _transfer_file {
  my $handle = shift;
  my $data_filename = shift;
  my $data;
  my @elements;

  unless(open(DATAHDL, $data_filename)) {
	carp "Unable to open data file: $data_filename for reading";
	return(0);
  }
  while (defined($data = <DATAHDL>)) {
	chomp($data);
	@elements = split(/\s+/, $data);
	foreach my $element (@elements) {
	  print $handle $element, "\t";
	}
	print $handle  "\n";
  }
  unless(close(DATAHDL)) {
    carp "Unable to close data file: $data_filename after reading";
  }
  return(1);
}

# 
# Subroutine: print_array()
# 
# Description:  print out each element of array, one per line
#

sub _print_array {
    my ($handle, @array) = @_;
    my $i;
    
    foreach $i (@array) {
	print $handle "$i\n";
    }
    return 1;
}

# 
# Subroutine: verify_ticks();
#   
# Description: check that the number of tick labels is the same
#              as the number of xy rows and columns. we can only have
#              as many ticks as the number of rows or columns
#              we make this subroutine so that the calling subroutine
#              is kept cleaner.

sub _verify_ticks {
    my ($cnt, $ticks_ref) = @_;

    # if no ticks are given then just
    # give the xrt binary "1, 2,..."
    if (not defined($ticks_ref)) {
	my @def_ticks;
	for (my $i = 0; $i < $cnt; $i++) {
	    $def_ticks[$i] = $i + 1;
        }
	$ticks_ref = \@def_ticks;
    }

    my $tick_cnt = @{$ticks_ref};

    if ($cnt ne $tick_cnt){
	carp "number of tick labels must equal the number of xy rows and columns";
	return 0;
    }
    return 1;
}

# 
# Subroutine: exec_xrt3d()
#
# Description: execute the xrt3d program on the command file.
#              xrt3d generates a xwd file.
# 
sub _exec_xrt3d {
    my ($command_file) = @_;
    my ($output);
    my ($childpid, $port);
    my $display_env = $ENV{DISPLAY};
    my $status;

    if ($Chart::Graph::use_xvfb) {
	# start the virtual X server
	($childpid, $port) = _exec_xvfb();
	$status = system("$Chart::Graph::xrt3d -display :$port.0 < $command_file");
    } else {
	# use the local X server
	# warning: colors might be messed up 
	# depending on your current setup
	$status = system("$Chart::Graph::xrt3d -display $display_env < $command_file");
    }

    #my $status = system("$xrt -display :$port.0 < $command_file");
    if (not _chk_status($status)) {
	return 0;
    }

    if ($Chart::Graph::use_xvfb) {
	kill('KILL', $childpid);
    }

    return 1;
}

# 
# Subroutine: exec_xrt2d()
#
# Description: execute the xrt2d program on the command file.
#              xrt2d generates a xwd file.
#
sub _exec_xrt2d {
    my ($command_file, $options) = @_;
    my ($output);
    my ($childpid, $port);
    my $display_env = $ENV{DISPLAY};
    my $status;

    if ($Chart::Graph::use_xvfb) {
	# start the virtual X server
	($childpid, $port) = _exec_xvfb();
	printf STDERR "\tXRT is $Chart::Graph::xrt2d\n";
	my $status = system("$Chart::Graph::xrt2d -display ipn:$port.0 < $command_file $options");
    } else {
	# use the local X server
	# warning: colors might be messed up 
	# depending on your current setup
	$status = system("$Chart::Graph::xrt2d -display $display_env < $command_file $options");
    }

    if (not _chk_status($status)) {
	return 0;
    }

    if ($Chart::Graph::use_xvfb) {
	kill('KILL', $childpid);
    }

    return 1;
}

# 
# Subroutine: exec_convert 
# 
#
# Description: Use the Imagemagick 'convert' utility to convert the xwd
#              file into any one of the other common raster image
#              formats used commonly in web page production.
#
sub _exec_convert {
    my ($convert, $FORMAT, $xwd_file, $target_file) = @_;
    my ($status);


    if ($Chart::Graph::debug) {
	$status = system(join('', "$convert -verbose $xwd_file ", 
                              $FORMAT, ":$target_file"));
    } else {
	$status = system(join('', "( $convert $xwd_file ",  $FORMAT, 
                                  ":$target_file; ) 2> /dev/null"));
    }
    
    if (not _chk_status($status)) {
	return 0;
    }
    return 1;
}
# 
# Subroutine: _exec_netpbm
# 
#
# Description: Convert a raster image using the older utilities now
#              collected under the name 'netpbm.'  Note that not all
#              conversions are commonly included wiht all UNIX
#              distributions so that while older conversions such as
#              'xwd' -> 'gif' are likely to work, others such as
#              conversions to 'png' may not without downloading new
#              software.
#
#              The conversion strategy always involves a pipe from the
#              X-windows 'xwd' format to some sort 'pbm' format and
#              then from that universal format into the target format.
#              For this reason, it is more prone to machine
#              architecture issues and other errors.
#
sub _exec_netpbm {
    my ($xwdtopbm, $pbmtotarget, $xwd_file, $target_file) = @_;
    my ($status);
    
    if ($Chart::Graph::debug) {
	$status = system("$xwdtopbm $xwd_file | $pbmtotarget > $target_file");
    } else {
	$status = system(join('', "( $xwdtopbm -quiet $xwd_file | ", 
                              "$pbmtotarget -quiet > $target_file; ) ",
                              "2> /dev/null"));
    }
    
    if (not _chk_status($status)) {
	return 0;
    }
    return 1;
}

# 
# Subroutine: exec_xvfb()
#
# Description:  this starts the vitualX server(X is required by xrt, so 
#               we fake out xrt with Xvfb, for speed and compatability)
#
#
sub _exec_xvfb {
    my $port = 99;
    my $childpid;
    my $sleep_time = 1;
    my $try_count = 0;
    my $trialnumber;
    my $childpid_status;

    # starting with port 100, we try to start
    # the virtual server until we find an open port
    # because of the nature of the virtual x server
    # we use, in order to know if we have found an 
    # open port, we have to sleep.
    # we check the pid of the virtual x process we started
    # and see if it died or not.

    while ($childpid_status = _childpid_dead($childpid)) {
	$port++;
	$try_count++;
	if ($try_count > 10) {
	    die "Error: Failed too many times\n";
	}
	$trialnumber = _number_to_eng($try_count);
	print STDERR "*** $trialnumber try ***" unless (not $Chart::Graph::debug);
	$childpid = _try_port($port);	
	sleep($sleep_time);
    }
    print STDERR " SUCCESS!!!\n" unless (not $Chart::Graph::debug);

    # save the childpid so we can stop the virtual server later
    # save the $port so we can tell xrt where the virtual server is.
    return ($childpid, $port);
}
# 
# Subroutine: try_port();
#
# Description:  will try to start Xvfb on specified port
sub _try_port {

    my ($port) = @_;
    my ($childpid);
   
    #fork a process
    if (not defined($childpid = fork())){
	# the fork failed
	carp "cannot fork: $!";
	return 0;
    } elsif ($childpid == 0) {
	# we are in the child process
	if ($Chart::Graph::debug) {
	    if (not exec "$xvfb :$port") {
		die "can't do $xvfb :$port: $!\n";
	  }	
	}
	else {
	    if (not exec "$xvfb :$port 2> /dev/null") {
		die "can't do $xvfb :$port 2> /dev/null: $!\n";
	    }
	} 

	die "should never reach here\n";

    } else {
	# we are in the parent, return the childpid
	# so we can kill it later.
	return $childpid;
    }
    
}

# 
# Subroutine: childpid_dead
# 
# Description: check to see if a PID has died or not
#
#
sub _childpid_dead {
    my ($childpid) = @_;
    
    if (not defined($childpid)) {
	return 1;
    }

    # WNOHANG: waitpid() will not suspend execution  of
    # the  calling  process  if  status is not
    # immediately available  for  one  of  the
    # child processes specified by pid.
    return waitpid($childpid, &WNOHANG);
}

1;

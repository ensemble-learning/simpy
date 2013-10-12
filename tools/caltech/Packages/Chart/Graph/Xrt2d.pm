## Xrt2d.pm is a sub-module of Graph.pm. It has all the subroutines 
## needed for the Xrt2d part of the package.
##
## $Id: Xrt2d.pm,v 1.26 2006/06/07 21:09:33 emile Exp $ $Name:  $
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
package Chart::Graph::Xrt2d;
use Exporter ();

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(&xrt2d);

use FileHandle;			# to create generic filehandles
use Carp;			# for carp() and croak()
#use POSIX ":sys_wait_h";	# for waitpid()
use Chart::Graph::Utils qw(:UTILS);	# get global subs and variables
use Chart::Graph::XrtUtils qw(:UTILS); # get Xrt subroutines

$cvs_Id = '$Id: Xrt2d.pm,v 1.26 2006/06/07 21:09:33 emile Exp $';
$cvs_Author = '$Author: emile $';
$cvs_Name = '$Name:  $';
$cvs_Revision = '$Revision: 1.26 $';

$VERSION = 3.2;

use strict;

#
# xrt graphing package
#

my %def_xrt_global_opts = (
			"output file" => "untitled-xrt2d.gif",
			"output type" => "gif",
			"x-axis title" => "x-axis",
			"y-axis title" => "y-axis",
			"set labels" => undef,
			"point labels" => undef,
			"header" => [],
			"footer" => [],
			"misc labels" => undef,
			"invert" => 0,
			"x time" => 0, # x axis labels are time indicators
			"style" => "bar",
		    );

my %def_xrt_data_opts = (
			"color" => undef,
		    );

#
#
# Subroutine: xrt2d()
#
# Description: this is the main function you will be calling from
#              our scripts. please see 
#              www.caida.org/Tools/Graph/ for a full description
#              and how-to of this subroutine
#

sub xrt2d {
    my ($user_global_opts_ref, @data_sets) = @_;
    my (%global_opts, %data_opts);
    my $data_set_ref;
    
    # variables to be written to the command file
    my ($plot_file, $x_axis, $y_axis);
    my ($header, $footer, $misc_labels, $set_labels, $point_labels);
    my ($x_cnt, $y_cnt, $hdr_cnt, $ftr_cnt, $label_cnt);
    my ($output_file, $output_type);
    my $xrt_options = "";
    my $x_timestamps = 0;

    _make_tmpdir("_Xrt2d_");

    # set paths for external programs
    if (not _set_xrtpaths("xrt2d")) {
	_cleanup_tmpdir();
	return 0;
    }
    
    
    # check first arg for hash
    if (ref($user_global_opts_ref) ne "HASH") {
	carp "Global options must be a hash.";
	_cleanup_tmpdir();
	return 0;
    }

    # call to combine user options with default options 
    %global_opts = _mesh_opts($user_global_opts_ref, \%def_xrt_global_opts);
    
    # check for values in command file
    while (my ($key, $value) = each %global_opts) {
	
	if ($key eq "output file") {
	    if (defined($value)) {
		$output_file = $value;
		unless (defined $global_opts{"output type"}) {
		    carp "Must have an output type defined";
		    _cleanup_tmpdir();
		    return 0;
		}

		# If the file is PostScript ... what XRT makes is PostScript
		if ($global_opts{"output type"} eq "ps") {
		    $plot_file = _make_tmpfile("plot", "ps");
		  } 
		# For all raster formats XRT starts out with X-Windows XWD format.
		elsif (($global_opts{"output type"} eq "gif") or
			   ($global_opts{"output type"} eq "xwd") or
			   ($global_opts{"output type"} eq "png") or
			   ($global_opts{"output type"} eq "jpg")
			  ) {
		  $plot_file = _make_tmpfile("plot", "xwd");
		} else {
		  # Default is XWD
		  carp "Unknown output type, defaulting to xwd";
		  $plot_file = _make_tmpfile("plot", "xwd");
		}
	  }
	}
	
	if ($key eq "x-axis title") {
	    if(defined($value)) {
		$x_axis = $value;
	    }
	}
	
	if ($key eq "y-axis title") {
	    if(defined($value)) {
		$y_axis = $value;
	    }
	}
	
	if ($key eq "header") {
	    if(defined($value)) {
		$header = $value;
	    }
	}
	
	if ($key eq "footer") {
	    if(defined($value)) {
		$footer = $value;
	    }
	}

	if ($key eq "misc labels") {
	    if(defined($value)) {
		$misc_labels = $value;
	    }
	}

	if ($key eq "set labels") {
	    if(defined($value)) {
		$set_labels = $value;
	    }
	}

	if ($key eq "point labels") {
	    if(defined($value)) {
		$point_labels = $value;
	    }
	}

	if ($key eq "invert") {
	    if(defined($value) and $value) {
		$xrt_options .= " -xrm '*xrtInvertOrientation: True'";
	    }
	}

	if ($key eq "x time") {
	    if(defined($value) and $value) {
		$xrt_options .= " -xrm '*xrtXAnnotationMethod: ANNOTIMELABELS'";
		$x_timestamps = 1;
	    }
	}

	if ($key eq "style") {
	    if(defined($value)) {
		if ($value eq "stackedbar") {
		   $xrt_options .= " -xrm '*xrtType: TYPESTACKINGBAR'";
		} elsif ($value eq "bar") {
		   $xrt_options .= " -xrm '*xrtType: TYPEBAR'";
		} elsif ($value eq "pie") {
		   $xrt_options .= " -xrm '*xrtType: TYPEPIE'";
		} elsif ($value eq "area") {
		   $xrt_options .= " -xrm '*xrtType: TYPEAREA'";
		} elsif ($value eq "stackedarea") {
		   $xrt_options .= " -xrm '*xrtLegendReversed: True'";
		   $xrt_options .= " -xrm '*xrtType: TYPEAREA'" .
			" -xrm '*xrtIsStacked: True'";
		} else {
		   carp "No graph type defined, defaulting to bar";
		   $xrt_options .= " -xrm '*xrtType: TYPEBAR'";
		}
	    }
	}
    }
    # Only reverse legend if NOT inverting bars.
    if (not $global_opts{"invert"} and $global_opts{"style"} eq "stackedbar") {
	$xrt_options .= " -xrm '*xrtLegendReversed: True'";
    }
    
    # get the number of columns and number of rows
    
    $x_cnt = @{$point_labels};
    $y_cnt = @{$set_labels};
    
    # because xrt allows multiline headers
    # get the length of the header array
    # each line of the header is one index
    # in the array
    $hdr_cnt = $#{$header} + 1;
    $ftr_cnt = $#{$footer} + 1;
    $label_cnt = $#{$misc_labels} + 1;

    my @colors;
    my @all_data;
    foreach $data_set_ref (@data_sets) {
	if (ref($data_set_ref) ne "ARRAY") {
	    carp "Data set must be an array";
	    _cleanup_tmpdir();
	    return 0;
	}
	my ($user_data_opts_ref, $data_ref) = @$data_set_ref;
	my (%data_opts);
	
	## check first arg for hash
	if (ref($user_data_opts_ref) ne "HASH") {
	    carp "Data options must be a hash.";
	    _cleanup_tmpdir();
	    return 0;
	}

	push @all_data, $data_ref;

	# call to combine user options with default options   
	%data_opts = _mesh_opts($user_data_opts_ref, \%def_xrt_data_opts);

	# write data options to command file
	while (my ($key, $value) = each %data_opts) {
	    if ($key eq "color") {
		if (defined $value) {
		    push @colors, $value;
		}
	    }
	}
    }	

    if (@colors) {
	my $resource_string;
	foreach my $color (@colors) {
	    $resource_string .=
		    "(LpatSolid FpatSolid \"$color\" 1 PointNone \"$color\" 7)";
	}
	$xrt_options .= " -xrm '*xrtDataStyles: ($resource_string)'";
    }
    ##
    ## print command file using this format (for bar graphs)
    ## data separated by tabs
    ##
    
    # output.file						
    # x_cnt (number of points of data)
    # y_cnt (number of sets of info per point)
    # point1title set1 set2 ... sety
    # point2title set1 set2 ...
    # .
    # .
    # pointxtitle set1 set2 ... sety
    # Number of header lines (multiple header lines available)
    # header1
    # header2
    # ...
    # Number of header lines (multiple header lines available)
    # foot1
    # foot2
    # ...
    # Title of x-axis
    # Title of y-axis
    # label_cnt (number of extra data labels)
    # label1title point# set#
    
    # create command file and open file handle 
    my $command_file = _make_tmpfile("command");
    my $handle = new FileHandle;
    if (not $handle->open(">$command_file")) {
	carp "could not open $command_file";
	_cleanup_tmpdir();
	return 0;
    }

    print $handle "$plot_file\n";
    print $handle "$x_timestamps\n";
    print $handle "$x_cnt\n";
    print $handle "$y_cnt\n";
    _print_array($handle, @{$set_labels});
    _print_array($handle, @{$point_labels});
    _print_matrix($handle, @all_data);
    print $handle "$hdr_cnt\n";
    _print_array($handle, @{$header});
    print $handle "$ftr_cnt\n";
    _print_array($handle, @{$footer});	
    print $handle "$x_axis\n";
    print $handle "$y_axis\n";
    print $handle "$label_cnt\n";
    _print_matrix($handle, @$misc_labels);
    $handle->close();

    # call xrt and convert file to gif
    if (not _exec_xrt2d($command_file, $xrt_options)) {
	_cleanup_tmpdir();
	return 0;
    }

    my $graph_format = $global_opts{"output type"};
    if ($graph_format eq "ps") {
	  if (not _chk_status(system("cp $plot_file $output_file"))) {
	    _cleanup_tmpdir();
	    return 0;
	  }
    } elsif ($graph_format eq "xwd") {
	  if (not _chk_status(system("cp $plot_file $output_file"))) {
	    _cleanup_tmpdir();
	    return 0;
	  }
	} else {
	  if(not _convert_raster($graph_format, $plot_file, $output_file)) {
	    _cleanup_tmpdir();
	    return 0;
	  }
    }

    _cleanup_tmpdir();
    return 1;
}

sub _xrt_data_check {
    my ($user_data_opts_ref, @data) = @_;
    my (%data_opts);
    
    my ($result, $color);
    
    ## check first arg for hash
    if (ref($user_data_opts_ref) ne "HASH") {
	carp "Data options must be a hash.";
	return 0;
    }

    # call to combine user options with default options   
    %data_opts = _mesh_opts($user_data_opts_ref, \%def_xrt_data_opts);

    # write data options to command file
    while (my ($key, $value) = each %data_opts) {
	if ($key eq "color") {
	    $color = $value;
	}
    }
}

# # 
# # Subroutine: set_xrtpaths()
# # 
# # Description: set paths for external programs required by xrt()
# #              if they are not defined already
# #
# sub _set_xrtpaths {

#     my $xrt = shift;

#     if (not defined($ppmtogif)) {
# 	if (not $ppmtogif = _get_path("ppmtogif")) {
# 	    return 0;
# 	}
#     }

#     if (not defined($xrt)) {
# 	if (not $xrt = _get_path("graph_2d")) {
# 	    return 0;
# 	}
#     }

#     if (not defined($xwdtopnm)) {
# 	if (!($xwdtopnm = _get_path("xwdtopnm"))) {
# 	    return 0;
# 	}
#     }

#     if (not defined($xvfb)) {
# 	if (not $xvfb = _get_path("Xvfb")) {
# 	    return 0;
# 	}
#     }

#     # make sure /usr/dt/lib is in the library path
#     _set_ldpath("/usr/dt/lib");

#     return 1;
# }

# #
# # Subroutine: set_ldpath()
# #
# # Description: Xvfb has trouble finding libMrm, so we have to add
# #              /usr/dt/lib to LD_LIBRARY_PATH
# #

# sub _set_ldpath {
#     my ($libpath) = @_;
    
#     if (not defined($ENV{LD_LIBRARY_PATH})) {
# 	$ENV{LD_LIBRARY_PATH} = "$libpath";
# 	return 1;
#     }

#     my @ldpath = split (/:/, $ENV{LD_LIBRARY_PATH});

#     # make sure library path isn't already defined
#     foreach my $i(@ldpath){
# 	if ($i eq $libpath) {
# 	    return 1;
# 	}
#     }

#     # add library path to LD_LIBRARY_PATH
#     $ENV{LD_LIBRARY_PATH} = "$libpath:$ENV{LD_LIBRARY_PATH}";
#     return 1;
# }

# # 
# # Subroutine: print_matrix() 
# # 
# # Description: print out all the elements 
# #              in a X by Y  matrix, row by row
# #

# sub _print_matrix {
#     my ($handle, @matrix) = @_;
    
#     foreach my $row (@matrix){
# 	foreach my $i (@{$row}){
# 	    print $handle "$i\t";
# 	}
# 	print $handle "\n";
#     }
#     return 1;
# }
# # 
# # Subroutine: print_array()
# # 
# # Description:  print out each element of array, one per line
# #

# sub _print_array {
#     my ($handle, @array) = @_;
#     my $i;
    
#     foreach $i (@array) {
# 	print $handle "$i\n";
#     }
#     return 1;
# }

# # 
# # Subroutine: verify_ticks();
# #   
# # Description: check that the number of tick labels is the same
# #              as the number of xy rows and columns. we can only have
# #              as many ticks as the number of rows or columns
# #              we make this subroutine so that the calling subroutine
# #              is kept cleaner.

# sub _verify_ticks {
#     my ($cnt, $ticks_ref) = @_;

#     # if no ticks are given then just
#     # give the xrt binary "1, 2,..."
#     if (not defined($ticks_ref)) {
# 	my @def_ticks;
# 	for (my $i = 0; $i < $cnt; $i++) {
# 	    $def_ticks[$i] = $i + 1;
#         }
# 	$ticks_ref = \@def_ticks;
#     }

#     my $tick_cnt = @{$ticks_ref};

#     if ($cnt ne $tick_cnt){
# 	carp "number of tick labels must equal the number of xy rows and columns";
# 	return 0;
#     }
#     return 1;
# }

# # 
# # Subroutine: exec_xrt()
# #
# # Description: execute the xrt program on the command file.
# #              xrt generates a xwd file.
# # 
# sub _exec_xrt {
#     my ($command_file, $options) = @_;
#     my ($output);
#     my ($childpid, $port);
#     my $display_env = $ENV{DISPLAY};
#     my $status;

#     if ($use_xvfb) {
# 	# start the virtual X server
# 	($childpid, $port) = _exec_xvfb();
# 	printf STDERR "\tXRT is $xrt\n";
# 	my $status = system("$xrt -display ipn:$port.0 < $command_file $options");
#     } else {
# 	# use the local X server
# 	# warning: colors might be messed up 
# 	# depending on your current setup
# 	$status = system("$xrt -display $display_env < $command_file $options");
#     }

#     if (not _chk_status($status)) {
# 	return 0;
#     }

#     kill('KILL', $childpid);
#     return 1;
# }
	
# # 
# # Subroutine: exec_xwdtogif 
# # 
# # Description: convert the xwd file generated by xrt into a gif
# #              this is a 2-step process. the xwd must be converted into 
# #              a pnm and then into a gif.
# sub _exec_xwdtogif {
#     my ($xwd_file, $gif_file) = @_;
#     my ($status);
    
#     if ($Chart::Graph::debug) {
# 	$status = system("$xwdtopnm $xwd_file | $ppmtogif > $gif_file");
#     } else {
# 	$status = system("( $xwdtopnm $xwd_file | $ppmtogif > $gif_file; ) 2> /dev/null");
#     }
    
#     if (not _chk_status($status)) {
# 	return 0;
#     }
#     return 1;
# }

# # 
# # Subroutine: exec_xvfb()
# #
# # Description:  this starts the vitualX server(X is required by xrt, so 
# #               we fake out xrt with Xvfb, for speed and compatability)
# #
# #
# sub _exec_xvfb {
#     my $port = 99;
#     my $childpid;
#     my $sleep_time = 1;
    

#     # starting with port 100, we try to start
#     # the virtual server until we find an open port
#     # because of the nature of the virtual x server
#     # we use, in order to know if we have found an 
#     # open port, we have to sleep.
#     # we check the pid of the virtual x process we started
#     # and see if it died or not.
#     while (_childpid_dead($childpid)) {
# 	$port++;
# 	$childpid = _try_port($port);	
# 	sleep($sleep_time);
#     }
  
#     # save the childpid so we can stop the virtual server later
#     # save the $port so we can tell xrt where the virtual server is.
#     return ($childpid, $port);
# }
# # 
# # Subroutine: try_port();
# #
# # Description:  will try to start Xvfb on specified port
# sub _try_port {

#     my ($port) = @_;
#     my ($childpid);
    
#     #fork a process
#     if (not defined($childpid = fork())){
# 	# the fork failed
# 	carp "cannot fork: $!";
# 	return 0;
#     } elsif ($childpid == 0) {
# 	# we are in the child process
# 	if ($Chart::Graph::debug) {
# 	    exec "$xvfb :$port";
# 	}
# 	else {
# 	    exec "exec $xvfb :$port 2> /dev/null";
# 	}

# 	die "should never reach here\n";
#     } else {
# 	# we are in the parent, return the childpid
# 	# so re can kill it later.
# 	return $childpid;
#     }
    
# }

# # 
# # Subroutine: childpid_dead
# # 
# # Description: check to see if a PID has died or not
# #
# #
# sub _childpid_dead {
#     my ($childpid) = @_;
    
#     if (not defined($childpid)) {
# 	return 1;
#     }

#     # WNOHANG: waitpid() will not suspend execution  of
#     # the  calling  process  if  status is not
#     # immediately available  for  one  of  the
#     #  child processes specified by pid.
#     return waitpid($childpid, &WNOHANG);
# }

1;


__END__

=head1 NAME

Chart::Graph::Xrt2d

=head1 SYNOPSIS

 #Include module
 use Chart::Graph::Xrt2d qw(xrt2d);

 # Function call
 xrt2d(\%options,
       [\%data_options1, \@data_set1],
       [\%data_options2, \@data_set2],
	.
	.
      );


=head1 DESCRIPTION

This module is unmaintained, it worked with Sitraka's XRT, and hasn't been
tested against newer versions.

Sitraka (now Quest) makes a number of graphics packages for UNIX systems.  XRT is
a Motif-based commercial software product that has been adapted by
CAIDA using a combination of C drivers and Perl function I<xrt2d()>.
The Perl function I<xrt2d()> provides access to the two dimensional
graphing capabilities of XRT from Perl.  To access the three dimensional
graphing using XRT, use I<xrt3d()> also supplied in the
I<Chart::Graph> package.

=head1 ARGUMENTS

The options to I<xrt2d()> are listed below.  Additional control over the
resulting graph is possible by using the XRT application itself once
the graph has been created.

 +--------------------------------------------------------------------------+
 |                                OPTIONS                                   |
 +----------------+--------------------------+------------------------------+
 | Name           |  Options                 | Default                      |
 +----------------+--------------------------+------------------------------+
 |"output file"   |  (set your own)          | "untitled-xrt2d.gif"         |
 |"output type"   |  "ps","xwd", "png", "jpg"| "xwd"                        |
 |"x-axis title"  |  (set your own)          | "x-axis"                     |
 |"y-axis title"  |  (set your own)          | "y-axis"                     |
 |"set labels"    |  (set your own bar chart | none                         |
 |                |   labels for each set of |                              |
 |                |   data)                  |                              |
 |"point labels"  |  (set your own labels for| none                         |
 |                |   bars themselves)       |                              |
 |"misc labels"   |  (misc annotation that   | none                         |
 |                |   can be added to bar    |                              |
 |                |   chart)                 |                              |
 |"x time"        |  timescale if appropriate| "0"                          |
 |"invert"        |  run bars horizontally   | "0" (vertically)             |
 |                |  instead of vertically.  |                              |
 |                |  1 = horizonally.        |                              |
 |"header"        |  (set your own - can     | "header"                     |
 |                |   be multiple lines)     |                              |
 |"footer"        |  (set your own - can     | "footer"                     |
 |                |   be multiple lines)     |                              |
 |"style"         |  Style of chart - types  | "bar"                        |
 |                |  include: bar, pie, area |                              |
 |                |  stackedbar, stackedarea |                              |
 +----------------+--------------------------+------------------------------+

The I<xrt2d()> function only accepts data in one form.  C<\@data_sets>:
a one dimensional array with a prefix of data options.: 
C<[[\%data1_opts, \@data1], [\%data2_opts, \@data2], [\%data3_opts, \@data3]]>
-- see example for more details on the syntax.

The data options are listed below.

 +--------------------------------------------------------------------------+
 |                             DATA OPTIONS                                 |
 +----------------+--------------------------+------------------------------+
 | Name           |  Options                 | Default                      |
 +----------------+--------------------------+------------------------------+
 | "color"        | Any valid web page color | none                         |
 +----------------+--------------------------+------------------------------+

=head2 DETAILS ON GRAPHICS CONVERTER OPTIONS

The xrt package supports only two graphics formats internally:
Postscript and the X windows format XWD.  Additional raster graphics
formats are supported with Chart::Graph by using one of two graphics
converter packages: I<Imagemagick> and I<Netpbm>.

If you need to install a converter package, I<Imagemagick>
I<http://www.imagemagick.org/> is probably preferable
simply for its comparatively simplicity.  It uses one program
I<convert> for all of it's conversion needs, so it is easy to manage
and simple for Chart::Graph to use.  Many UNIX systems come with some
collection of the I<Netpbm> utilities already installed, thus users
may be able to start using Chart::Graph without adding any additional
converters.  Alas, it is unlikely any distributions would include all
the converters for the newest graphics formats used by Chart::Graph.
In that case it may still preferable to use I<Imagemagick> simply for
the sake of avoiding installing over 80 utilities that come with
current distributions of I<Netpbm>.  For more information on the
current distribution of I<Netpbm> go to the current website at:
I<http://netpbm.sourceforge.net/>

The xrt package also allows for multiple header and footers with each
graph.  As a result, instead of just the usual string, an array
reference containing the multiple strings for the header and footer
text.

=head1 EXAMPLES

The following two examples show Chart::Graph::Xrt2d in different roles
and producing different styles of output.

=head2 HORIZONTAL BAR GRAPH

The following example creates a two, two dimensional bars charts of
fictitious stock data from rival restaurants that is displayed in a
single graphic file F<xrt2d-1.jpg>.

 #make sure to include Chart::Graph
 use Chart::Graph::Xrt2d qw(xrt2d);

 # Call to xrt2d with two data sources.
 xrt2d({"output file" => "xrt2d-1.jpg",
        "output type" => "jpg",
        "set labels"=> ["Joe's", "Ralph's"],
        # Flip graph from vertical to horizontal
	"invert" => 1,
	"point labels" => ["Jan/Feb", "Mar/Apr", "May/Jun", "Jul/Aug",
                         "Sep/Oct", "Nov/Dec"],
        "x-axis title" => "Month's tracked",
        "y-axis title" => "Stock prices for Rival restaurant chains"
       },
       [{"color" => "MistyRose"}, ["8", "13", "20", "45", "50", "100"]],
       [{"color" => "#000000"},   ["75", "50", "25", "25", "50", "75"]]
    )
 )

=for html
<p><center><img src="http://www.caida.org/tools/utilities/graphing/xrt2d-1.jpg"></center></p>
<p><center><em>xrt2d-1.jpg</em></center></p>

=head2 VERTICAL BAR GRAPH

The following example is simply a more elaborate instance of the first
example. which produces the file F<xrt2-2.gif>.  Note the relationship
between sets and points.  Accidentally reversing the order will cause
unpredicable results.

 #make sure to include Chart::Graph
 use Chart::Graph::Xrt2d qw(xrt2d);

 xrt2d({"output file" => "xrt2d-2.gif",
			 "output type" => "gif",
		"set labels" => ["set1", "set2", "set3", "set4"],
		"point labels" => ["point1", "point2", "point3"]},
                # Each entry here corresponds to a set
		[{"color" => "MistyRose"}, ["15", "23", "10"]],
		[{"color" => "#0000FF"}, ["13", "35", "45"]],
		[{"color" => "#00FF00"}, ["15", "64", "24"]],
		[{"color" => "Navy"}, ["18", "48", "32"]],
      );

=for html
<p><center><img src="http://www.caida.org/tools/utilities/graphing/xrt2d-2.gif"></center></p>
<p><center><em>xrt2d-2.gif</em></center></p>


=head1 MORE INFO

For more information on XRT:

 http://www.quest.com/xrt_pds/

=head1 CONTACT

Send email to graph-dev@caida.org is you have problems, questions,
or comments. To subscribe to the mailing list send mail to
graph-dev-request@caida.org with a body of "subscribe your@email.com"

=head1 AUTHOR

 CAIDA Perl development team (cpan@caida.org)

=cut

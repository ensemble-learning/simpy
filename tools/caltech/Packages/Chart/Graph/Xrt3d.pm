## Xrt3d.pm is a sub-module of Graph.pm. It has all the subroutines 
## needed for the Xrt3d part of the package.
##
## $Id: Xrt3d.pm,v 1.30 2006/06/07 21:09:33 emile Exp $ $Name:  $
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
package Chart::Graph::Xrt3d;
use Exporter ();

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(&xrt3d);

use FileHandle;			# to create generic filehandles
use Carp;			# for carp() and croak()
use Chart::Graph::Utils qw(:UTILS);	# get global subs and variables
use Chart::Graph::XrtUtils qw(:UTILS);

$cvs_Id = '$Id: Xrt3d.pm,v 1.30 2006/06/07 21:09:33 emile Exp $';
$cvs_Author = '$Author: emile $';
$cvs_Name = '$Name:  $';
$cvs_Revision = '$Revision: 1.30 $';

$VERSION = 3.2;

use strict;

#
# xrt graphing package
#


my %def_xrt_global_opts;		# xrt specific globals

%def_xrt_global_opts = (
			"output file" =>  "untitled-xrt3d.gif",
			"output type" => "gif",
			"x-axis title" => "x-axis",
			"y-axis title" => "y-axis",
			"z-axis title" => "z-axis",
			"x-min" => "0",
			"y-min" => "0",
			"x-step" => "1",
			"y-step" => "1",
			"x-ticks" => undef,
			"y-ticks" => undef,
			"header" => ["header"],
			"footer" => ["footer"],
		       );
#
#
# Subroutine: xrt()
#
# Description: this is the main function you will be calling from
#              our scripts. please see 
#              www.caida.org/Tools/Graph/ for a full description
#              and how-to of this subroutine
#

sub xrt3d {
    my $user_global_opts_ref = shift; 
    my $data_set_ref = shift;
	my $matrix_data_ref;
	my $data_filename;
    my (%global_opts);
    
    # variables to be written to the command file
    my ($plot_file, $x_axis, $y_axis, $z_axis, $x_step, $y_step);
    my ($x_min, $y_min, $x_ticks, $y_ticks, $header, $footer);
    my ($x_cnt, $y_cnt, $hdr_cnt, $ftr_cnt);
    my ($output_file);

    if (@_) {
	carp 'Too many arguments. Usage: xrt3d(\%options, \@data_set)';
	return 0;
    }

    _make_tmpdir("_Xrt3d_");

    # set paths for external programs
    if (not _set_xrtpaths("xrt3d")) {
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
		$output_file = $value;
		unless (defined $global_opts{"output type"}) {
		    carp "Must have an output type defined";
		    _cleanup_tmpdir();
		    return 0;
		  }
	  }
		
	  # If the file is PostScript ... what XRT makes is PostScript
	  if ($global_opts{"output type"} eq "ps") {
		$plot_file = _make_tmpfile("plot", "ps");
	  } 
	  # For all raster formats XRT starts out with 
	  # X-Windows XWD format.
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
		
		if ($key eq "z-axis title") {
		    if(defined($value)) {
			$z_axis = $value;
		    }
		}
		
		if ($key eq "x-min") {
		    if(defined($value)) {
			$x_min = $value;
		    }
		}
		
		if ($key eq "y-min") {
		    if(defined($value)) {
			$y_min = $value;
		    }
		}
		
		if ($key eq "x-step") {
		    if(defined($value)) {
			$x_step = $value;
		    }
		}
		
		if ($key eq "y-step") {
		    if(defined($value)) {
			$y_step = $value;
	    }
		}
		
		if ($key eq "x-ticks") {
		    if(defined($value)) {
			$x_ticks = $value;
		    }
		}
		
		if ($key eq "y-ticks") {
		    if(defined($value)) {
			$y_ticks = $value;
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
	}
	

	# Extract options for data.
	my $data_opts = shift @{$data_set_ref};
	while (my ($key, $value) = each %{$data_opts}) {
	  if ($key eq "type") {
		if ($value eq "matrix") {
		  $matrix_data_ref = $data_set_ref;
		} 
		elsif ($value eq "file") {
		  $data_filename = pop @{$data_set_ref};
		} else {
		  carp "Unsupported or unknown format for data";
		}
	  }
	}
    
    # because xrt allows multiline headers
    # get the length of the header array
    # each line of the header is one index
    # in the array

    $hdr_cnt = $#{$global_opts{"header"}} + 1;
    $ftr_cnt = $#{$global_opts{"footer"}} + 1;
    

	if (defined($matrix_data_ref)) {
	  # get the number of columns and number of rows
	  # and verify that each row has same number of 
	  # columns
    
	  $x_cnt = $#{$matrix_data_ref} + 1;
	  my $tmp = $#{$matrix_data_ref->[0]} + 1;
    
	  foreach my $i (@{$matrix_data_ref}) {
		if ($tmp != $#{$i} + 1) {
		  carp "each row must have the same number of columns";
		  _cleanup_tmpdir();
		  return 0;
		} 
	  }
    
	  $y_cnt = $tmp;
	  

	  # verify that number of tick marks == corresponds
	  # to that of xy matrix. One tick mark for every x
	  # y.
	  
	  if (not _verify_ticks($x_cnt, $global_opts{"x-ticks"})) {
		_cleanup_tmpdir();
		return 0;
	  }
    
	  if (not _verify_ticks($y_cnt, $global_opts{"y-ticks"})) {
		_cleanup_tmpdir();
		return 0;
	  }
    } else {
	  # XXX
	  # Poor man's hack to compute rows and columns in data file.  This will
	  # make a second pass through file, but is probably faster than doing it
	  # in Perl.  
	  my ($lead, $words, $bytes);
	  ($lead, $x_cnt, $words, $bytes)  = split(/\D+/, `wc $data_filename`);
	  if (($x_cnt > 0) and ($words > 0)) {
		$x_cnt++;
		$y_cnt = $words/$x_cnt;
	  } else {
		$x_cnt = 0;
		$y_cnt = 0;
		carp "Cannot compute number of rows and/or columns in file data";
	  }
	}

    ##
    ## print command file using this format
    ##
    
    # output.file						
    # x_min (normally 0)
    # y_min (normally 0)
    # x_step (normally 1)
    # y_step (normally 1)
    # x_cnt (number of rows of input)
    # y_cnt (number of columns of input)
    # data11 data12 data13 data14 .... (x by y matrix of doubles)
    # data21 data22 data23 ....
    # .
    # .
    # .
    # datax1 datax2 ... dataxy
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
    # Title of z-axis
    # xlabel0 (x_cnt number of labels for ticks along x-axis)
    # ...
    # xlabelx
    # ylabel0 (y_cnt number of labels for ticks along y-axis)
    # ...
    # ylabely
    
    # create command file and open file handle 
    my $command_file = _make_tmpfile("command");
    my $handle = new FileHandle;
    if (not $handle->open(">$command_file")) {
	carp "could not open $command_file";
	_cleanup_tmpdir();
	return 0;
    }
    
    print $handle "$plot_file\n";
    print $handle "$x_min\n";
    print $handle "$y_min\n";
    print $handle "$x_step\n";
    print $handle "$y_step\n";
    print $handle "$x_cnt\n";
    print $handle "$y_cnt\n";
	if (defined($matrix_data_ref)) {
	  _print_matrix($handle, @{$matrix_data_ref});
	} else {
	   _transfer_file($handle, $data_filename)
	}
    print $handle "$hdr_cnt\n";
    _print_array($handle, @{$header});
    print $handle "$ftr_cnt\n";
    _print_array($handle, @{$footer});	
    print $handle "$x_axis\n";
    print $handle "$y_axis\n";
    print $handle "$z_axis\n";
    _print_array($handle, @{$x_ticks});	
    _print_array($handle, @{$y_ticks});
    $handle->close();

    # call xrt and convert file to gif
    if (not _exec_xrt3d($command_file)) {
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
 
1;


__END__

=head1 NAME

Chart::Graph::Xrt3d

=head1 SYNOPSIS

 #Include module
 use Chart::Graph::Xrt3d qw(xrt3d);

 # Function call
 xrt3d(\%options,
       \@data_set
      );


=head1 DESCRIPTION

This module is unmaintained, it worked with Sitraka's XRT, and hasn't been
tested against newer versions.

Sitraka (now Quest) makes a number of graphics packages for UNIX systems.  XRT is
a Motif-based commercial software product that has been adapted by
CAIDA using a combination of C drivers and Perl function I<xrt3d()>.
The Perl function I<xrt3d()> provides access to the three dimensional
graphing capabilities of XRT from Perl.  To access the two dimensional
graphing using XRT, use I<xrt2d()> also supplied in the
I<Chart::Graph> package.

=head1 ARGUMENTS

The options to I<xrt3d()> are listed below.  Additional control over the
resulting graph is possible by using the XRT application itself once
the graph has been created.

 +--------------------------------------------------------------------------+
 |                                OPTIONS                                   |
 +----------------+--------------------------+------------------------------+
 | Name           |  Options                 | Default                      |
 |"output file"   |  (set your own)          | "untitled-xrt3d.gif"         |
 |"output type"   |  "ps","xwd", "png", "jpg"| "xwd"                        |
 |"x-axis title"  |  (set your own)          | "x-axis"                     |
 |"y-axis title"  |  (set your own)          | "y-axis"                     |
 |"z-axis title"  |  (set your own)          | "z-axis"                     |
 |"x-min"         |  "0" or "1"(normally 0)  | "0"                          |
 |"y-min"         |  "0" or "1"(normally 0)  | "0"                          |
 |"x-step"        |  "0" or "1"(normally 1)  | "1"                          |
 |"y-step"        |  "0" or "1"(normally 1)  | "1"                          |
 |"x-ticks"       |  (set your own)          | none                         |
 |"y-ticks"       |  (set your own)          | none                         |
 |"header"        |  (set your own)          | Array ref of "header" text   |
 |"footer"        |  (set your own)          | Array ref of "footer" text   |
 +----------------+--------------------------+------------------------------+

The I<xrt3d> function only accepts data in one of two forms.  The
choices are: either C<[\%data1_opts, \@data_matrix]> or
C<[\%data1_opts, "filename"]> The data options are listed below.

 +--------------------------------------------------------------------------+
 |                             DATA OPTIONS                                 |
 +----------------+--------------------------+------------------------------+
 | Name           |  Options                 | Default                      |
 +----------------+--------------------------+------------------------------+
 | "type"         | Data format: "matrix" or | none                         |
 |                | "file"                   |                              |
 +----------------+--------------------------+------------------------------+

=head2 DETAILS ON GRAPHICS CONVERTER OPTIONS

The xrt package supports only two graphics formats internally:
Postscript and the X windows format XWD.  Additional raster graphics
formats are supported with Chart::Graph by using one of two graphics
converter packages: I<Imagemagick> and I<Netpbm>.

If you need to install a converter package,I<Imagemagick> 
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

The following four examples show Chart::Graph::Xrt3d in different roles
and producing different styles of output.

=head2 EXAMPLE: STOCK PRICES FOR JOE'S RESTAURANT

The first example creates a three dimensional bar chart of
fictitious stock data that is displayed in the graphic file
F<xrt3d-1.gif>.  Note that I<xrt3d()> uses the older gif file format,
but can use others as noted above if you have the available converters
provided.

 #make sure to include Chart::Graph
 use Chart::Graph::Xrt3d qw(xrt3d);

 #using a 3 by 6 matrix for the data set
 xrt3d({"output file" => "xrt3d-1.gif",
			   "output type" => "gif",
				   "header" => 
			   ["Stock prices for Joe's restaurant chain",
				"Compiled from local records"
				],
			   "footer" =>
			   ["Joe's Restaurant"],
			   "y-ticks"=>["Jan/Feb", "Mar/Apr", "May/Jun", "Jul/Aug",
						   "Sep/Oct", "Nov/Dec"],
			   "x-axis title" => "Years monitored",
			   "y-axis title" => "Month's tracked",
			   "z-axis title" => "Stock prices",
			  },
			  [{"type" => "matrix"},
               ["4", "5", "3", "6", "6", "5"],
			   ["8", "13", "20", "45", "100", "110" ],
			   ["70", "45", "10", "5", "4", "3"]])

=for html
<p><center><img src="http://www.caida.org/tools/utilities/graphing/xrt3d-1.jpg"></center></p>
<p><center><em>xrt3d-1.jpg</em></center></p>

=head2 EXAMPLE: EARLY GROWTH OF THE INTERNET

The following example creates a three dimensional bar chart of data
collected on the early growth of the Internet (URL and corporate
source included on graph.)  The result in this case is display in one
of the newest graphics formats the PNG format: F<xrt3d-2.png>.

 #make sure to include Chart::Graph
 use Chart::Graph::Xrt3d qw(xrt3d);

  xrt3d({"output file" => "xrt3d-2.png",
         "output type" => "png",
	 "header" => 
	 ["Growth of Early Internet", 
	 "(according to Internet Wizards - http://www.nw.com/)",
	 ],
	 "footer" =>
	 ["http://www.mit.edu/people/mkgray/net/internet-growth-raw-data.html"],
	 "y-ticks"=>["Jan 93", "Apr 93", "Jul 93",
		     "Oct 93", "Jan 94", "Jul 94",
		     "Oct 94", "Jan 95", "Jul 95",
		     "Jan 96"
		    ],
	 "x-ticks"=>["Hosts", "Domains", "Replied to Ping"],},
	 [{"type" => "matrix"},
      ["1.3e6", "1.5e6", "1.8e6", "2.1e6", "2.2e6", "3.2e6", 
	   "3.9e6","4.9e6", "6.6e6", "9.5e6"
	  ],
	  ["21000","22000", "26000", "28000", "30000", "46000", 
	   "56000", "71000", "120000", "240000"
	  ],
	  ["NA", "0.4e6", "NA", "0.5e6", "0.6e6", "0.7e6", 
	   "1.0e6", "1.0e6", "1.1e6", "1.7e6" 
	  ]
	 ]
        );

=for html
<p><center><img src="http://www.caida.org/tools/utilities/graphing/xrt3d-2.png"></center></p>
<p><center><em>xrt3d-2.png</em></center></p>

=head2 EXAMPLE: USING A DATA FILE FOR INPUT

The next example uses a file instead of a array for it's data source.
The file is listed below the Perl code.

 #make sure to include Chart::Graph
 use Chart::Graph::Xrt3d qw(xrt3d);

    if (xrt3d({"output file" => "xrt3d-3.gif",
               "output type" => "gif",
               "x-ticks"=>["a", "b", "c"],
	       "y-ticks"=>["w", "x", "y", "z"],},
	    [{"type" => "file"},
		 "xrt3d_data.txt"])) {
	print "ok\n";
    } else {
	print "not ok\n";
    }

The data file used in the above example is as follows.

 10 15 23  10
 4 13 35 45
 29 15 64 24

=for html
<p><center><img src="http://www.caida.org/tools/utilities/graphing/xrt3d-3.gif"></center></p>
<p><center><em>xrt3d-3.gif</em></center></p>


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

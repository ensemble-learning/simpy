## Gnuplot.pm is a sub-module of Graph.pm. It has all the subroutines 
## needed for the gnuplot part of the package.
##
## $Id: Gnuplot.pm,v 1.48 2006/06/07 21:09:33 emile Exp $ $Name:  $
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
package Chart::Graph::Gnuplot;
use Exporter ();

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(&gnuplot);

use Carp;			# for carp() and croak()
use Chart::Graph::Utils qw(:UTILS);	# get global subs and variable
use POSIX 'strftime';
use FileHandle;

$cvs_Id = '$Id: Gnuplot.pm,v 1.48 2006/06/07 21:09:33 emile Exp $';
$cvs_Author = '$Author: emile $';
$cvs_Name = '$Name:  $';
$cvs_Revision = '$Revision: 1.48 $';

$VERSION = 3.2;

use strict;

use vars qw($show_year $show_seconds);

# these variables  hold default options for gnuplot
my %def_gnu_global_opts = (
			   "title" => "untitled",
			   "output type" => "png",
			   "output file" => "untitled-gnuplot.png",
			   "x-axis label" => "x-axis",
			   "y-axis label" => "y-axis",
			   "x2-axis label" => undef,
			   "y2-axis label" => undef,
			   "logscale x" => "0",
         		   "logscale y" => "0",
			   "logscale x2" => "0",
			   "logscale y2" => "0",
			   "xtics" => undef,
			   "ytics" => undef,
			   "x2tics" => undef,
			   "y2tics" => undef,
			   "xdata" => undef,
			   "ydata" => undef,
			   "x2data" => undef,
			   "y2data" => undef,
			   "timefmt" => undef,
			   "format" => undef,
			   "xrange" => undef,
			   "yrange" => undef,
                           "extra_opts" => undef,
			   "uts" => undef,
			   "uts_normalize" => undef,
			   "size"=> undef,
			  );


my %def_gnu_data_opts = (
			 "title" => "untitled data",
			 "style" => "points", # points, lines...
			 "axes" => "x1y1", 
			 "type" => undef,
                         "using" => "1:2",
			);



#
#
# Subroutine: gnuplot()
#
# Description: this is the main function you will be calling from
#              our scripts. please see 
#              www.caida.org/Tools/Graph/ for a full description
#              and how-to of this subroutine
#

sub gnuplot {
    my ($user_global_opts_ref, @data_sets) = @_;
    my (%data_opts, %global_opts,);
    my ($plottype, $output_file, $plot_file, $output_type, $data_set_ref);
    # create a new filehandle to be used throughout package
    my $handle = new FileHandle;

    # create tmpdir
    _make_tmpdir("_Gnuplot_");

    # set paths for external programs
    if (not _set_gnupaths()) {
	_cleanup_tmpdir();
	return 0;
    }
    
    # check first arg for hash
    if (ref($user_global_opts_ref) ne "HASH") {
	carp "Global options must be a hash";
	_cleanup_tmpdir();
	return 0;
    }

    # check for data sets
    if (not @data_sets) {
	carp "no data sets";
	$handle->close;
	_cleanup_tmpdir();
	return 0;
    }
    
    # call to combine user options with default options
    %global_opts = _mesh_opts($user_global_opts_ref, \%def_gnu_global_opts);
    
    my $command_file = _make_tmpfile("command");

    #remember to close the file if we return
    if (not $handle->open(">$command_file")) {
	carp "could not open file: $command_file";
	_cleanup_tmpdir();
	return 0;
	}

    # check if uts option is chosen and process first
    if (my $value = $global_opts{uts}) {
	if (defined($value) and ref($value) eq "ARRAY") {
	    if (@{$value} < 2 || @{$value} > 4) {
		carp "out of range for 'uts': [start, end, <scale>, <use_local_timezone>]\n";
		_cleanup_tmpdir();
		return 0;
	    }		
	    # set x tics to human readable time stamps
	    _gnuplot_date_utc($value->[0], $value->[1], $value->[2], $value->[3], \%global_opts);
	} else {
	    carp "Invalid value for 'uts', give [start, end, <scale>, <use_local_timezone>]\n";
	}
    }

    # uts_normalize will be removed in a future version, so don't use it
    if (my $value = $global_opts{uts_normalize}) {
	if (defined($value) and ref($value) eq "ARRAY") {
	    if (@{$value} < 2 || @{$value} > 3) {
		carp "out of range for 'uts_normalize': [start, end, <scale> ]\n";
		_cleanup_tmpdir();
		return 0;
	    }		
	    # set x tics to human readable time stamps
	    _gnuplot_date_utc_normalize($value->[0], $value->[1], $value->[2], \%global_opts);
	} else {
	    carp "Invalid value for 'uts_normalize', give [start, end]";
	}
    }

    # Check if we have options for reading data as date/time formats
    # these must be in command file before any others.
    foreach my $time_set ("xdata", "ydata", "x2data", "y2data") {
	if (defined($global_opts{$time_set})) {
	    print $handle "set $time_set time\n";
	}
    }
	
    # Set the format for reading data/time data. Only one format for all axes.
    if (my $value = $global_opts{timefmt}) {
	if (defined($value)) {
	    print $handle "set timefmt \"$value\"\n";
	}
    }

    # Now write remain global options to command file
    while (my ($key, $value) = each %global_opts) {
	
	## Generic pass-thru, stuff random Gnuplot commands in this key
        if ($key eq "extra_opts") {
	    if (defined $value) {
		if (ref($value) eq 'ARRAY') {
		    # arrayref
		    print $handle join("\n", @$value);
		    print $handle "\n"; 
		} else {
		    # assume it's a string and user provided \n's
		    print $handle "$value\n";
		}
	    }
        }

	if ($key eq "title") {
	    print $handle "set title \"$value\"\n";
	}
	
	if ($key eq "x-axis label") {
	    print $handle "set xlabel \"$value\"\n";
	}
	
	if ($key eq "y-axis label") {
	    print $handle "set ylabel \"$value\"\n";
	}
	
	if ($key eq "x2-axis label") {
	    if (defined($value)) {
		print $handle "set x2label \"$value\"\n";
	    }
	}
	
	if ($key eq "y2-axis label") {
	    if (defined($value)) {
		print $handle "set y2label \"$value\"\n";
	    }
	}
	
	if ($key eq "logscale x") {
	    if ($value == 1) {
		print $handle "set logscale x\n";
	    }
	}
	
	if ($key eq "logscale y") {
	    if ($value == 1) {
		print $handle "set logscale y\n";
			}
	}
	if ($key eq "logscale x2") {
	    if ($value == 1) {
		print $handle "set logscale x2\n";
	    }
	}
	
	if ($key eq "logscale y2") {
	    if ($value == 1) {
		print $handle "set logscale y2\n";
	    }
	}
	
	# tics are not required so we can fall through if we want 
	if ($key eq "xtics") {
	    if (defined($value) and ref($value) eq "ARRAY") {
		_print_tics($handle, $key, @{$value});
	    }
	}

	if ($key eq "ytics") {
	    if (defined($value) and ref($value) eq "ARRAY") {
		_print_tics($handle, $key, @{$value});
	    }
	}
	
	if ($key eq "x2tics") {
	    if (defined($value)) { 
		if (ref($value) eq "ARRAY") {
		    _print_tics($handle, $key, @{$value});
		}
		if ($value eq "on") {
		    print $handle "set $key\n";
		}
	    }		
	}
	if ($key eq "y2tics") {
	    if (defined($value)) { 
		if (ref($value) eq "ARRAY") {
		    _print_tics($handle, $key, @{$value});
		}
		if ($value eq "on") {
		    print $handle "set $key\n";
		}
	    }		
	}
	if ($key eq 'xrange' || $key eq 'yrange' ) {
	    if (defined($value)) { 
		if (ref($value) eq 'ARRAY') {
		    # arrayref
		    print $handle "set $key [$value->[0] : $value->[1]]\n";
		} else {
		    # assume string
		    print $handle "set $key $value\n";
		}
	    }
	}
	# The only time related Gnuplot code that doesn't need to be
	# output first.
	if ($key eq "format") {
	    if (defined($value)) {
		if(ref($value) eq "ARRAY") {
		    # Print value supplying quotes for time format.
		    print $handle "set $key ", $$value[0]," \"",
			$$value[1], "\" \n";
		} else {
		    carp "Invalid setting for format option";
		}
	    }
	}
	# Date/Time and UTS keys which are processed first.
	if ($key eq "timefmt") {
	    # already processed
	}
	if ($key eq "xdata") {
	    # already processed
	}
	if ($key eq "ydata") {
	    # already processed
	}
	if ($key eq "x2data") {
	    # already processed
	}
	if ($key eq "y2data") {
	    # already processed
	}
	if ($key eq "uts") {
	    # already been processed 	    
	}
	if ($key eq "uts_normalize") {
	    # already been processed 	    
	}
	if ($key eq "size" && defined $value) {
            if (ref($value) eq 'ARRAY' && @{$value} == 2) {
		print $handle "set $key ". @{$value}[0] . "," . @{$value}[1] . "\n";
	    }
	    else {
		
	       print STDERR "option `size' must be given a two element array\n";
	    }
	}
	    
	if ($key eq "output file") {
	    $output_file = $value;
	}
	
	if ($key eq "output type") {
	    if (!($value =~ /^(pbm|gif|tgif|png|svg|eps(:? .*)?)$/)) {
		carp "invalid output type: $value";
		$handle->close();
		_cleanup_tmpdir();
		return 0;
	    }
	    $output_type = $value;
	}
    }

    # create the data file
    if ($output_type =~ /^eps( .*)?$/) {
	my $options = $1 || "";
	if (defined $output_file) {
	    $plot_file = _make_tmpfile("plot", "eps");
	    print $handle "set output \"$plot_file\"\n";
	}
	#print $handle "set terminal postscript eps color \"Arial\" 18\n";
	print $handle "set terminal postscript eps $options\n";
    } elsif ($output_type eq "pbm" ) {
	if (defined $output_file) {
	    $plot_file = _make_tmpfile("plot", "pbm");
	    print $handle "set output \"$plot_file\"\n";
	}
	print $handle "set terminal pbm small color\n";
    } elsif ($output_type eq "gif") {
	# always needs the tempfile because of conversion later on
	$plot_file = _make_tmpfile("plot", "pbm");
	print $handle "set output \"$plot_file\"\n";
	print $handle "set terminal pbm small color\n";
    } elsif ($output_type eq "png") {
	if (defined $output_file) {
	    $plot_file = _make_tmpfile("plot", "png");
	    print $handle "set output \"$plot_file\"\n";
	}
	print $handle "set terminal png small\n";
    } elsif ($output_type eq "tgif") {
	if (defined $output_file) {
	    $plot_file = _make_tmpfile("plot", "obj");
	    print $handle "set output \"$plot_file\"\n";
	}
	print $handle "set terminal tgif\n";
    } elsif ($output_type eq 'svg') {
	if (defined $output_file) {
	    $plot_file = _make_tmpfile("plot", "svg");
	    print $handle "set output \"$plot_file\"\n";
	}
	print $handle "set terminal svg\n";
    }

    # process data sets
    print $handle "plot ";
    while (@data_sets) {
	
	$data_set_ref = shift @data_sets;
	
	if (ref($data_set_ref) ne "ARRAY") {
	    carp "Data set must be an array";
	    $handle->close();
	    _cleanup_tmpdir();
	    return 0;
	}
	
	if (not _gnuplot_data_set($handle, @{$data_set_ref})) {
	    ## already printed error message
	    $handle->close();
	    _cleanup_tmpdir();
	    return 0;
	}
	
	if (@data_sets) {
	    print $handle ", ";
	}
    }	
    
    $handle->close();
	
    # gnuplot and convert pbm file to gif
    if (not _exec_gnuplot($command_file)) {
	_cleanup_tmpdir();
	return 0;
    }
    
    if ($output_type eq "gif") {
	if(not _exec_pbmtogif($plot_file, $output_file)) {
	    _cleanup_tmpdir();
	    return 0;
	}
    } elsif (defined $output_file && $output_type =~ /^(pbm|eps(?: .*)?|png|tgif)$/) {
	#try to get rid of the ugly warnings when moving a file on freebsd
	if ($^O eq 'freebsd') {
	    eval { #this is opportunistic, if it fails it doesn't really matter
		my $gid = `/usr/bin/id -g`;
		chomp($gid);
		system('/usr/bin/chgrp',$gid,$plot_file);
	    };
	}

	my $status = system("mv", "$plot_file", "$output_file");

        if (not _chk_status($status)) {
	    if ($Chart::Graph::debug) {
		print STDERR "Couldn't mv $plot_file to $output_file: $!\n";
	    }
	    _cleanup_tmpdir();
	    return 0;
	}
    } 
    _cleanup_tmpdir();
    return 1;
}

#
#
# Subroutine: gnuplot_data_set()
# 
# Description: this functions processes the X number
#              of data sets that a user gives as 
#              arguments to gnuplot(). Again, please
#              see http://www.caida.org/Tools/Graph/
#              for the format of the dataset.
#


 
sub _gnuplot_data_set {
    my ($handle, $user_data_opts_ref, @data) = @_;
    my (%data_opts);
    
    # set these values with empty string because we print them out later 
    # we don't want perl to complain of uninitialized value. 
    my ($title, $style, $axes, $ranges, $type,) = ("", "", "", "", "");
    my ($using) = ("");
    my $result;
    my $filename = _make_tmpfile("data");
    
    ## check first arg for hash
    if (ref($user_data_opts_ref) ne "HASH") {
	carp "Data options must be a hash.";
	return 0;
    }

    # call to combine user options with default options   
    %data_opts = _mesh_opts($user_data_opts_ref, \%def_gnu_data_opts);

    # write data options to command file
    while (my ($key, $value) = each %data_opts) {
	
        if ($key eq "using") {
            $using = "using $value";
        }

	if ($key eq "title") {
	    $title = "title \"$value\"";
	}
	
        if ($key eq "style") {
	    $style = "with $value"
	}
		
	if ($key eq "axes") {
	    $axes = "axes $value";
	}
	
	if ($key eq "type") {
	    $type = $value;
	}
    }
    
    if ($type eq "function") {
        #$ranges = "[t=:]";	# XXX ?
        print $handle "$ranges " . $data[0] . " $axes $title $style";
	return 1;	
    } else {
        print $handle "$ranges \"$filename\" $using $axes $title $style";

	# we give the user 3 formats for supplying the data set
	# 1) matrix
	# 2) column
	# 3) file
	# please see the online docs for a description of these 
	# formats
	if ($type eq "matrix") {
	    $result = _matrix_to_file($filename, @data);
	} elsif ($type eq "columns") {
	    $result = _columns_to_file($filename, @data);
	} elsif ($type eq "file") {
	  $result = _file_to_file($filename, @data);
	} elsif ($type eq "") {
	    carp "Need to specify data set type";
	    return 0;
	} else {
	    carp "Illegal data set type: $type"; 
	    return 0;
	}
    }
    return $result;
}

# 
# Subroutine: set_gnupaths()
# 
# Description: set paths for external programs required by gnuplot()
#              if they are not defined already
#

sub _set_gnupaths {

    if (not defined($gnuplot)) {
	if (not $gnuplot = _get_path("gnuplot")) {
	    return 0;
	}
    }
   
    if (not defined($ppmtogif)) {
	if (not $ppmtogif = _get_path("ppmtogif")) {
	    return 0;
	}
    }
    return 1;
}

#
#
#  Subroutine: print_tics()
#  Description: this subroutine takes an array 
#               of graph tic labels and prints 
#               them to the gnuplot command file.
#               This subroutine is called by gnuplot().
#
#  Arguments: $tic_type: which axis to print the tics on
#             @tics: the array of tics to print to the file
#

sub _print_tics {
    my ($handle, $tic_type, @tics) = @_;
    my (@tic_array, $tics_formatted, $tic_label, $tic_index);
   
    # no tic set found, user entered empty tic array
    if (not @tics) {
	carp "Warning: empty tic set found";
	return 1;
    }

    foreach my $tic (@tics) {
	#tics can come in two formats
	#this one is [["label1", 10], ["label2", 20],...]
	if (ref($tic) eq "ARRAY") {
	    
	    if ($#{$tic} != 1) {
		carp "invalid tic format";
		return 0;
	    }
	    $tic_label = $tic->[0];
	    $tic_index = $tic->[1];
	    push (@tic_array, "\"$tic_label\" $tic_index");
	    # this one is [10, 20,...]
	} else {
	    push (@tic_array, "$tic");
	}
    }
    $tics_formatted = join(",", @tic_array);
    print $handle "set $tic_type ($tics_formatted)\n";
    return 1;
}

#
#
# Subroutine: matrix_to_file()
#
# Description: converts the matrix data input into a the gnuplot
#              data file format. See www for the specific on the 
#              matrix format
# 
# 
sub _matrix_to_file {
    my ($file, $matrix_ref) = @_;
    my $entry_ref;
    my $matrix_len;
    
    if (ref($matrix_ref) ne "ARRAY") {
	carp "Matrix data must be a reference to an array";
	return 0;
    }
    
    open (DATA, ">$file");
    
    $matrix_len = @{$matrix_ref};
    for (my $i = 0; $i < $matrix_len; $i++) {
	$entry_ref = $matrix_ref->[$i];
	
	if (ref($entry_ref) ne "ARRAY") {
	    carp "Matrix entry must be a reference to an array";
	    close DATA;
	    return 0;
	}
	
        # prints blank lines for blank entries, this allows 
        # the user to tell gnuplot to not connect lines between 
        # all points when displaying data with lines.
        if (@{$entry_ref} == 0)
        {
            print DATA "\n";
        }
        else
        {
if (0) {
            # check that each entry ONLY has two entries
            if (@{$entry_ref} != 2) {
                carp "Each entry must be an array of size 2";
                return 0;
            }
            print DATA $entry_ref->[0], "\t", $entry_ref->[1], "\n";
}
	    # XXX
            print DATA join("\t", @{$entry_ref}), "\n";
        }
    }
    close DATA;
    return 1;
}

#
#
# Subroutine: columns_to_file()
#
# Description: converts the column data input into a the gnuplot
#              data file format. please see www page for specifics
#              on this format.
#

sub _columns_to_file {
    my ($file, @columns) = @_;

    foreach my $dataset ( @columns ) {
        if (!(ref($dataset) eq "ARRAY")) {
            carp "Column data must be a reference to an array";
            return 0;
        }
   
        if ($#{$dataset} != $#{$columns[$[]}) {
            carp "All columns must be of same length";
            return 0;
        }
    }

    if ($#{$columns[$[]} == 0) {
        carp "Warning: Columns have no data!";
    }

    open (DATA, ">$file");
   
    for (my $i = 0; $i <= $#{$columns[$[]}; $i++) {
        foreach my $dataset ( @columns ) {
            print DATA "$dataset->[$i]\t";
        }
        print DATA "\n";
    }

    close DATA;
    return 1;
}


#
# Subroutine: file_to_file()
#
# Description: If a gnuplot data set was given in
#              file format, we simply copy the data 
#              and read it into 
#

sub _file_to_file {
    my ($file_out, $file_in) = @_;
    
    if (not $file_in) {
	carp "Data set file missing";
	return 0;
    }

    if (not -f $file_in) {
	carp "Data set file, '$file_in', does not exist.";
	return 0;
    }
    
    my $status = system("cp", "$file_in", "$file_out");
    
    if (not _chk_status($status)) {
	return 0;
    }
    
    return 1;
}

# 
# Subroutine: exec_gnuplot()
#
# Description: this executes gnuplot on the command file 
#              and data sets that we have generated.
#

sub _exec_gnuplot {
    my ($command_file) = @_;
    my $status = system("$gnuplot", "$command_file");
    
    if (not _chk_status($status)) {
	return 0;
    }
    
    return 1;
}	

#
# Subroutine: exec_pbmtogif()
# 
# Description: convert pbm file that gnuplot makes into
#              a gif. usually used for web pages
#
sub _exec_pbmtogif {
    my ($pbm_file, $gif_file) = @_;
    my $status;
    my $cmd = "$ppmtogif $pbm_file ";
    if ($gif_file) {
	$cmd .= "> $gif_file ";
    } 
    unless ($Chart::Graph::debug) {
	$cmd .= "2> /dev/null ";
    }
    $status = system($cmd);
    
    if (not _chk_status($status)) {
	return 0;
    }
    
    return 1;
}

#
# Subroutine: gnuplot_date_utc()
#
# Description: wrapper function that handles UNIX
#              time stamps as x values nicely
#
# Author: Ryan Koga - rkoga@caida.org
#

sub _gnuplot_date_utc {
    my ($start, $end, $samp_scale, $use_local_tz, $global_options) = @_;

    my $min_len = 60;
    my $hour_len = $min_len*60;
    my $day_len = $hour_len*24;
    my $interval = $end - $start;
    my $min_samp;
    my @tics;

    if (!defined($samp_scale)) {
        $samp_scale = 1;
    }

    if ($interval < 10) {
	$min_samp = 1;
    } elsif ($interval < 30) {
	$min_samp = 4;
    } elsif ($interval < $min_len) {
	$min_samp = 10;
    } elsif ($interval < 3*$min_len) {
	$min_samp = 30;
    } elsif ($interval < 10*$min_len) {
	$min_samp = $min_len;
    } elsif ($interval < $hour_len) {
	$min_samp = 5*$min_len;
    } elsif ($interval < 2*$hour_len) {
	$min_samp = 10*$min_len;
    } elsif ($interval < 3*$hour_len) {
	$min_samp = 15*$min_len;
    } elsif ($interval < 4*$hour_len) {
	$min_samp = 20*$min_len;
    } elsif ($interval < 5*$hour_len) {
	$min_samp = 30*$min_len;
    } elsif ($interval < 12*$hour_len) {
	$min_samp = $hour_len;
    } elsif ($interval < $day_len) {
	$min_samp = 2*$hour_len;
    } elsif ($interval < 2*$day_len) {
	$min_samp = 4*$hour_len;
    } elsif ($interval < 5*$day_len) {
	$min_samp = 12*$hour_len;
    } elsif ($interval < 7*$day_len) {
	$min_samp = $day_len;
    } elsif ($interval < 15*$day_len) {
	$min_samp = 2*$day_len;
    } elsif ($interval < 30*$day_len) {
	$min_samp = 7*$day_len;
    } elsif ($interval < 365*$day_len) {
	$min_samp = 30*$day_len;
    } elsif ($interval < 2*365*$day_len) {
	$min_samp = 60*$day_len;
    } else {
	$min_samp = 120*$day_len;
    }
    $min_samp /= $samp_scale;
    my $start_min = int($start/$min_samp);
    my $end_min = int($end/$min_samp);

    for (my $curr_min = $start_min; $curr_min <= $end_min; $curr_min++) {
	my $bucket = $curr_min*$min_samp;
	my ($bucket_str,@time_data);
	if ( $use_local_tz ) {
	    @time_data = localtime($bucket);
	} else {
	    @time_data = gmtime($bucket);
	}
	$time_data[$#time_data] = -1;  # unset dst data, broken strftime

	# keep compatibility with the undocumented 'utc_seconds' global var
	$Chart::Graph::Gnuplot::show_seconds = $Chart::Graph::Gnuplot::utc_seconds;

	if ($min_samp >= $min_len && !$Chart::Graph::Gnuplot::show_seconds) {
	    $bucket_str = strftime("%H:%M", @time_data);
	} else {
	    $bucket_str = strftime("%H:%M:%S", @time_data);
	}
	if ($bucket_str =~ /^00:00(:00)?$/ || $curr_min == $start_min + 1) {
	    $bucket_str .= strftime('\n%m/%d', @time_data);
	    if ($Chart::Graph::Gnuplot::show_year) {
		$bucket_str .= strftime('\n%Y', @time_data);
	    }
	}
	push @tics, [ $bucket_str, $bucket ];
    }

    # must check to see if xtics were previously set in the globals
    # if they are, we'll append them to the time stamp tics
    # note: collisions are handled by the user

    if (defined ($global_options->{"xtics"})) {
	push @tics, @{$global_options->{"xtics"}};
    }

    $global_options->{"xtics"} = \@tics;
    return 1;
}


sub _gnuplot_date_utc_normalize {
    my ($start, $end, $samp_scale, $global_options) = @_;
    ### this code used to be used as a workaround for an old gnuplot bug
    ### newer versions don't need it
    carp "'uts_normalize' is going to be depreciated in a future release\n";

    my $min_len = 60;
    my $hour_len = $min_len*60;
    my $day_len = $hour_len*24;
    my $interval = $end - $start;
    my $min_samp;

    if (!defined($samp_scale)) {
        $samp_scale = 1;
    }

    if ($interval < 10) {
	$min_samp = 1;
    } elsif ($interval < 30) {
	$min_samp = 4;
    } elsif ($interval < $min_len) {
	$min_samp = 10;
    } elsif ($interval < 3*$min_len) {
	$min_samp = 30;
    } elsif ($interval < 10*$min_len) {
	$min_samp = $min_len;
    } elsif ($interval < $hour_len) {
	$min_samp = 5*$min_len;
    } elsif ($interval < 2*$hour_len) {
	$min_samp = 10*$min_len;
    } elsif ($interval < 3*$hour_len) {
	$min_samp = 15*$min_len;
    } elsif ($interval < 4*$hour_len) {
	$min_samp = 20*$min_len;
    } elsif ($interval < 5*$hour_len) {
	$min_samp = 30*$min_len;
    } elsif ($interval < 12*$hour_len) {
	$min_samp = $hour_len;
    } elsif ($interval < $day_len) {
	$min_samp = 2*$hour_len;
    } elsif ($interval < 2*$day_len) {
	$min_samp = 4*$hour_len;
    } elsif ($interval < 5*$day_len) {
	$min_samp = 12*$hour_len;
    } elsif ($interval < 7*$day_len) {
	$min_samp = $day_len;
    } elsif ($interval < 15*$day_len) {
	$min_samp = 2*$day_len;
    } elsif ($interval < 30*$day_len) {
	$min_samp = 3*$day_len;
    } else {
	$min_samp = 30*$day_len;
    }
    $min_samp /= $samp_scale;
    my $start_min = int($start/$min_samp);
    my $end_min = int($end/$min_samp);
    my @tics;

    my $first_date_shown = 0;
    for (my $curr_min = $start_min; $curr_min <= $end_min; $curr_min++) {
	my $bucket = $curr_min*$min_samp;
	my $bucket_str;
	my @time_data = gmtime($bucket);
	$time_data[$#time_data] = -1;  # unset dst data, broken strftime
	if ($min_samp >= $min_len) {
	    $bucket_str = strftime('%H:%M', @time_data);
	} else {
	    $bucket_str = strftime('%H:%M:%S', @time_data);
	}
	my $show_date = 0;
	if ($bucket_str =~ /^00:00(:00)?$/) {
	    $show_date = 1;
	    $first_date_shown = 1;
	}
	if ($curr_min == $start_min + 1 && !$first_date_shown) {
	    $show_date = 1;
	}
	if ($show_date) {
	    $bucket_str .= strftime('\n%m/%d', @time_data);
	}
	push @tics, [ $bucket_str, ($bucket - $start) / $end ];
    }

    # must check to see if xtics were previously set in the globals
    # if they are, we'll append them to the time stamp tics
    # note: collisions are handled by the user

    if (defined ($global_options->{"xtics"})) {
	push @tics, @{$global_options->{"xtics"}};
    }
    
    $global_options->{"xtics"} = \@tics;    
    return 1;
}

1;

__END__

=head1 NAME

Chart::Graph::Gnuplot

=head1 SYNOPSIS

    use Chart::Graph::Gnuplot qw(&gnuplot);



 gnuplot(\%global_options, [\%data_set_options, \@matrix],
                           [\%data_set_options, \@x_column, \@y_column],
                           [\%data_set_options, < filename >], ... );

=head1 DESCRIPTION

I<gnuplot()> is a function in module Chart::Graph that lets you
generate graphs on the fly in perl. It was written as a front-end
application to gnuplot for hassle-free generation of
graphs. I<gnuplot()> can be supplied with many of the same options and
arguments that can be given to gnuplot. For more information on
I<gnuplot> see the end of this section.

=head1 OPTIONS

I<gnuplot()> has a very large number of options corresponding to
options available with the gnuplot application itself.  This Perl
wrapper provides a large subset of the functionality of the
application.

 +----------------------------------------------------------------------------+
 |                             GLOBAL OPTIONS:                                |
 +----------------+-----------------------------+-----------------------------+
 | NAME           |  OPTIONS                    |        DEFAULT              |
 +----------------+-----------------------------+-----------------------------+
 |'title'         |  set your own title         |     'untitled'              |
 |'output type'   |  'pbm','gif','tgif','png',  |     'png'                   |
 |                |   'svg' or "eps $epsoptions"|                             |
 |'output file'   |  set your own output file,  |     'untitled-gnuplot.png'  |
 |                |   undef to output to STDOUT |                             |
 |'x-axis label'  |  set your own label         |     'x-axis'                |
 |'y-axis label'  |  set your own label         |     'y-axis'                |
 |'x2-axis label' |  set your own label         |     none                    |
 |'y2-axis label' |  set your own label         |     none                    |
 |'logscale x'    |  0 or 1                     |     0                       |
 |'logscale y'    |  0 or 1                     |     0                       |
 |'logscale x2'   |  0 or 1                     |     0                       |
 |'logscale y2'   |  0 or 1                     |     0                       |
 | 'xtics'        | set your own tics on x-axis |     none                    |
 |                |   (see example below)       |                             |
 | 'x2tics'       | set your own tics on x2-axis|     none                    |
 |                |   (see example below)       |                             |
 | 'ytics'        | set your own tics on y-axis |     none                    |
 |                |   (see example below)       |                             |
 | 'y2tics'       | set your own tics on y2-axis|     none                    |
 |                |   (see example below)       |                             |
 | 'xrange'       | set xrange, accepts both    |     none                    |
 |                |  string '[$xmin:$xmax]'     |                             |
 |                |  or arrayref [$xmin,$xmax]  |                             |
 | 'yrange'       | set yrange, see xrange      |     none                    |
 |                |                             |                             |
 | 'uts'          | set your own range in unix  |     none                    |
 |                |  timestamps, array ref:     |                             |
 |                |  [start_ts,end_ts,<scale>,  |                             |
 |                |   <use_local_tz> ]          |                             |
 |                |  see UNIX TIMESTAMPS example|                             |
 | 'xdata'        | 'time' to indicate that     |     none                    |
 |                |  x-axis is date/time data   |                             |
 | 'ydata'        | 'time' to indicate that     |     none                    |
 |                |  y-axis is date/time data   |                             |
 | 'x2data'       | 'time' to indicate that     |     none                    |
 |                |  x2-axis is date/time data  |                             |
 | 'y2data'       | 'time' to indicate that     |     none                    |
 |                |  y2-axis is date/time data  |                             |
 | 'timefmt'      | "Input date/time string"    |     none                    |
 |                |  see Gnuplot manual for info|                             |
 | 'format'       | array ref: First element is |                             |
 |                |  axis: 'x', 'y', 'x2', 'y2'.|                             |
 |                |  Second element is          |                             |
 |                |  'output date/time string"  |                             |
 |                |  see Gnuplot manual for info|                             |
 | 'extra_opts'   | set your own Gnuplot        |     none                    |
 |                |  options, either an arrayref|                             |
 |                |  or string ("\n"-separated) |                             |
 | 'size'         | scale the display size of   |     none                    |
 |                |  the plot, arrayref [$x, $y]|                             |
 +----------------+-----------------------------+-----------------------------+

 +----------------------------------------------------------------------------+
 |                       Data Set Options:                                    |
 +----------------+-----------------------------+-----------------------------+
 |      Name      |          Options            |           Default           |
 +----------------+-----------------------------+-----------------------------+
 | 'type'         | 'matrix', 'columns', 'file',|      none                   |
 |                |  'function', see examples   |                             |
 |                |  below                      |                             |
 | 'title'        | set your own title          |     'untitled data'         |
 | 'style'        | 'points','lines','impulses' |     'points'                |
 |                |  'errorbars', etc...        |                             |
 |                |  see ERRORBARS example      |                             |
 | 'axes'         | 'x1y1', 'x2y2', 'x1y2', etc.|      'x1y1'                 |
 | 'using'        | map data to what will be    |      '1:2'                  |
 |                |  plotted, see ERRORBARS     |                             |
 |                |  example                    |                             |
 +----------------+-----------------------------+-----------------------------+

Data can be presented to I<Chart::Graph::Gnuplot> in one of 3 formats for
the convenience of the user:

 \@matrix: an array reference of [x,y] pairs of data

Alternatively:

 \@x_column, \@y_column: two array references of data of equal length.
 \@x_column is the x-axis data. \@y_column is the y-axis data.

Finally, data can be stored in a file.

=head2 USING GNUPLOT TO READ AND PLOT DATE/TIME DATA DIRECTLY

I<Gnuplot> now has the capability to read date/time data and to create
graphs which display date/time on any axis.  Unfortunately, mechanism
for reading data is less sophisticated than the mechanism for writing
data.  I<Chart::Graph::Gnuplot> implements date/time data in the same
way as I<Gnuplot> itself is presently implemented for consistency with
the application.

Any axis can be set to read date/time data instead of numerical
data. This is done by setting the options C<xdata>, C<ydata>,
C<x2data>, or C<y2data> to the value C<time>.  Unfortunately, you can
set only one format to read your data; therefore, consistency is
advised.  The input format is set using the C<timefmt> command noted
above.  The C<timefmt> command takes a string consisting of the
elements noted below.

I<Gnuplot> uses the same format codes for date/time input and output so
the following table applies to both situations.

 +---------+------------------------------------------------------------------+
 |  Format |                    Explanation                                   |
 +---------+------------------------------------------------------------------+
 |    %d   |       day of the month, 1--31                                    |
 |    %m   |       month of the year, 1--12                                   |
 |    %y   |       year, 0--99                                                |
 |    %Y   |       year, 4-digit                                              |
 |    %j   |       day of the year, 1--365                                    |
 |    %H   |       hour, 0--24                                                |
 |    %M   |       minute, 0--60                                              |
 |    %S   |       second, 0--60                                              |
 |    %b   |       three-character abbreviation of the name of the month      |
 |    %B   |       name of the month                                          |
 +---------+------------------------------------------------------------------+

In addition there are some additional special cases for reading
date/time data. To quote from I<Gnuplot> manual: "Any character is
allowed in the string, but must match exactly. C<\t> (tab) is
recognized. Backslash-octals (C<\nnn>) are converted to char. If there is
no separating character between the time/date elements, then %d, %m,
%y, %H, %M and %S read two digits each, %Y reads four digits and %j
reads three digits. %b requires three characters, and %B requires as
many as it needs."  I<Gnuplot> uses the space character as white space
pattern match - essentially the same as Perl's: C<\s*>.

I<Gnuplot> normally uses whitespace to separate datasets.  However,
I<Gnuplot> does recognize white space specified in the C<timefmt>
string.  So for example, x-y data can be specified in columns like
this:

 25/01/2001 03:05:39 2.05e-5

The C<timefmt> string required would be: C<"%d/%m/%y %H:%M:%S">.  Note
that while the month and month abbreviation can be accepted, any other
text must be matched (excluded) in the timefmt string.  Certainly,
representing dates as numerically is probably the most conservative.

Creating date/time labels for any of the axes is basically analogous.
The I<Chart::Graph:Gnuplot> global option is C<format>, and it takes a
two element array reference: the axis to be plotted and the format
string.  In addition to the time options listed above, C<format>
supports the following additional codes for formatting the numerical
data on the axes.

 +-------------+--------------------------------------------------------------+
 |  Format     |                Explanation                                   |
 +-------------+--------------------------------------------------------------+
 |    %f       |   floating point notation                                    |
 |    %e or %E |   exponential notation; an "e" or "E" before the power       |
 |    %g or %G |   the shorter of %e (or %E) and %f                           |
 |    %x or %X |   hex                                                        |
 |    %o or %O |   octal                                                      |
 |    %t       |   mantissa to base 10                                        |
 |    %l       |   mantissa to base of current logscale                       |
 |    %s       |   mantissa to base of current logscale; scientific power     |
 |    %T       |   power to base 10                                           |
 |    %L       |   power to base of current logscale                          |
 |    %S       |   scientific power                                           |
 |    %c       |   character replacement for scientific power                 |
 |    %P       |   multiple of pi                                             |
 +-------------+--------------------------------------------------------------+

As in the case of input there some additional options related to these
output formats.  Again to quote the I<Gnuplot> manual "Other
acceptable modifiers (which come after the I<%> but before the format
specifier) are I<->, which left-justifies the number; I<+>, which
forces all numbers to be explicitly signed; I<#>, which places a
decimal point after floats that have only zeroes following the decimal
point; a positive integer, which defines the field width; I<0> (the
digit, not the letter) immediately preceding the field width, which
indicates that leading zeroes are to be used instead of leading
blanks; and a decimal point followed by a non-negative integer, which
defines the precision (the minimum number of digits of an integer, or
the number of digits following the decimal point of a float)."

I<Gnuplot> also provides more flexibility in terms of the output format
codes available for date/time.  In addition to those shared with
input, the following codes can be used for formatting output date/time
axes only.

 +-------------+--------------------------------------------------------------+
 |  Format     |                Explanation                                   |
 +-------------+--------------------------------------------------------------+
 |    %a       |   abbreviated name of day of the week                        |
 |    %A       |   full name of day of the week                               |
 |    %b or %h |   abbreviated name of the month                              |
 |    %B       |   full name of the month                                     |
 |    %D       |   shorthand for "%m/%d/%y"                                   |
 |    %H or %k |   hour, 0--24                                                |
 |    %I or %l |   hour, 0--12                                                |
 |    %p       |   "am" or "pm"                                               |
 |    %r       |   shorthand for "%I:%M:%S %p"                                |
 |    %R       |   shorthand for %H:%M"                                       |
 |    %T       |   shorthand for "%H:%M:%S"                                   |
 |    %U       |   week of the year (week starts on Sunday)                   |
 |    %w       |   day of the week, 0--6 (Sunday = 0)                         |
 |    %W       |   week of the year (week starts on Monday)                   |
 +-------------+--------------------------------------------------------------+

Finally, I<Chart::Graph::Gnuplot> has an extension to support UNIX
timestamps.  Note this B<not> built into I<Gnuplot> itself.
Users can access this option by setting the C<xrange> using the C<uts> 
option instead.  UNIX timestamps are only available on the x-axis at this 
time.  They cannot be used on y, x2, or y2.  See the last example for more 
details on using UNIX timestamps.

=head1 EXAMPLES

The following are four examples on how to use I<Chart::Graph::Gnuplot> in
a variety of settings.

=head2 GENERAL EXAMPLE

The following example illustrates most of the general capabilities of
I<Chart::Graph::Gnuplot>. It creates the output file F<gnuplot1.png>.
in the I<png> graphics format.  The data is coming from all three
sources.  The first data source is a matrix, the second is a column,
and the last is an external data file.

  use Chart::Graph::Gnuplot qw(gnuplot);

  gnuplot({'title' => 'foo',
           'x2-axis label' => 'bar',
           'logscale x2' => '1',
           'logscale y' => '1',
           'output type' => 'png',
           'output file' => 'gnuplot1.png',
           'xtics' => [ ['small\nfoo', 10], ['medium\nfoo', 20], ['large\nfoo', 30] ],
           'ytics' => [10,20,30,40,50],
           'extra_opts' => 'set key left top Left'},
          [{'title' => 'data1',
            'type' => 'matrix'}, [[1, 10],
                                  [2, 20],
                                  [3, 30]] ],
          [{'title' => 'data2',
            'style' => 'lines',
            'type' => 'columns'}, [8, 26, 50, 60, 70],
                                  [5, 28, 50, 60, 70] ],
          [{'title' => 'data3',
            'style' => 'lines',
            'type' => 'file'}, './samplefile'],);

=for html
<p><center><img src="http://www.caida.org/tools/utilities/graphing/gnuplot1.png"></center></p>
<p><center><em>gnuplot1.png</em></center></p>

=head2 ERRORBARS

I<Gnuplot> supports errorbars to aid in data interpretation.  To use an
arbitrary number of data columns (for errorbars), set C<style> to
C<errorbars> and include extra data columns.  The example below
produces the file F<gnuplot2.png>

Note the following: These columns MUST be the the following formats:
C<(x, y, ydelta)>, C<(x, y, ylow, yhigh)>, C<(x, y, xdelta)>, C<(x, y,
xlow, xhigh)>, C<(x, y, xdelta, ydelta)>, or C<(x, y, xlow, xhigh,
ylow, yhigh)> This will only work with data type C<columns>. Also, you
MUST use the C<using> option to specify how columns of the data file
are to be assigned to C<x>, C<y>, C<ydelta>, C<ylow> and C<yhigh>,
C<xdelta>, C<xlow> and C<xhigh>.

     use Chart::Graph::Gnuplot qw(gnuplot);

     gnuplot({"title" => "Examples of Errorbars",
              "xrange" => "[:11]",
              "yrange" => "[:45]",
              "output file" => "gnuplot2.gif",
	      "output type" => "gif",
             },
             # dataset 1
             [{"title" => "yerrorbars",
               "style" => "yerrorbars",
               "using" => "1:2:3:4",
               "type" => "columns"},
              [ 1, 2, 3, 4, 5, 6 ], # x
              [ 5, 7, 12, 19, 28, 39 ], # y
              [ 3, 5, 10, 17, 26, 38 ], # ylow
              [ 6, 8, 13, 20, 30, 40 ] ], # yhigh
             # dataset 2
             [{"title" => "xerrorbars",
               "style" => "xerrorbars",
               "using" => "1:2:3:4",
               "type" => "columns"},
              [ 4, 5, 6, 7, 8, 9 ], # x
              [ 1, 4, 5, 6, 7, 10 ], # y
              [ 3.3, 4.4, 5.5, 6.6, 7.7, 8.8 ], # xlow
              [ 4.1, 5.2, 6.1, 7.3, 8.1, 10 ] ], # xhigh
             # dataset 3
             [{"title" => "xyerrorbars",
               "style" => "xyerrorbars",
               "using" => "1:2:3:4:5:6",
               "type" => "columns"},
              [ 1.5, 2.5, 3.5, 4.5, 5.5, 6.5 ], # x
              [ 2, 3.5, 7.0, 14, 15, 20 ], # y
              [ 0.9, 1.9, 2.8, 3.7, 4.9, 5.8 ], # xlow
              [ 1.6, 2.7, 3.7, 4.8, 5.6, 6.7 ], # xhigh
              [ 1, 2, 3, 5, 7, 8 ], # ylow
              [ 5, 7, 10, 17, 18, 24 ] ], # yhigh
             # dataset 4
             [{"title" => "xerrorbars w/ xdelta",
               "style" => "xerrorbars",
               "using" => "1:2:3",
               "type" => "columns"},
              [ 4, 5, 6, 7, 8, 9 ], # x
              [ 2.5, 5.5, 6.5, 7.5, 8.6, 11.7 ], # y
              [ .2, .2, .1, .1, .3, .3 ] ], # xdelta
             # dataset 5
             [{"title" => "yerrorbars w/ ydelta",
               "style" => "yerrorbars",
               "using" => "1:2:3",
               "type" => "columns"},
              [ .7, 1.7, 2.7, 3.7, 4.7, 5.7 ], # x
              [ 10, 15, 20, 25, 30, 35 ], # y
              [ .8, 1.2, 1.1, 2.1, 1.3, 3.3 ] ], # ydelta
             # dataset 6
             [{"title" => "dummy data",
               "type" => "matrix"},
              [ [1,1] ]],
             # dataset 7
             [{"title" => "xyerrorbars w/ xydelta",
               "style" => "xyerrorbars",
               "using" => "1:2:3:4",
               "type" => "columns"},
               [ 7.5, 8.0, 8.5, 9.0, 9.5, 10.0 ], # x
               [ 30, 27, 25, 23, 27, 33 ], # y
               [ .2, .1, .3, .6, .4, .3 ], # xdelta
              [ .8, .7, .3, .6, 1.0, .3 ] ], # ydelta
           );

=for html
<p><center><img src="http://www.caida.org/tools/utilities/graphing/gnuplot2.gif"></center></p>
<p><center><em>gnuplot2.gif</em></center></p>

=head2 PLOTTING DATES - CUSTOM GNUPLOT OPTIONS

As noted above, I<Chart::Graph::Gnuplot> includes support for plotting
date and times as source data.  the following shows how to plot data,
where the x-axis contains dates, and the y-axis contains stock prices
from a major computer major during the "dot-com meltdown." For date
and time data that requires high precision, using UNIX time stamps is
probably the best solution (see below.)  As used in the first example,
any option available to I<Gnuplot> can be passed to I<Gnuplot> using
the C<extra_opts> option. This example uses this feature to enable two
options: a grid over the graph and a timestamp for when the graph was
created.

  use Chart::Graph::Gnuplot qw(gnuplot);

  #Debugging aid - save the temporary files if desired
  #$Chart::Graph::save_tmpfiles = 1;
  #Debugging aid - turn on extra debugging messages
  #$Chart::Graph::debug = 1; 

  # Call and "usual" global parameters

  gnuplot({'title' => 'Corporate stock values for a major computer maker',
           'x-axis label' => 'Month and Year',
           'y-axis label' => 'Stock price',
           'output type' => 'png',
           'output file' => 'gnuplot3.png',
           # Setting date/time specific options.
           'xdata' => 'time',
           'timefmt' => '%m/%d/%Y',
           'format' => ['x', '%m/%d/%Y'],
           # Set output range - note quoting of date string
           'xrange' => '["06/01/2000":"08/01/2001"]',
           'extra_opts' => join("\n", 'set grid', 'set timestamp'),
          },
          # Data for when stock opened
          [{'title' => 'open',
            'type' => 'matrix',
            'style' => 'lines',
           },
           [
            ['06/01/2000',  '81.75'],
            ['07/01/2000', '52.125'],
            ['08/01/2000', '50.3125'],
            ['09/01/2000', '61.3125'],
            ['10/01/2000', '26.6875'],
            ['11/01/2000', '19.4375'],
            ['12/01/2000', '17'],
            ['01/01/2001', '14.875'],
            ['02/01/2001', '20.6875'],
            ['03/01/2001', '17.8125'],
            ['04/01/2001', '22.09'],
            ['05/01/2001', '25.41'],
            ['06/01/2001', '20.13'],
            ['07/01/2001', '23.64'],
            ['08/01/2001', '19.01'],
           ]
          ],

          # Data for stock high
          [{'title' => 'high',
            'type' => 'matrix',
            'style' => 'lines',
           },
           [
            ['06/01/2000', '103.9375'],
            ['07/01/2000', '60.625'],
            ['08/01/2000', '61.50'],
            ['09/01/2000', '64.125'],
            ['10/01/2000', '26.75'],
            ['11/01/2000', '23'],
            ['12/01/2000', '17.50'],
            ['01/01/2001', '22.50'],
            ['02/01/2001', '21.9375'],
            ['03/01/2001', '23.75'],
            ['04/01/2001', '27.12'],
            ['05/01/2001', '26.70'],
            ['06/01/2001', '25.10'],
            ['07/01/2001', '25.22'],
            ['08/01/2001', '19.90'],
           ]
          ],

          # Data for stock close
          [{'title' => 'close',
            'type' => 'matrix',
            'style' => 'lines',
           },
           [

            ['06/01/2000', '52.375'],
            ['07/01/2000', '50.8125'],
            ['08/01/2000', '60.9375'],
            ['09/01/2000', '25.75'],
            ['10/01/2000', '19.5625'],
            ['11/01/2000', '16.50'],
            ['12/01/2000', '14.875'],
            ['01/01/2001', '21.625'],
            ['02/01/2001', '18.25'],
            ['03/01/2001', '22.07'],
            ['04/01/2001', '25.49'],
            ['05/01/2001', '19.95'],
            ['06/01/2001', '23.25'],
            ['07/01/2001', '18.79'],
            ['08/01/2001', '18.55'],
           ]
          ]
);



=for html
<p><center><img src="http://www.caida.org/tools/utilities/graphing/gnuplot3.png"></center></p>
<p><center><em>gnuplot3.png</em></center></p>

=head2 UNIX TIMESTAMPS

I<Chart::Graph::Gnuplot> can convert Unix timestamps into normal dates
for x-axis values. Collisions with existing user x-tics are can be
remedied by prepending a literal '\n' (or "\\n") to their tic-labels.
The 'uts' option takes an array ref with 2 to 4 elements:
[ start_timestamp, end_timestamp, <scale>, <use_local_timezone> ]
If the optional element 'scale' is > 1 the number of tics will be reduced.
If the optional element 'use_local_timezone' is set to non-zero value
the local timezone is used, UTC is assumed otherwise.
The variables I<$Chart::Graph::Gnuplot::show_year> and 
I<$Chart::Graph::Gnuplot::show_seconds> influence the formatting of the x-tics.

    [...]

    %options = (
                 'title' => 'uts example',
                 'output file' => 'gnuplot4.gif',
                 'output type' => 'gif',
                 'x2-axis label' => 'time',
                 'xtics' => [ ['\n9pm UTC', 954795600] ],
                 'ytics' => [10,20,30,40,50],
                 'extra_opts' => 'set nokey',
                 'uts' => [954791100, 954799300],
               );

    $plot = [{'title' => 'Your title',
              'type' => 'matrix'},
              [
                [954792100, 10],
                [954793100, 18],
                [954794100, 12],
                [954795100, 26],
                [954795600, 13], # 21:00
                [954796170, 23],
                [954797500, 37],
                [954799173, 20],
                [954799300, 48],
              ],
            ];

    gnuplot(\%options, $plot);

=for html
<p><center><img src="http://www.caida.org/tools/utilities/graphing/gnuplot4.gif"></center></p>
<p><center><em>gnuplot4.gif</em></center></p>


B<Note:> The present implementation of UNIX time stamps only supports
assigning I<xtics> for x-axis labels.  Using the Gnuplot directive:
C<format> is not supported.

=head2 FUNCTIONS

I<Chart::Graph::Gnuplot> supports the plotting of functions, this can
be mixed with other data-types:

   my %options = (
                   'title' => 'plot functions example',
                   'output file' => 'gnuplot5.png',
                 );

   my $data = [{ 'title' => 'data 1',
                 'style' => 'lines',
                 'type' => 'matrix',
               },
               [
                 [0,10],
                 [3,30],
                 [6,0],
                 [9,-10],
                 [12,-0],
               ]
              ];

   my $fnc1 = [{ 'title' => 'function 1',
                 'style' => 'lines',
                 'type' => 'function',
               },
               '10*sin(x)+2*cos(1.1 * x)+.5*tan(x)'
              ];

   my $fnc2 = [{ 'title' => 'function 2',
                 'style' => 'lines',
                 'type' => 'function',
               },
              '20*sin(sqrt(2**x))/sqrt(2**x)'
              ];

    gnuplot(\%options, $data, $fnc1, $fnc2);


=for html
<p><center><img src="http://www.caida.org/tools/utilities/graphing/gnuplot5.png"></center></p>
<p><center><em>gnuplot5.png</em></center></p>


=head1 MORE INFO

This version of I<Chart::Graph::Gnuplot> was tested against
Gnuplot Version 4.0 patchlevel 0, some features might not work
on older versions of gnuplot. 

For more information on gnuplot, please see the gnuplot web page:

 http://www.gnuplot.org/

=head1 CONTACT

Send email to graph-dev@caida.org is you have problems, questions,
or comments. To subscribe to the mailing list send mail to
graph-dev-request@caida.org with a body of "subscribe your@email.com"

=head1 AUTHOR

 CAIDA Perl development team (cpan@caida.org)

=head1 SEE ALSO

 gnuplot(1).
 

=cut

## Xmgrace.pm is a sub-module of Graph.pm. It has all the subroutines 
## needed for the gnuplot part of the package.
##
## $Id: Xmgrace.pm,v 1.34 2006/06/07 21:09:33 emile Exp $ $Name:  $
##
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
package Chart::Graph::Xmgrace;
use Exporter();

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(&xmgrace);

use Carp;			# for carp() and croak()
use Chart::Graph::Utils qw(:UTILS);	# get global subs and variable
use Chart::Graph::XrtUtils qw(:UTILS);
use Chart::Graph::Xmgrace::Grace;
use Chart::Graph::Xmgrace::Graph_Presentation_Type;
use Chart::Graph::Xmgrace::Dataset;
use FileHandle;

$cvs_Id = '$Id: Xmgrace.pm,v 1.34 2006/06/07 21:09:33 emile Exp $';
$cvs_Author = '$Author: emile $';
$cvs_Name = '$Name:  $';
$cvs_Revision = '$Revision: 1.34 $';

$VERSION = 3.2;

use strict;

#for debugging purposes.
my $stdout = 0;

# these variables hold default options for xmgrace
my %def_xmgrace_global_opts = (
			       "title" => "untitled", # @description
			       "subtitle" => "",
			       "type of graph" => "XY graph", # XY graph, XY chart, Polar Graph, Smith Chart, Fixed
			       "output type" => "PNG", # png, jpeg, ps
   			       "output file" => "untitled-grace",
			       "grace output file" => "untitled-grace.agr",
			       "xrange" => undef,
			       "yrange" => undef,
			       "x-axis label" => "x-axis",
			       "y-axis label" => "y-axis",
			       "alt x-axis label" => undef,
			       "alt y-axis label" => undef,
			       "logscale x" => undef,
			       "logscale y" => undef,
			       "xtics" => undef,
			       "ytics" => undef,
			       "alt xtics" => undef,
			       "alt ytics" => undef,
			       "stacked" => "false",
			       "extra opts" => undef,
			      ); 


my $def_graph_appearance = new Chart::Graph::Xmgrace::Graph_Presentation_Type;

my %def_xmgrace_data_opts = (
			     "set presentation" => "XY",
			     "options" => $def_graph_appearance->{"XY graph"},
			     "title" => "untitled data set", # comment, legend
			     "data format" => undef, # columns, matrix, or file
			    );

# New Hash of array references to hold list of options that might need to
# be checked against options.  If Xmgrace program changes.  Update here
my %def_xmgrace_available_options = 
    ( "type of graph" => ["XY graph", "XY chart", "Bar chart",
			  "Polar graph", "Smith chart", "Fixed",
			  # "Pie Chart", # Not yet implemented
			  ],
    );


#
#
# Subroutine: xmgrace()
#
# Description: this is the main function you will be calling from
#              our scripts. please see 
#              www.caida.org/Tools/Graph/graph_xmgrace.html for a 
#              full description and how-to of this subroutine
#

sub xmgrace {
    my ($user_global_opts_ref, @data_sets) = @_;
    my (%data_opts, %global_opts,);
    my ($plottype, $output_file, $plot_file, $output_type, $data_set_ref);
    #my $autoscale = 1;		# by default use autoscale
    #my ($xscale, $yscale) = (0,0);
    my $autoscale = "";

    # create a new filehandle to be used throughout package
    my $handle = new FileHandle;
    my $grace_output_file;

    my $grace = new Chart::Graph::Xmgrace::Grace("g0");

    # create tmpdir
    _make_tmpdir("_Xmgrace_");	# grace files should be saved for user tweaking
    
    # set paths for external programs
    if (not _set_xmgrace_paths()) {
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
    
    my $command_file = _make_tmpfile("command", "agr");
    # remember to close the file if we return
    if (not $handle->open(">$command_file")) {
	carp "could not open file: $command_file";
	_cleanup_tmpdir();
	return 0;
    }
    
    # process global options
    # call to combine user options with default options
    %global_opts = _mesh_opts($user_global_opts_ref, \%def_xmgrace_global_opts);
    
    # crunch on the global options
    while (my ($key, $value) = each %global_opts) {
	
	if ($key eq "title") {
	    $grace->title_options->title($value);
	}
	
	if ($key eq "subtitle") {
	    $grace->subtitle_options->title($value);
	}
	
	if ($key eq "type of graph") {
	    unless(_is_available_option($key, $value,
					\%def_xmgrace_available_options)) {
		carp "Not a valid Xmgrace graph type";
		$handle->close();
		_cleanup_tmpdir();
		return 0;
	    }
	    if ($value =~ "graph") { 
		$grace->graph_global_options->type("XY");
	    } elsif ($value =~ "chart") {
		$grace->graph_global_options->type("Chart");
	    }
	}
	
	if ($key eq "stacked") {
	    if ($value eq "true" || $value eq "false") {
		$grace->graph_global_options->stacked("$value");
	    } else {		# anything else is considered "false"
		carp join ('',
			   "Warning: ambiguous setting for stacking",
			   "bar chart\nAssuming no stacking"
			  );
		$grace->graph_global_options->stacked("false");
	    }
	}
	
	if ($key eq "x-axis label") {
	    $grace->axes_options->xaxis->label_options->label($value);
	}
	
	if ($key eq "y-axis label") {
	    $grace->axes_options->yaxis->label_options->label($value);
	}
	
	if ($key eq "alt x-axis label") {
	    if (defined($value)) {
		$grace->axes_options->altxaxis->status("on");
		$grace->axes_options->altxaxis->label_options->label($value);
	    }
	    # else altxaxis  off
	}
	
	if ($key eq "alt y-axis label") {
	    if (defined($value)) {
		$grace->axes_options->altyaxis->status("on");	
		$grace->axes_options->altyaxis->label_options->label($value);
	    }
	    # else altyaxis  off
	}


	if ($key eq "xrange") {
	    if (defined($value)) {
		if (ref($value) ne "ARRAY") {
		    carp "xrange values must be an ARRAY\n";
		    $handle->close();
		    _cleanup_tmpdir();
		    return 0;
		}
		$grace->world_options->xmin($value->[0]);
		$grace->world_options->xmax($value->[1]);
		$autoscale = "x";
	    }
	}

	if ($key eq "yrange") {
	    if (defined($value)) {
		if (ref($value) ne "ARRAY") {
		    carp "yrange values must be an ARRAY\n";
		    $handle->close();
		    _cleanup_tmpdir();
		    return 0;
		}
		$grace->world_options->ymin($value->[0]);
		$grace->world_options->ymax($value->[1]);
		$autoscale = "y";
	    }
	}
	
	if ($key eq "xtics") {
	    if (defined($value)) {
		$grace->axes_options->xaxis->tick_options->type("spec");
		$grace->axes_options->xaxis->ticklabel_options->type("spec");
		$grace->axes_options->xaxis->ticklabel_options->ticklabels($value);
	    }
	}
	
	
	if ($key eq "ytics") {
	    if (defined($value)) {
		$grace->axes_options->yaxis->tick_options->type("spec");
		$grace->axes_options->yaxis->ticklabel_options->type("spec");
		$grace->axes_options->yaxis->ticklabel_options->ticklabels($value);
	    }
	}
	
	if ($key eq "alt xtics") {
	    if (defined($value)) {
		$grace->axes_options->altxaxis->tick_options->type("spec");
		$grace->axes_options->altxaxis->ticklabel_options->type("spec");
		$grace->axes_options->altxaxis->ticklabel_options->ticklabels($value);
	    }
	}
	
	
	if ($key eq "alt ytics") {
	    if (defined($value)) {
		$grace->axes_options->altyaxis->tick_options->type("spec");
		$grace->axes_options->altyaxis->ticklabel_options->type("spec");
		$grace->axes_options->altyaxis->ticklabel_options->ticklabels($value);
	    }
	}

	if ($key eq "logscale x") {
	    if (defined($value)) {
		$grace->axes_options->xaxes->scale("Logarithmic");
	    }
	}
	
	if ($key eq "logscale y") {
	    if (defined($value)) {
		$grace->axes_options->yaxes->scale("Logarithmic");
	    }
	}

	if ($key eq "output type") {
	    $output_type = uc($value);
	    if (not _check_output_type($output_type)) {
		$output_type = $def_xmgrace_global_opts{"output type"};
	    }
	    #print "output type = $value\n";
	    $grace->output_type($output_type);
	}
	
	if ($key eq "output file") {
	    #print "output file = $value\n";
	    $grace->output_file($value);	    
	}	
	
	if ($key eq "grace output file") {
	    if (defined($value)) {
		$grace->grace_output_file(_make_gracefile("$value"));
		#print "grace output file = $value\n";
	    }
	}

	if ($key eq "extra opts") {
	    if (defined($value)) {
		$grace->extra_options->extras($value);
	    }
	}
    }
    
    # process data sets
    my @cum_data_set_objects;
    my $dataset_count = 0;
    
    while (@data_sets) {
	$data_set_ref = shift @data_sets;
	
	if (ref($data_set_ref) ne "ARRAY") {
	    carp "Data set must be an array";
	    $handle->close();
	    _cleanup_tmpdir();
	    return 0;
	}
	
	# create a new Dataset object with each new datase
	my $data_set_object = new Chart::Graph::Xmgrace::Dataset; 
	my ($user_data_opts_ref, @data) = @{$data_set_ref};
	
	# set values in the dataset object
	$data_set_object->set_number($dataset_count);
	$data_set_object->data(\@data);
	
	# process data
	$data_set_object->data_format($user_data_opts_ref->{"data format"});
	my $formatted_data = _xmgrace_data_set($data_set_object);
	
	if (not $formatted_data) {
	    # error message already printed
	    $handle->close();
	    _cleanup_tmpdir();
	    return 0;
	}
	
	# stick formatted data back in the dataset object
	$data_set_object->data($formatted_data);

	# set set_type
	my $tog = $global_opts{"type of graph"}; # type of graph
	if ($tog =~ "XY") {
	    $data_set_object->set_type("XY");
	} elsif ($tog =~ m/BAR/i) {
	    $data_set_object->set_type("BAR");
	} else {
	    carp "Untrapped graph type - Xmgrace.pm internal error";
	}
		
	# need to create a new Presentation Type object for each dataset
	my $gpt = new Chart::Graph::Xmgrace::Graph_Presentation_Type($dataset_count);
	my $graph_appearance = $gpt->{$global_opts{"type of graph"}};
	$def_xmgrace_data_opts{options} = $graph_appearance;	
	
	# mesh data options
	my $data_options_ref = _mesh_xmgrace_opts($user_data_opts_ref, 
						  \%def_xmgrace_data_opts);
	
	if (not $data_options_ref) {
	    # error message already printed
	    $handle->close();
	    _cleanup_tmpdir();
	    return 0;
	}

	# change the "title" into "comment" and "legend"
	_set_title($data_options_ref);

	# set data options in the dataset object
	$data_set_object->options($data_options_ref);

	# create new set presention for each data set
	my $sp = $user_data_opts_ref->{"set presentation"};
	if ($sp) {
	    if ($sp eq "XY") {
		_set_XY($data_set_object);
	    } elsif ($sp =~ m/BAR/i) {
		_set_BAR($data_set_object);
	    } else {
		carp "Untrapped graph type - Xmgrace.pm internal error";
	    }
	}

	# process data
	# accumulate data objects into an array.
        push @cum_data_set_objects, $data_set_object;
	$dataset_count++;
    }

    # print commandfile
    if ($stdout) { 
	$handle = \*STDOUT; 
    }
    $grace->print($handle);
 
    # we're done crunching on the data now print something!    
    # print data options


     foreach my $set_object (@cum_data_set_objects) {
#       XXX This code breaks under Perl 5.6 - not sure why XXX
# 	$set_object->{options}->{options}->print($handle,"s" . 
# 						 "$set_object->{set_number}");
     }

    # print data sets
    foreach my $set_object (@cum_data_set_objects) {
	_printline($handle, "target G0.S$set_object->{\"set_number\"}\n");
	_printline($handle, "type $set_object->{\"set_type\"}\n");
	if (not _print_data_set($handle, $set_object)) {
	    # already printed error message
	    return 0;
	}
    }

    # commandfile successfully completed    
    $handle->close unless ($stdout);

    # if user chooses to, can save the .agr file
    if (defined $grace->{grace_output_file}) {
	_copy_commandfile($command_file, $grace->{grace_output_file});
    }

    # run xmgrace with commandline args to create the specified graph type
    # $xmgrace needs to be defined in the user's script    

    _exec_xmgrace($command_file, $grace, $autoscale); # could probably put the autoscale flag inside the grace object


}				# end of xmgrace()


############################################################
#
# The following functions are internal functions used by this package. The
# user won't know about these functions.
#
############################################################

#
#
# Subroutine _is_available_option($$$)
#
# Description: Modest helper subroutine to loop through available
#              options for more complex Xmgrace options and return
#              true if a match is found.
#
#
sub _is_available_option($$$) {
    my $key = shift;
    my $to_match = shift;
    my $avail_options_ref = shift;

    if (exists($avail_options_ref->{$key})) {
	foreach my $option (@{$avail_options_ref->{$key}}) {
	    if ($to_match =~ m/$option/i) { 
		return($option);
	    }
	}
    }
    return(undef);
}



sub _exec_xmgrace ($$$) {
    my ($command_file, $grace, $autoscale) = @_;
    my ($childpid, $port);
    my $display_env = $ENV{DISPLAY};
    my $status;
    my $framebuffer;

    #if ($Chart::Graph::use_xvfb) {
	# start the virtual X server
	#($childpid, $port) = _exec_xvfb();
#	$status = system("$xmgrace -display :$port.0 < $command_file");
	
    #} else {
	# use the local X server
	# warning: colors might be messed up 
	# depending on your current setup
	#$status = system("$Chart::Graph::xrt3d -display $display_env < $command_file");

    #}

    #my $status = system("$xrt -display :$port.0 < $command_file");
    #if (not _chk_status($status)) {
#	return 0;
#    }

    #if ($Chart::Graph::use_xvfb) {
#	kill('KILL', $childpid);
 #   }


    if ($autoscale eq "") {
	_run_app("$xmgrace", "-hardcopy -hdevice $grace->{output_type}", 
		 "-printfile $grace->{output_file} $command_file -pexec autoscale") unless $stdout;
	if ($Chart::Graph::debug) {
	    print STDERR "$xmgrace -hardcopy -hdevice $grace->{output_type} -printfile $grace->{output_file} $command_file -pexec autoscale\n";
	    print STDERR "done.\n";
	    ### fill in tmp files...
	}
	1;
    } elsif ($autoscale eq "x") {
	_run_app("$xmgrace", "-autoscale y", "-hardcopy -hdevice $grace->{output_type}", 
		 "-printfile $grace->{output_file} $command_file -pexec autoscale") unless $stdout;
	print STDERR "xscale\n";
    } elsif ($autoscale eq "y") {
	_run_app("$xmgrace", "-autoscale x", "-hardcopy -hdevice $grace->{output_type}", 
		 "-printfile $grace->{output_file} $command_file -autoscale") unless $stdout;
	print STDERR "yscale\n";
    } else {
	_run_app("$xmgrace", "-autoscale xy", "-hardcopy -hdevice $grace->{output_type}", 
		 "-printfile $grace->{output_file} $command_file -autoscale") unless $stdout;
    	print STDERR "both\n";
    }
}


# 
# Subroutine: set_xmgrace_paths()
# 
# Description: set paths for external programs required by gnuplot()
#              if they are not defined already
#

sub _set_xmgrace_paths {

    if (not defined($xmgrace)) {
	if (not $xmgrace = _get_path("xmgrace")) {
	    return 0;
	}
    }

    if (not defined($xvfb)) {
	if (not $xvfb = _get_path("Xvfb")) {
	    return 0;
	}
    }
    
    # make sure /usr/dt/lib is in the library path
    _set_ldpath("/usr/dt/lib");	# defined in Xrtutils.pm

    return 1;
}


#
# Subroutine: _xmgrace_data_set()
#
# Description: this function sets up the data sets into a columns format
#              and returns a formatted data reference
#
#

sub _xmgrace_data_set {
    my $data_set_object = shift;

    my $data_format = $data_set_object->{"data_format"};
    my @data = @{$data_set_object->data};
    my $formatted_data_ref;
    
    # we give the user 3 formats for supplying the data set
    # 1) matrix
    # 2) column
    # 3) file
    # please see the online docs for a description of these 
    # formats

    if ($data_format eq "matrix") {
	$formatted_data_ref = _matrix_to_columns(@data);
    } elsif ($data_format eq "columns") {
	$formatted_data_ref = _columns_to_columns(@data);
    } elsif ($data_format eq "file") {
	$formatted_data_ref = _file_to_columns(@data); # single item list (filename)
    } elsif ($data_format eq "") {
	carp "Need to specify data set type";
	return 0;
    } else {
	carp "Illegal data format: $data_format"; 
        return 0;
    }
    return $formatted_data_ref;
}

#
#
# Subroutine: matrix_to_columns()
#
# Description: converts the matrix data input into a an array 
#              of pairs of input. See www for the specific on 
#              the matrix format
# 
# 

sub _matrix_to_columns (\@ ) {
    my ($matrix_ref) = @_;
    my $entry_ref;
    my $matrix_len;
    my @pairs = ();
    my $single_pair = "";
    my (@x_col, @y_col);

    if (ref($matrix_ref) ne "ARRAY") {
	carp "Matrix data must be a reference to an array";
	return 0;
    }
    
    $matrix_len = @{$matrix_ref};
    for (my $i = 0; $i < $matrix_len; $i++) {
	$entry_ref = $matrix_ref->[$i];
	
	if (ref($entry_ref) ne "ARRAY") {
	    carp "Matrix entry must be a reference to an array";
	    return 0;
	}
	
	# check that each entry ONLY has two entries
	if (@{$entry_ref} != 2) {
	    carp "Each entry must be an array of size 2";
	    return 0;
	}
	
	push @x_col, $entry_ref->[0];
	push @y_col, $entry_ref->[1];

	@pairs = (\@x_col, \@y_col);
    }

    return \@pairs;
}

#
#
# Subroutine: columns_to_columns()
#
# Description: converts the column data input into a the gnuplot
#              data file format. please see www page for specifics
#              on this format.
#

sub _columns_to_columns (\@ ) {
    my ($x_col, $y_col) = @{$_[0]};
    my ($x_len, $y_len);
    my @pairs = ();
    my $single_pair = "";
      
    if ((ref($x_col) ne "ARRAY") or (ref($y_col) ne "ARRAY")) {
	unless($x_col) {
	    carp "No data in X column";
	}
	unless($y_col) {
	    carp "No data in Y column";
	}
	carp "Column data must be a reference to an array";
	return 0;
    }
    
    $x_len = @{$x_col};
    $y_len = @{$y_col};
  
    if ($x_len != $y_len) {
	carp "x and y columns must be of same length";
	return 0;
    }

    if ($x_len == 0) {
	carp "Warning: x column has no data";
    }
    
    if ($y_len == 0) {
	carp "Warning: y column has no data";
    }
    
    @pairs = ($x_col, $y_col);   
    return \@pairs;
}

#
# Subroutine: file_to_columns()
#
# Description: If a gnuplot data set was given in
#              file format, we simply copy the data 
#              and read it into 
#

sub _file_to_columns (\@ ) {
    my ($file_in) = @_;
    my $single_pair = "";
    my @pairs = ();
    my (@x_col, @y_col);
    my ($x_len, $y_len);

    if (not $file_in) {
	carp "File data format selected but data set file missing";
	return 0;
    }

    if (not -f $file_in) {
	carp "Data set file, '$file_in', does not exist.";
	return 0;
    }
    
    my $fh = new FileHandle;
    
    if (not $fh->open("<$file_in")) {
	carp "Could not open file: '$file_in'";
	cleanup_tmpdir();
	return 0;
    }

    while (<$fh>) {
	chomp;
	my @data = split /\s+/;
	push @x_col, $data[0];
	push @y_col, $data[1];
    }
    
    $x_len = @x_col;
    $y_len = @y_col;

    if ($x_len != $y_len) {
	carp "x and y columns must be of same length";
	return 0;
    }
    
    if ($x_len == 0) {
	carp "Warning: x column has no data";
    }
    
    if ($y_len == 0) {
	carp "Warning: y column has no data";
    }

    @pairs = (\@x_col, \@y_col); 
    $fh->close;
    return \@pairs;
}

#
# Subroutine: _printline()
#
# Description: prints out an Xmgrace formatted line
#
#
 
sub _printline ($$$ ) { 
    my ($handle, $string, $length) = @_;

    unless ($length) {$length = 1;}
    print $handle "@";
    print $handle ' ' x $length;
    print $handle "$string"; 

    return 1;                   # just for fun 
 
} 

#
# Subroutine: _copy_commandfile()
#
# Description: copies the temp commandfile into a grace 
#              output file (.agr)
#

sub _copy_commandfile ($$ ) {
    
    my ($commandfile, $grace_output_file) = @_;
    my $status = system("cp", "$commandfile", "$grace_output_file");
    if (not _chk_status($status)) {
	return 0;
    }
}

#
# Subroutine: _make_gracefile()
#
# Description: constructs the gracefile.agr name 
#

sub _make_gracefile ($ ) {
    
    my $file = shift;
    my $ext = ".agr";
    
    if ($file !~ m|\.agr$|) {
	return "$file$ext";
    } else {
	return "$file";
    }
}

#
# Subroutine: _check_output_type()
#
# Description: checks user output type to known output types. 
#

sub _check_output_type ($ ) {
    
    my $user_type = shift;
    my @types = ("PNG", "PS", "JPEG", "PNM");
    my $seen = 0;
    
    foreach my $known_type (@types) {
	if ($user_type =~ m|$known_type|i)  {
	    $seen = 1;
	}
	if ($seen) {
	    return 1;
	}
    }
    carp "Unknown output type: \'$user_type\'\nUsing default output type: \'PNG\'";    
    return 0;
}


# 
# Subroutine: _run_app() 
#
# Description: to simplify calls to external applications
# 

sub _run_app($@ ) {
  my $application = shift;
  my @arguments = @_;
  my $command_str;

        # Run application
        # Test if an executable file exists as specified location
        if (-x "$application") {
          # Try to run application.  Quit if that fails.
          $command_str = join(' ', $application, @arguments);
          my $rc = 0xffff & system ($command_str);
          if ($rc != 0) {
	      die "Execution error: $application";
	  }
        } else {
          die "No such application: $application";
        }
}

#
# Subroutine: _mesh_xmgrace_opts
#
# Description: meshes the user's data_options with the default 
#              data options these "options" include: title, set 
#              presentation, options, and data format             
#

sub _mesh_xmgrace_opts {
    my ($user_opts_ref, $default_opts_ref) = @_;

    my %user_opts = %{$user_opts_ref};
    my %default_opts = %{$default_opts_ref};
    my %opts;

    # check user opts against defaults and mesh
    # the basic algorithm here is to override the 
    # the default options against the ones that
    # the user has passed in. 
    while (my ($key, $value) = each %default_opts) {
	
	if ($key eq "options") {
	    my $opts_obj = $def_xmgrace_data_opts{options};
	    my $def_options = $opts_obj->{options};
	    
	    if (defined($user_opts{options})) {
		# if user provides options for the datasets we mesh them
		my $data_options_obj = _mesh_option_opts($user_opts_ref->{options}, 
							 $def_options);
		$opts_obj->{options} = $data_options_obj;
		$opts{options} = $opts_obj;
		delete $user_opts{$key}; # remove options
	    } else {
		$opts_obj->{options} = $def_options;
		$opts{options} = $opts_obj;	    
	    }
	} else {
	    if (defined($user_opts{$key})) {
		$opts{$key} = $user_opts{$key};
		delete $user_opts{$key}; # remove options 
		# that are matching
	    } else {
		$opts{$key} = $default_opts{$key};
	    }
	}
    }
    
    # any left over options in the table are unknown
    # if the user passes in illegal options then we 
    # warn them with an error message but still 
    # proceed.
    while (my ($key, $value) = each %user_opts) {
	carp "unknown option: $key";
    }
    
    return \%opts;
}

#
# Name: _mesh_option_opts
#
# Description: if the user introduces a "option" in their dataset 
#              hash, we use this function to mesh those options 
#              with the default options of a particular dataset. 
#              note: default options depend on the type of
#              graph
#

sub _mesh_option_opts {
    my ($user_opts_ref, $default_opts_ref) = @_;
    my %user_opts = %{$user_opts_ref};
    my %def_options = %{$default_opts_ref};
    
    # check user opts against defaults and mesh
    # the basic algorithm here is to override the 
    # the default options against the ones that
    # the user has passed in. 
    while (my ($key, $value) = each %def_options) {
	if (defined($user_opts{$key})) {
	    
	    if (!ref($value)) {	# scalar
		$def_options{$key} = $user_opts{$key};
		delete $user_opts{$key}; # remove options 	   
	    }
	    
	    else {		# blessed object
		my $return_obj;
		my $class = ref($def_options{$key});
		$return_obj = _mesh_option_opts($user_opts_ref->{$key}, 
						$def_options{$key}->{options});
		bless($return_obj, $class);
		$def_options{$key}->{options} = $return_obj;
		delete $user_opts{$key};
	    }
	}
    }
    
    # any left over options in the table are unknown
    # if the user passes in illegal options then we 
    # warn them with an error message but still 
    # proceed.
    while (my ($key, $value) = each %user_opts) {
	carp "unknown option: $key";
    }
    
    my $ret_obj = \%def_options;
    return $ret_obj;
}

#
# Name: _print_data_set
#
# Description: used to print out the data in the dataset in a 
#              formatted manner
#

sub _print_data_set ($\@ ) {
    my ($handle, $set_object) = @_;
    my ($x_col, $y_col);
    my $length;
    my $i;
    my $delimeter = "\t\t"; # set to 15 spaces, nice one xmgrace!
    my $pairs = $set_object->data;

    ($x_col, $y_col) = ($pairs->[0], $pairs->[1]);
    $length = $#{$x_col};

    #
    # data set format
    # @target G<graph #>.S<set #>
    # @type $type_of_graph
    #                x1                y1
    #                 .                 .
    #                 .                 .
    #                xn                yn
    # &
    #

    # we don't print out the hidden data set
    # XXX Protection against trying to access nonexistentant 
    # XXX accessor function or OVERLOADED function?
    if (exists($set_object->{options}->{options}->{hidden})) {
	# This statement now breaks in 5.6 - unblessed reference.
	if ($set_object->{options}->{options}->hidden eq "true") {
	    return 1;
	}
    }

    for ($i = 0; $i <= $length; $i++) {
        print $handle "$delimeter", "$x_col->[$i]", "$delimeter",
                      "$y_col->[$i]\n";
    }

    # each data set is separated by a "&"
    print $handle "&\n";
    return 1;
}                                          


#
# Name: _set_BAR
#
# Description: sets up the necessary characteristics of a BAR
#              dataset type
#

sub _set_BAR ($ ) {
    my $ds_object = shift;
    my $data_options_ref = $ds_object->{options};
    $ds_object->set_type("BAR");
    $data_options_ref->{options}->type("BAR");
    $data_options_ref->{options}->symbol->fill_pattern("1");
    $data_options_ref->{options}->symbol->color("1");
    $data_options_ref->{options}->line->type("0");
    return 1;
}                                       

#
# Name: _set_XY
#
# Description: sets up the necessary characteristics of an XY
#              dataset type
#

sub _set_XY ($ ) {
    my $ds_object = shift;
    my $data_options_ref = $ds_object->{options};
    $ds_object->set_type("XY");
    $data_options_ref->{options}->type("XY");
    $data_options_ref->{options}->symbol->fill_pattern("0");
    $data_options_ref->{options}->line->type("1");
    return 1;
}                  

# Name: _set_title
#
# Description: sets up the necessary characteristics for the
#              title of the graph
#

sub _set_title ($ ) {
    my $ds_ref = shift;
    my $title = $ds_ref->{title};
    $ds_ref->{options}->{options}->{comment} = $title;
    $ds_ref->{options}->{options}->{legend} = $title;
    return 1;
}                                  

####################################################################
#
# The following functions are used to map an option value from 
# prose to it's corresponding number value which grace understands
# Note: these routines are not used in the current version. the user
#       must use numbers (NOT prose) for now.
#
###################################################################


# 
# Subroutine: _get_color()
# 
# Description: simply looks up a hash table and returns
#              a number that xmgrace understands that 
#              corresponds to the color
#

sub _get_color ($ ) {
    my $color = shift;
    my %colortable = (
		      "white" => "0",
		      "black" => "1",
		      "red" => "2",
		      "green" => "3",
		      "blue" => "4",
		      "yellow" => "5",
		      "brown" => "6",
		      "grey" => "7",
		      "violet" => "8",
		      "cyan" => "9",
		      "magenta" => "10",
		      "orange" => "11",
		      "indigo" => "12",
		      "maroon" => "13",
		      "turquoise" => "14",
		      "green4" => "15",
		      );
    
    $color = lc($color);	
    # "\L$color\E"
    #color =~ y/A-Z/a-z/;
    my $retval = $colortable{$color};
    
    if ($retval) {
	return $retval;
    } else {
	return 0;
    }
}

# 
# Subroutine: _get_symbol()
# 
# Description: simply looks up a hash table and returns
#              a number that xmgrace understands that 
#              corresponds to the symbol
#

sub _get_symbol ($ ) {
    my $symbol = shift;
    my %symtab = (
		  "none" => "0",
		  "circle" => "1",
		  "square" => "2",
		  "diamond" => "3",
		  "triangle up" => "4",
		  "triangle left" => "5",
		  "triangle down" => "6",
		  "triangle right" => "7",
		  "plus" => "8",
		  "x" => "9",
		  "star" => "10",
		  "char" => "11",
		 );
    
    $symbol = lc($symbol);  
    my $retval = $symtab{$symbol};
    
    if ($retval) {
	return $retval;
    } else {
	return 0;
    }    
}


# 
# Subroutine: _get_linetype()
# 
# Description: simply looks up a hash table and returns
#              a number that xmgrace understands that 
#              corresponds to the linetype
#

sub _get_linetype ($ ) {
    my $linetype = shift;
    my %linetypetab = (
		  "none" => "0",
		  "straight" => "1",
		  "left stairs" => "2",
		  "right stairs" => "3",
		  "segments" => "4",
		  "3-segments" => "5",
		 );
    
    $linetype = lc($linetype);  
    my $retval = $linetypetab{$linetype};
    
    if ($retval) {
	return $retval;
    } else {
	return 0;
    }    
}

# 
# Subroutine: _get_linestyle()
# 
# Description: simply looks up a hash table and returns
#              a number that xmgrace understands that 
#              corresponds to the linestyle
#

sub _get_linestyle ($ ) {
    my $linestyle = shift;
    my %linestyletab = (
		  "none" => "0",
		  "solid" => "1",
		  "dotted" => "2",
		  "en dash" => "3",
		  "em dash" => "4",
		  "dot-en dash" => "5",
		  "dot-em dash" => "6",
		  "dot-en-dot dash" => "7",
		  "en-dot-en dash" => "8",
		 );
    
    $linestyle = lc($linestyle);  
    my $retval = $linestyletab{$linestyle};
    
    if ($retval) {
	return $retval;
    } else {
	return 0;
    }    
}

# 
# Subroutine: _get_baselinetype()
# 
# Description: simply looks up a hash table and returns
#              a number that xmgrace understands that 
#              corresponds to the baselinetype
#

sub _get_baselinetype ($ ) {
    my $baselinetype = shift;
    my %baselinetypetab = (
		  "zero" => "0",
		  "set min" => "1",
		  "set max" => "2",
		  "graph min" => "3",
		  "graph max" => "4",
		 );
    
    $baselinetype = lc($baselinetype);  
    my $retval = $baselinetypetab{$baselinetype};
    
    if ($retval) {
	return $retval;
    } else {
	return 0;
    }    
}


# 
# Subroutine: _get_filltype()
# 
# Description: simply looks up a hash table and returns
#              a number that xmgrace understands that 
#              corresponds to the filltype
#

sub _get_filltype ($ ) {
    my $filltype = shift;
    my %filltypetab = (
		  "none" => "0",
		  "as polygon" => "1",
		  "to baseline" => "2",
		 );
    
    $filltype = lc($filltype);  
    my $retval = $filltypetab{$filltype};
    
    if ($retval) {
	return $retval;
    } else {
	return 0;
    }    
}

# 
# Subroutine: _get_avaluetype()
# 
# Description: simply looks up a hash table and returns
#              a number that xmgrace understands that 
#              corresponds to the avaluetype
#

sub _get_avaluetype ($ ) {
    my $avaluetype = shift;
    my %avaluetypetab = (
		  "none" => "0",
		  "x" => "1",
		  "y" => "2",
		  "xy" => "3",
		  "string" => "4",
		  "z" => "5",
	 );
    
    $avaluetype = lc($avaluetype);  
    my $retval = $avaluetypetab{$avaluetype};
    
    if ($retval) {
	return $retval;
    } else {
	return 0;
    }    
}


1;

__END__

=head1 NAME

Chart::Graph::Xmgrace

=head1 SYNOPSIS

 use Chart::Graph::Xmgrace qw(xmgrace);
 xmgrace(\%global_options, [\%data_set_options, \@matrix],
                           [\%data_set_options, \@x_column, \@y_column],
                           [\%data_set_options, < filename >], ... );

=head1 DESCRIPTION

The function xmgrace() is part of the module Chart::Graph that lets
you generate graphs on the fly in perl. It was written as a front-end
application to Xmgrace for hassle-free generation of graphs. xmgrace()
can be supplied with many of the same options and arguments that can
be given to Xmgrace (the UNIX program that evolved from xmgr). For
more information on Xmgrace see the end of this documentation.


=head1 ARGUMENTS

Xmgrace has a great deal of options for the overall appearance of a
graph.  Chart::Graph::Xmgrace provides control over an essential
subset of them.  Others can be easily changed by saving the file using
the C<grace output file> option and then manipulating the file
directing in Xmgrace.

 +----------------------------------------------------------------------------+
 |                             GLOBAL OPTIONS:                                |
 +-------------------+--------------------------------+-----------------------+
 |     NAME          |         OPTIONS                |        DEFAULT        |
 +-------------------+--------------------------------+-----------------------+
 |"title"            |   (set your own title)         | "untitled"            |
 |"subtitle"         |   (set your own subtitle)      | "untitled"            |
 |"type of graph"    |   "XY chart", "XY graph",      | "XY graph"            |
 |                   |   "Bar chart", "Bar graph"     |                       |
 |"output type"      |   "png"                        | "png"                 |
 |"output file"      |   (set your own output file)   | "untitled-grace.png"  |
 |"grace output file"|   (set your own grace output   | "untitled-grace.agr"  |
 |                   |    file)                       |                       |
 |"x-axis label"     |   (set your own label)         | "x-axis"              |
 |"y-axis label"     |   (set your own label)         | "y-axis"              |
 |"x2-axis label"    |   (set your own label)         | undefined             |
 |"y2-axis label"    |   (set your own label)         | undefined             |
 |"logscale x"       |   "0" or "1"                   | undefined             |
 |"logscale y"       |   "0" or "1"                   | undefined             |
 |"xtics"            |   (set your own ticks) look at | undefined             |
 |                   |    example                     |                       |
 |"ytics"            |   (set your own ticks) look at | undefined             |
 |                   |    example                     |                       |
 +-------------------+--------------------------------+-----------------------+

In Xmgrace each set of data has it's own options.  Because Xmgrace is
so complex.  a sub-hash of options is needed for all the options
associated with each data set.  For that reason, only a few data
options are noted here.

 +----------------------------------------------------------------------------+
 |                           DATA SET OPTIONS:                                |
 +-------------------+--------------------------------+-----------------------+
 |     NAME          |         OPTIONS                |        DEFAULT        |
 +-------------------+--------------------------------+-----------------------+
 |"title"            |    (set your own title)        | ""                    |
 |"options"          |    \%sub_options               |(depends on graph type)|
 |"data format"      |    "matrix","columns","file"   | ""                    |
 |"hidden"           |    "true" or "false"           | "false"               |
 +-------------------+--------------------------------+-----------------------+

Data can be presented to Chart::Graph::Xmgrace Gnuplot in one of three
formats for the convenience of the user:

 \@matrix: an array reference of [x,y] pairs of data

Alternatively:

 \@x_column, \@y_column: two array references of data of equal length.
 \@x_column is the x-axis data. \@y_column is the y-axis data.

Finally, data can be stored in a file as a parameter to be read into
C<Chart::Graph::Xmgrace>.

since xmgrace allows for many data set options, options is a hash of
suboptions (%sub_options below).

    %sub_options = (
                     "symbol" => \%symbol_options,
                     "line" => \%line_options,
                     "baseline" => \%baseline_options,
                     "dropline" => \%dropline_options,
                     "fill" => \%fill_options,
                     "avalue" => \%avalue_options
                     "errorbar" => \%errorbar_options,
                   );

There are seven types of suboptions as listed below and described in detail in the following tables.

=over 1

=item *

symbol options

=item *

line options

=item *

baseline options

=item *

dropline options

=item *

fill options

=item *

annotated value options

=item *

errorbar options

=back


 +----------------------------------------------------------------------------+
 |                           SYMBOL SUBOPTIONS:                               |
 +-------------------+--------------------------------+-----------------------+
 |     NAME          |         OPTIONS                |        DEFAULT        |
 +-------------------+--------------------------------+-----------------------+
 |  "type"           |      "0"..."11"                |    "none"             |
 |                   |      (look at symbol table)    |                       |
 |  "size"           |      (set own size)            |    "1.000000"         |
 |  "color"          |      (look at color table)     |    "auto"             |
 |  "pattern"        |      "0"..."31"                |    "1"                |
 |  "fill color"     |      "0"..."31"                |    "1"                |
 |  "fill pattern"   |      "0"..."31"                |    "1"                |
 |  "linewidth"      |      (set own linewidth)       |    "1.0" (max. value) |
 |  "linestyle"      |      "0"..."9"                 |    "1"                |
 |  "symbol char"    |      (not implemented)         |    "65"               |
 |  "char font"      |      (not implemented)         |    "0"                |
 |  "skip"           |      "0" or "1"                |    "0"                |
 +-------------------+--------------------------------+-----------------------+

 +----------------------------------------------------------------------------+
 |                            LINE SUBOPTIONS:                                |
 +-------------------+--------------------------------+-----------------------+
 |     NAME          |         OPTIONS                |        DEFAULT        |
 +-------------------+--------------------------------+-----------------------+
 |  "type"           |      (look at line type)       |    "1"                |
 |  "linestyle"      |      (look at line style)      |    "1"                |
 |  "linewidth"      |      (set own width)           |    "1.0" (max. value) |
 |  "color"          |      (look at color table)     |    "auto"             |
 |  "pattern"        |      "0"..."31"                |    "1"                |
 +-------------------+--------------------------------+-----------------------+

 +----------------------------------------------------------------------------+
 |                        BASELINE SUBOPTIONS:                                |
 +-------------------+--------------------------------+-----------------------+
 |     NAME          |         OPTIONS                |        DEFAULT        |
 +-------------------+--------------------------------+-----------------------+
 |  "type"           |      (look at baseline table)  |     "0"               |
 |  "status"         |      "on" or "off"             |     "off"             |
 +----------------------------------------------------------------------------+

 +----------------------------------------------------------------------------+
 |                        DROPLINE SUBOPTIONS:                                |
 +-------------------+--------------------------------+-----------------------+
 |     NAME          |         OPTIONS                |        DEFAULT        |
 +-------------------+--------------------------------+-----------------------+
 |  "status"         |      "on" or "off"             |      "off"            |
 +-------------------+--------------------------------+-----------------------+

 +----------------------------------------------------------------------------+
 |                            FILL SUBOPTIONS:                                |
 +-------------------+--------------------------------+-----------------------+
 |     NAME          |         OPTIONS                |        DEFAULT        |
 +-------------------+--------------------------------+-----------------------+
 |  "type"           | "as polygon" or "to baseline"  |      "as polygon"     |
 |  "rule"           | "winding" or "even-odd"        |      "winding"        |
 |  "color"          | (look at color table)          |      "auto"           |
 |  "pattern"        | "0"..."31"                     |      "1"              |
 +----------------------------------------------------------------------------+

 +----------------------------------------------------------------------------+
 |                          AVALUE SUBOPTIONS:                                |
 +-------------------+--------------------------------+-----------------------+
 |     NAME          |         OPTIONS                |        DEFAULT        |
 +-------------------+--------------------------------+-----------------------+
 |  "status"         |  "on" or "off"                 |      "off"            |
 |  "type"           |  "X","Y","XY","string","Z"     |      "XY"             |
 |  "char size"      |  (set your own size)           |      "1.000000"       |
 |  "font"           |  "0".."13"                     |      "0"              |
 |  "color"          |  (look at color table)         |      "auto"           |
 |  "rot"            |  (set own angle)               |      "0"              |
 |  "format"         |  "0"..."31"                    |      "1"              |
 |  "prec"           |  "0"..."9"                     |      "3"              |
 |  "prepend"        |  (set your own prepend)        |      ""               |
 |  "append"         |  (set your own apppend)        |      ""               |
 |  "offset"         |  ["own value", "own value"]    |      "[0.00, 0.00]"   |
 +----------------------------------------------------------------------------+

 +----------------------------------------------------------------------------+
 |                          ERRORBAR SUBOPTIONS:                              |
 +-------------------+--------------------------------+-----------------------+
 |     NAME          |         OPTIONS                |        DEFAULT        |
 +-------------------+--------------------------------+-----------------------+
 |  "status"         |  "on" or "off"                 |  "off"                |
 |  "place"          |  "normal","opposite","both"    |  "normal"             |
 |  "color"          |  (look at color table)         |  "auto"               |
 |  "pattern"        |  "0"..."31"                    |  "1"                  |
 |  "size"           |  (set your own size)           |  "1.000000"           |
 |  "font"           |  "0".."13"                     |  "0"                  |
 |  "linewidth"      |  (set own width)               |  "1.0" (max. value)   |
 |  "linestyle"      |  (look at line type)           |  "1"                  |
 |  "riser linewidth"|  (set own riser linewidth)     |  "1.0"                |
 |  "riser linestyle"|  (look at line type)           |  "1"                  |
 |"riser clip status"|  "on" or "off"                 |  "off"                |
 |"riser clip length"|  (set own clip length)         |  "0.100000"           |
 +----------------------------------------------------------------------------+

The suboptions above use the arguments listed below.

 +----------------------------------------------------------------------------+
 |                              SYMBOL TYPE:                                  |
 +--------+-------+--------+------+-------+--------+--------------------------+
 | SYMBOL | VALUE | SYMBOL | TYPE | VALUE | SYMBOL | VALUE                    |
 +--------+-------+--------+------+-------+--------+--------------------------+
 |  none  |  "0"  |triangle|  up  |  "4"  |  plus  |  "8"                     |
 | circle |  "1"  |triangle| left |  "5"  |   x    |  "9"                     |
 | square |  "2"  |triangle| down |  "6"  |  star  |  "10"                    |
 | diamond|  "3"  |triangle| right|  "7"  |  char  |  "11"                    |
 +--------+-------+--------+------+-------+--------+--------------------------+

 +-----------------------------------------------------------------+
 |                             LINE TYPE                           |
 +------------------------+-------+------------------------+-------+
 |  LINE TYPE             | VALUE |          LINE TYPE     | VALUE |
 +------------------------+-------+------------------------+-------+
 |     none               |  "0"  |          right stairs  |  "3"  |
 |   straight             |  "1"  |            segments    |  "4"  |
 | left stairs            |  "2"  |           3-segments   |  "5"  |
 +------------------------+-------+------------------------+-------+

 +-----------------------------------------------------------------+
 |                             LINE STYLE                          |
 +------------------------+-------+------------------------+-------+
 |  LINE STYLE            | VALUE |         LINE STYLE     | VALUE |
 +------------------------+-------+------------------------+-------+
 |       none             |  "0"  |     solid              |  "1"  |
 |       dotted           |  "2"  |     en-dash            |  "3"  |
 |       em-dash          |  "4"  |     dot-en dash        |  "5"  |
 |       dot-em dash      |  "6"  |     dot-en-dot dash    |  "7"  |
 |       en-dot-en dash   |  "8"  |                        |       |
 +------------------------+-------+------------------------+-------+

 +-----------------------------------------------------------------+
 |                             COLORS                              |
 +-------+-----+-------+-----+--------+-----+-----------+----------+
 | COLOR |VALUE| COLOR |VALUE| COLOR  |VALUE| COLOR     |  VALUE   |
 | white | "0" | blue  | "4" | violet | "8" | indigo    |   "12"   |
 | black | "1" | yellow| "5" | cyan   | "9" | maroon    |   "13"   |
 | red   | "2" | brown | "6" | magenta| "10"| turquoise |   "14"   |
 | green | "3" | grey  | "7" | orange | "11"| dark green|   "15"   |
 +-------+-----+-------+-----+--------+-----+-----------+----------+

=head1 EXAMPLES

The following three examples show the various capabilities of the
Chart::Graph interface to the Xmgrace program.

=head2 GENERAL EXAMPLE

The following example produces the file F<xmgrace1.png> and contains
three kinds of data plots.  The first plot is an XY plot using
triangles for the presentation style and rightstairs lines.  The
second plot is also an XY plot using lines andtriangle symbols.  The
last plot is a bar graph.

  # Include modules
  use Chart::Graph::Xmgrace qw(xmgrace);

  xmgrace( { "title" => "Example of a XY Chart",
             "subtitle" =>"optional subtitle",
             "type of graph" => "XY chart",
             "output type" => "png",
             "output file" => "xmgrace1.png",
             "x-axis label" => "my x-axis label",
             "y-axis label" => "my y-axis label",
             "logscale y" => "1",
             "xtics" => [ ["one", "1"], ["two", "2"], ["three", "3"] ],
             "ytics" => [ ["one", "1"], ["two", "2"], ["three", "3"] ],
             "grace output file" => "xmgrace1.agr",
           },

           [ { "title" => "XY presentation data1",
               "set presentation" => "XY",
               "options" => {
                           "line" => {
                                      "type" => "1",
                                      "color" => "8",
                                      "linewidth" => "1",
                                      "linestyle" => "3",
                                     },
                           "symbol" => {
                                        "symbol type" => "6",
                                        "color" => "1",
                                        "fill pattern" => "1",
                                        "fill color" => "1",
                                       },
                           "fill" => {
                                      "type" => "0",
                                     },
                          },
               "data format" => "matrix",
             },

             [ [1,2],
               [2,4],
               [3,6],
               [4,8],
               [5,10],
               [6,12],
               [7,14],
               [8,16],
               [9,18],
               [10,20] ]
           ],

           [ { "title" => "XY presentation data2",
               "options" => {
                           "line" => {
                                      "type" => "2",
                                      "color" => "4",
                                     },
                           "symbol" => {
                                        "symbol type" => "1",
                                        "color" => "1",
                                        "fill pattern" => "3",
                                        "fill color" => "5",
                                       },
                           "fill" => {
                                      "type" => "0",
                                     }
                          },
               "data format" => "columns",
             },
	     [
              [1,2,3,4,5,6,7,8,9,10],
              [3,6,9,12,15,18,21,24,27,30],
	     ]  
           ],

           [ { "title" => "BAR presentation data3",
               "set presentation" => "BAR",
               "data format" => "file"}, "sample"],

       );

=for html
<p><center><img src="http://www.caida.org/tools/utilities/graphing/xmgrace1.png"></center></p>
<p><center><em>xmgrace1.png</em></center></p>

=head2 NON-STACKING REGIONS

The following shorter example shows how Xmgrace handles regions
without stacking the graphs (the default for Xmgrace is to not stack
data.)

  # Include modules
  use Chart::Graph::Xmgrace qw(xmgrace);

  xmgrace({"title" => "Example of a XY graph",
           "subtitle" => "optional subtitle",
           "type of graph" => "XY graph",
           "output type" => "png",
           "output file" => "xmgrace2.png",
	   "grace output file" => "xmgrace2.agr",
           "x-axis label" => "my x-axis label",
           "y-axis label" => "my y-axis label"
	  },
	  [{"title" => "data",
	    "options" => {
                          "fill" => { "type" => "2" },
			 },
            "data format" => "file"
	   },
	   "sample"
	  ],
	 );

=for html
<p><center><img src="http://www.caida.org/tools/utilities/graphing/xmgrace2.png"></center></p>
<p><center><em>xmgrace2.png</em></center></p>


=head2 MULTIPLE DATA SETS IN MATRIX FORM

The following example shows how to graph more complicated datasets
using the Chart-Graph interface to Xmgrace.  It produces the file
F<xmgrace3.png>.The numbers from this example were generated from the
script that created it and saved using the standard Perl module
Data-Dumper.

  # Include modules
  use Chart::Graph::Xmgrace qw(xmgrace);

	xmgrace({'y-axis label' => 'Percent of widgets',
		 'output file' => 'xmgrace3.png',
		 'type of graph' => 'Bar chart',
		 'output type' => 'png',
		 'title' => 'Percent of widgets',
		 'grace output file' => 'xmgrace3.agr',
		 'subtitle' => 'Data collected from 07/24/2001 to 08/01/2001',
		 'x-axis label' => 'Date of data sampling'
		},
		[{'data format' => 'matrix',
		  'title' => 'Widget A'
		 },
		 [
		  [ '2001-07-24',  '32.58' ],
		  [ '2001-07-25',  '30.4291287386216'  ],
		  [ '2001-07-26',  '34.4106463878327'  ],
		  [ '2001-07-27',  '34.44'	  ],
		  [ '2001-07-28',  '37.4482270936458' ],
		  [ '2001-07-29',  '37.8769479862376'  ],
		  [ '2001-07-30',  '34.9437860832574'  ],
		  [ '2001-07-31',  '36.0707388962293'  ],
		  [ '2001-08-01',  '40.0591353996737'  ]
		 ]
		],
		[{'data format' => 'matrix',
		  'title' => 'Widget B'
		 },
		 [
		  [ '2001-07-24',  '29.13'  ],
		  [ '2001-07-25',  '30.8192457737321'  ],
		  [ '2001-07-26',  '29.1775065039023'  ],
		  [ '2001-07-27',  '29.82'             ],
		  [ '2001-07-28',  '28.9221133447823'  ],
		  [ '2001-07-29',  '28.5772110908723'  ],
		  [ '2001-07-30',  '29.2109794388737'  ],
		  [ '2001-07-31',  '26.8624860250025'  ],
		  [ '2001-08-01',  '8.442088091354'    ]
		 ]
		],
		[
		 {
		  'data format' => 'matrix',
		  'title' => 'Widget C'
		 },
		 [
		  [ '2001-07-24', '15.42'        ],
		  [ '2001-07-25', '17.2251675502651' ],
		  [ '2001-07-26', '15.6093656193716' ],
		  [ '2001-07-27', '16.02'            ],
		  [ '2001-07-28', '14.526719870694'  ],
		  [ '2001-07-29', '15.1791135397693' ],
		  [ '2001-07-30', '16.8337891218475' ],
		  [ '2001-07-31', '16.3227970322187' ],
		  [ '2001-08-01', '17.7304241435563' ]
		 ]
		],
		[
		 {
		  'data format' => 'matrix',
		  'title' => 'Widget D'
		 },
		 [
		  [ '2001-07-24', '7.61'  ],
		  [ '2001-07-25', '7.80234070221066' ],
		  [ '2001-07-26', '7.82469481689013' ],
		  [ '2001-07-27', '7.57'            ],
		  [ '2001-07-28', '7.72805333872108'  ],
		  [ '2001-07-29', '7.34669095324833' ],
		  [ '2001-07-30', '7.95097741314697' ],
		  [ '2001-07-31', '10.7226344140665'  ],
		  [ '2001-08-01', '12.9282218597064'  ]
		 ]
		],
		[
		 {
		  'data format' => 'matrix',
		  'title' => 'Widget E'
		 },
		 [
		  [  '2001-07-24', '10.75'  ],
		  [  '2001-07-25', '9.53285985795739'  ],
		  [  '2001-07-26', '8.375025015009'    ],
		  [  '2001-07-27', '7.79'           ],
		  [  '2001-07-28', '6.32387109809072'  ],
		  [  '2001-07-29', '6.90143695608177'  ],
		  [  '2001-07-30', '6.26962422769169'  ],
		  [  '2001-07-31', '5.43754446590101'  ],
		  [  '2001-08-01', '14.8960032626427'  ]
		 ]
		],
		[
		 {
		  'data format' => 'matrix',
		  'title' => 'Widget F'
		 },
		 [
		  [  '2001-07-24', '3.16'         ],
		  [  '2001-07-25', '2.68080424127238'   ],
		  [  '2001-07-26', '3.08184910946568'   ],
		  [  '2001-07-27', '2.85'           ],
		  [  '2001-07-28', '2.78816042024447'  ],
		  [  '2001-07-29', '2.6006881198138'   ],
		  [  '2001-07-30', '3.0892332624329'   ],
		  [  '2001-07-31', '3.02876308567944'  ],
		  [  '2001-08-01', '3.02814029363785'  ]
		 ]
		],
		[
		 {
		  'data format' => 'matrix',
		  'title' => 'Widget G'
		 },
		 [
		  [ '2001-07-24',  '1.14'      ],
		  [ '2001-07-25',  '1.28038411523457'  ],
		  [ '2001-07-26',  '1.26075645387232'  ],
		  [ '2001-07-27',  '1.33'              ],
		  [ '2001-07-28',  '2.09112031518335'  ],
		  [ '2001-07-29',  '1.27504553734062'  ],
		  [ '2001-07-30',  '1.43826597791958'  ],
		  [ '2001-07-31',  '1.31110885252566'  ],
		  [ '2001-08-01',  '2.76305057096248'  ]
		 ]
		],
		[
		 {
		  'data format' => 'matrix',
		  'title' => 'Widget H'
		 },
		 [
		  [ '2001-07-24', '0.09'	  ],
		  [ '2001-07-25', '0.110033009902971'  ],
		  [ '2001-07-26', '0.150090054032419'  ],
		  [ '2001-07-27', '0.07'             ],
		  [ '2001-07-28', '0.111122335589453' ],
		  [ '2001-07-29', '0.121432908318154' ],
		  [ '2001-07-30', '0.121543603767852' ],
		  [ '2001-07-31', '0.111799979672731' ],
		  [ '2001-08-01', '0.0815660685154976']
		 ]
		],
		[
		 {
		  'data format' => 'matrix',
		  'title' => 'Widget I'
		 },
		 [
		  [  '2001-07-24', '0.04'  ],
		  [  '2001-07-25', '0.0500150045013504'  ],
		  [  '2001-07-26', '0.0500300180108065'  ],
		  [  '2001-07-27', '0.02'             ],
		  [  '2001-07-28', '0.0303060915243964' ],
		  [  '2001-07-29', '0.0607164541590771'  ],
		  [  '2001-07-30', '0.0709004355312468'  ],
		  [  '2001-07-31', '0.0203272690314056'  ],
		  [  '2001-08-01', '0.0101957585644372'  ]
		 ]
		],
		[
		 {
		  'data format' => 'matrix',
		  'title' => 'Widget J'
		 },
		 [
		  [ '2001-07-24', '0.03'  ],
		  [ '2001-07-25', '0.0600180054016205'  ],
		  [ '2001-07-26', '0.0400240144086452'  ],
		  [ '2001-07-27', '0.08' ],
		  [ '2001-07-28', '0.0202040610162643'   ],
		  [ '2001-07-29', '0.0303582270795386'   ],
		  [ '2001-07-30', '0.0607718018839259'   ],
		  [ '2001-07-31', '0.0609818070942169'   ],
		  [ '2001-08-01', '0.0203915171288744'   ]
		 ]
		],
		[
		 {
		  'data format' => 'matrix',
		  'title' => 'Widget K'
		 },
		 [
		  [ '2001-07-24', '0.05' ],
		  [ '2001-07-25','0.0100030009002701' ],
		  [ '2001-07-26','0.0200120072043226' ],
		  [ '2001-07-27', '0.01'             ],
		  [ '2001-07-28','0.0101020305081321' ],
		  [ '2001-07-29', '0.0303582270795386' ],
		  [ '2001-07-30',  '0.010128633647321'  ],
		  [ '2001-07-31',  '0.0508181725785141' ],
		  [ '2001-08-01',  '0.0407830342577488' ]
		 ]
		]
	       ) # xmgrace call


=for html
<p><center><img src="http://www.caida.org/tools/utilities/graphing/xmgrace3.png"></center></p>
<p><center><em>xmgrace3.png</em></center></p>



=head1 MORE INFO

For more information on Xmgrace, please see the Xmgrace web page:

 http://plasma-gate.weizmann.ac.il/Grace


=head1 CONTACT

Send email to graph-dev@caida.org is you have problems, questions,
or comments. To subscribe to the mailing list send mail to
graph-dev-request@caida.org with a body of "subscribe your@email.com"

=head1 AUTHOR

 CAIDA Perl development team (cpan@caida.org)

=head1 SEE ALSO

 xmgrace(1).


=cut

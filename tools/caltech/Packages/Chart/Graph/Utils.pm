## Graph.pm is a graphing package that supports on-the-fly graphing 
## from the gnuplot, xrt, and xmgrace  graphing packages.
##
## $Id: Utils.pm,v 1.24 2006/06/07 21:09:33 emile Exp $ $Name:  $
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
package Chart::Graph::Utils;
use Exporter ();

@ISA = qw(Exporter);
%EXPORT_TAGS = (UTILS => [qw($gnuplot $xmgrace $convert $ppmtogif 
                             $xwdtopnm $xrt $xvfb $ppmtojpg $pnmtopng
			     $use_xvfb 
			     $tmpcount 
			     &_make_tmpdir &_cleanup_tmpdir &_get_path 
			     &_chk_status &_mesh_opts &_make_tmpfile 
                             &_number_to_eng)],

                 	     # variables that user may set
		USER => [qw($gnuplot $convert $xmgrace $ppmtogif $pnmtopng 
                    $ppmtojpg $xwdtopnm $xrt $xvfb 
			    $use_xvfb)]
	       );

# add symbols from tags into @EXPORT_OK
Exporter::export_ok_tags('UTILS');

use Carp;			# for carp() and croak()
use File::Path;	                # for rmtree()

$cvs_Id = '$Id: Utils.pm,v 1.24 2006/06/07 21:09:33 emile Exp $';
$cvs_Author = '$Author: emile $';
$cvs_Name = '$Name:  $';
$cvs_Revision = '$Revision: 1.24 $';

$VERSION = 3.2;

use strict;

# default don't use xvfb by default, use their own display environment
use vars qw($use_xvfb);
$use_xvfb = 0;

# hold paths to programs 
# user may choose to set paths to these programs
# if paths are not set then we will attempt to search
# PATH for the programs                                    
use vars qw($gnuplot $xvfb $convert $xwdtopnm $ppmtogif $ppmtojpg $pnmtopng); 
# Store names of converters programs until either user sets path or
# path is automatically assigned to these variables.  These will be overridden
# by users if another path is desired.
$convert  = "convert";
$xwdtopnm = "xwdtopnm";   
$ppmtogif = "ppmtogif";
$ppmtojpg = "ppmtojpeg";
$pnmtopng = "pnmtopng";

#
# remove tmp files in case program exits abnormally
#
END {
    _cleanup_tmpdir();
}

#
# remove tmp files in case of a signal interrupt (ctrl-C)
#
$SIG{INT} = \&_handle_sigint;
#$SIG{INT} = 'DEFAULT';

#
# general purpose global variables
#								

use vars qw($tmpcount $tmpdir @clean_tmpdirs);
$tmpcount = 0;	# used to create unique tmp filenames
@clean_tmpdirs = (); # Storage for every tmp directory to be eventually 
                     # removed.

#
#
# general purpose subroutines - these are subroutines shared across
#                               all packages     
#

# 
# Subroutine: make_tmpdir()
# 
# Description: creates temporary dir for storage 
#              of temporary files with read, write,
#              and execute for user and group
sub _make_tmpdir {
    my $package = shift;
    die 'Too many arguments\n' if @_;

    my $PID = $$;
    if (not defined($ENV{TMPDIR})) {
	while(-d "/tmp/Graph$package$PID") {$PID++ }
	$tmpdir = "/tmp/Graph$package$PID";
    } else {
	while(-d "/tmp/Graph$package$PID") {$PID++ }
	$tmpdir = "$ENV{TMPDIR}/Graph$package$$";
    }

    if (not mkdir($tmpdir, 0770)) {
	croak "could not make temporary directory: `$tmpdir'";
	$tmpdir = undef;
    }

    # Add directory to list of directories to remove
    push @clean_tmpdirs, $tmpdir;
    return $tmpdir;
}

# 
# Subroutine: cleanup_tmpdir()
# 
# Description: remove the tmp dir we created for
#              tmp files

sub _cleanup_tmpdir {

    my $count = 1;     # Convenience variable for comments.

    # Loop through every temporary directory created and either
    # report to the user it was created or remove the tree structure.
    foreach my $tmp_dir (@clean_tmpdirs) {
	if ($Chart::Graph::save_tmpfiles) {
	    if (defined($tmp_dir) and -d $tmp_dir ) {
		print "Set $count of tmp files are in $tmp_dir\n";
	    }
	} elsif (defined($tmp_dir) and -d $tmp_dir ) {
	    # remove directory and associated files.
	    unless ( rmtree ($tmp_dir, 0, 1) ) {
		carp "Unable to successfully remove set $count of " .
		     "temporary files";
	    } # Careful not to change permissions but not unlink files.
	}
	$count++;
    }
}

# 
# Subroutine: make_tmpfile() 
# Description: create temporary filenames with unique extensions
#
#
sub _make_tmpfile {
    my ($file, $ext) = @_;
    $tmpcount++;
    if (not defined($ext)) {
	$ext = "";
    } elsif ($ext ne "") {
	$ext = ".$ext"
    };
    return "$tmpdir/$file.$tmpcount$ext";
}

# 
# Subroutine: get_path()
# 
# Description: searches PATH for specified program given as arg 
#
#
#
#
sub _get_path {
    my ($exe) = @_;
    my @path = split (/:/, $ENV{PATH});
    my $program;	
    
    foreach my $i(@path){
	  $program = "$i/$exe";
	  if (-x $program) {
		return $program;
	  }
    }

    carp "program not found in search path: $exe";
    return 0;
}

# 
# Subroutine: chk_status
# 
# Description: checks the exit status of system() calls for errors
#
#
sub _chk_status {
    my ($status) = @_;
    if ($status) { 
	my $exit_value = $? >> 8;
	my $signal_num = $? & 127;
	my $dumped_core = $? & 128;
	carp "exit value = $exit_value\nsignal number = $signal_num\ndumped core = $dumped_core\n";
	return 0;
    }
    return 1;
}

# 
# Subroutine: mesh_opts
#
# Description: merges user and default option for 
#              gnuplot and or xrt options
#

sub _mesh_opts {
    my ($user_opts_ref, $default_opts_ref) = @_;
    my %user_opts = %{$user_opts_ref};
    my %default_opts = %{$default_opts_ref};
    my %opts;

    # check user opts against defaults and mesh
    # the basic algorithm here is to override the 
    # the default options against the ones that
    # the user has passed in. 
    while (my ($key, $value) = each %default_opts) {
	if (exists($user_opts{$key})) {
	  $opts{$key} = $user_opts{$key};
	  delete $user_opts{$key}; # remove options 
	  # that are matching
	} else {
	  $opts{$key} = $default_opts{$key};
	}
  }
    
    # any left over options in the table are unknown
    # if the user passes in illegal options then we 
    # warn them with an error message but still 
    # proceed.
    while (my ($key, $value) = each %user_opts) {
	carp "unknown option: $key";
    }
    
    return %opts;
}


# 
# Subroutine: _number_to_eng
#
# Description: returns the correct suffix for  
#              numerics, ie, 1st, 2nd, 3rd, 4th..
#

sub _number_to_eng {
    my $number = shift @_;
    my $retval;

    if ($number == 1) {
	$retval = $number . "st";
    } elsif ($number == 2) {
	$retval = $number . "nd";
    } elsif ($number == 3) {
	$retval = $number . "rd";
    } else {
	$retval = $number . "th";
    }

    return $retval;
}

# 
# Subroutine: _handle_sigint
#
# Description: cleans up if ctrl c detected  
#             
#

sub _handle_sigint {
    warn "\nSIGINT detected: tmp files cleaned up\n";
    _cleanup_tmpdir();
    $SIG{INT} = 'DEFAULT';
    exit;
   
}

1;

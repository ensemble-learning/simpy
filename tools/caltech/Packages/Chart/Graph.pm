## Graph.pm is a graphing package that supports on-the-fly graphing 
## from the gnuplot, XRT, and Xmgrace graphing packages.
##
## $Id: Graph.pm,v 1.48 2006/06/07 21:09:32 emile Exp $ $Name:  $
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
package Chart::Graph;
use Exporter();

$VERSION = '3.2';
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(&gnuplot &xrt3d &xrt2d &xmgrace);

$cvs_Id = '$Id: Graph.pm,v 1.48 2006/06/07 21:09:32 emile Exp $';
$cvs_Author = '$Author: emile $';
$cvs_Name = '$Name:  $';
$cvs_Revision = '$Revision: 1.48 $';

use Chart::Graph::Utils qw(:USER);
use Chart::Graph::Gnuplot qw(&gnuplot);
use Chart::Graph::Xrt2d qw(&xrt2d);
use Chart::Graph::Xrt3d qw(&xrt3d);
use Chart::Graph::Xmgrace qw(&xmgrace);

# Create "global" variables for use in various parts of Chart-Graph.
# use vars qw($debug $save_tmpfiles);

# debug level, people can set this with $Chart::Graph::debug = 1; after the
# use Graph; to increase debug output.

$debug = 0;			# turn debug mesg off by default

# set this with $Chart::Graph::save_tmpfiles = 1; after use Graph; to keep
# temporary (command) files used by the graphing packages around

$save_tmpfiles = 0;		# turn tmpfiles off by default


1;
__END__

=head1 NAME

Chart::Graph - Perl extension for a front-end to gnuplot, XRT, and Xmgrace.

=head1 SYNOPSIS


 # EXAMPLE: gnuplot
 #make sure to include Chart::Graph
 use Chart::Graph:Gnuplot qw(gnuplot);

 gnuplot(\%global_options, [\%data_set_options, \@matrix],
                           [\%data_set_options, \@x_column, \@y_column],
                           [\%data_set_options, < filename >], ... );



 # EXAMPLE: Xmgrace
 #make sure to include Chart::Graph
 use Chart::Graph::Xmgrace qw(xmgrace);
 xmgrace(\%global_options, [\%data_set_options, \@matrix],
                           [\%data_set_options, \@x_column, \@y_column],
                           [\%data_set_options, < filename >], ... );

 # EXAMPLE: xrt2d
 #make sure to include Chart::Graph
 use Chart::Graph::Xrt2d qw(xrt2d);

 xrt2d(\%options, \@data_set);

 #say for example we have a 3 by 4 matrix -> dataxy
 xrt2d(\%options,
       [[data11, data12, data13, data14],
       [data21, data22, data23, data24],
       [data31, data32, data33, data34]])


 # EXAMPLE: xrt3d
 #make sure to include Chart::Graph
 use Chart::Graph::Xrt3d qw(xrt3d);

 xrt3d(\%options, \@data_set);

 #say for example we have a 3 by 4 matrix -> dataxy
 xrt3d(\%options,
       [[data11, data12, data13, data14],
       [data21, data22, data23, data24],
       [data31, data32, data33, data34]])


=head1 DESCRIPTION

 use Chart::Graph;

Graph.pm is a wrapper module that allows easy generation of graphs
within perl. Currently Graph.pm supports three graphing packages,
gnuplot, XRT, and Xmgrace.  These software packages must be obtained
separately from this Perl module.  Information on each graphing
package and it's availability is provided in the documentation on that
module.  Gnuplot and Xmgrace are freely available software pages for
UNIX systems.  XRT is a commercial product.

Currently the xrt3d and xrt2d package is not being supported,
although it works. It is still in the development stage. Feel free
to give it a try though.


=head1 INSTALLATION

Because Chart-Graph is a wrapper script, you need to install the
graphic package that you wish to use B<before> attempting to install
Chart-Graph.pm.  Unless the appropriate graphics software is
installed, the testing portions of the install will fail.

If you want to use the xrt2d/xrt3d package, you need to build the
respective "graph" binaries in the xrt2d/xrt3d directories.  Refer to
the F<README> files in the xrt2d and xrt3d directories for
instructions on creating the required binaries.

Chart-Graph.pm use the standard Perl module installation procedure.  To
install Graph.pm on your system you need to run:

         perl Makefile.PL
         make
         make test
         make install

The standard Perl options apply. For example you can specify the
location of you Perl Installation by the option:
C<PREFIX=/home/your/private/dir>.  Which results modifying the first
command as follows:

         perl Makefile.PL PREFIX=/your/private/dir

Running the Makefile.PL Perl script mades some additional preparation
before creating the the Makefile.  In particular, the script sets up
the testing routines for the various graphics modules.

        Enter (space separated) graphing drivers to test: [gnuplot xrt3d xrt2d xmgrace]

Enter the names of the graphical software packages that you have
installed.  The others will be ignored even if Chart-Graph will appear
to "test" them.  If you are using Xmgrace or XRT there are additional
options you will need to supply.  In order to permit Xmgrace to
perform its test without using the X server you should provide a path
to a X virtual frame buffer.

        Enter path to X virtual frame buffer(Xvfb):


Finally, If you are running XRT, you need to provide that path to the
XRT binaries:

        Enter path to xrt2d binary (built from xrt2d/):

Note that running the tests will create image files that are placed in
the directory F<test_results>.  These images are almost all identical
with the examples provided with the documentation and can be used to
check if there are subtile errors in your image creation software.

=head1 USAGE

Chart::Graph attempts as much as possible to provide a uniform
interface to these different graphics packages.  Unfortunately, the
functionality of each program is sufficiently different that the
interface cannot be entirely uniform.

=head2 GENERAL DIAGNOSTICS AND TOOLS

Currently Chart::Graph supports two levels of debug, C<0> (no debug msgs)
and C<1>(debug msgs). You need to set the C<$Chart::Graph::debug> flag
accordingly. If you are having problems with Graph.pm set the debug
flag to C<1>. Also Graph.pm will check C<$ENV{TMPDIR}> for the temporary
file storage. If you do not specify, it will be set to F</tmp>
automatically. Temporary files can also be saved for further debugging
purposes by setting C<$Chart::Graph::save_tmpfiles> flag accordingly, C<0>
(no tmp files saved) or C<1> (save tmp files specified in C<$ENV{TMPDIR)> or
/tmp by default)

Note: Currently, XRT and Xmgrace use the local x server to draw it's
graphics by default. With XRT, if you are having problems with color
or speed is an issue, set C<$Chart::Graph::use_xvbf> to C<1> to use the
virtual x frame buffer. With Xmgrace, you MUST set
C<$Chart::Graph::use_xvbf> to C<1> if you are not using a local x server.

All the documentation is also provided in HTML with the sample graphic
files for Graph.pm are located in the F<doc> directory.

=head1 CONTENT SUMMARY

 Graph.pm        - top level file of Chart::Graph
 Graph/          - sub modules of Chart::Graph
 Graph/Xmgrace/  - sub modules of Chart::Graph::Xmgrace
 doc/            - documentation in HTML 
 xrt2d/          - xrt2d wrapper executable code
 xrt3d/          - xrt3d wrapper executable code
 test_Graph.pl   - the test script used for debugging

=head1 MORE INFO

For more information on gnuplot, please see the gnuplot web page:

 http://www.gnuplot.org/

For more information on Xmgrace, please see the Xmgrace web page:

 http://plasma-gate.weizmann.ac.il/Grace

For more information on XRT, please contact Sitraka (now Quest):

 http://www.quest.com/xrt_pds/

=head1 CONTACT

Send email to graph-dev@caida.org is you have problems, questions,
or comments. To subscribe to the mailing list send mail to
graph-dev-request@caida.org with a body of "subscribe your@email.com"

=head1 AUTHOR

 CAIDA Perl development team (cpan@caida.org)

=head1 SEE ALSO

 gnuplot(1).
 xmgrace(1).
 

=cut

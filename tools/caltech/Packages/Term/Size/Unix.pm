package Term::Size::Unix;

use strict;
use Carp;
use vars qw(@EXPORT_OK @ISA $VERSION);

use AutoLoader ();
use DynaLoader ();
use Exporter ();

@ISA = qw(Exporter DynaLoader);
@EXPORT_OK = qw(chars pixels);

$VERSION = '0.204';

=head1 NAME

Term::Size::Unix - Perl extension for retrieving terminal size (Unix version)

=head1 SYNOPSIS

    use Term::Size::Unix;

    ($columns, $rows) = Term::Size::Unix::chars *STDOUT{IO};
    ($x, $y) = Term::Size::Unix::pixels;

=head1 DESCRIPTION

  THIS IS AN UNOFFICIAL PATCH AGAINST Term-Size 0.2 DISTRIBUTION 
  FOUND ON CPAN (http://search.cpan.org/~timpx/Term-Size-0.2/).
  IT IS UNOFFICIAL IN THE SENSE THAT THE AUTHOR Tim Goodwin 
  HASN'T APPROVED IT (YET, I HOPE). BECAUSE OF THIS, THIS 
  DISTRIBUTION IS NOT INDEXED AND AVAILABLE VIA cpan OR cpanp 
  SHELLS UNLESS YOU EXPLICITLY SAY 
  "install FERREIRA/Term-Size-0.203.tar.gz". 
  
  THIS IS UNDELICATE? I THINK IT IS IN A CERTAIN SENSE. BUT IT 
  IS A WAY TO UNFREEZE THE CURRENT DISTRIBUTION STATUS. IF TIM 
  DISAPPROVES, I WILL REMOVE THIS DISTRIBUTION RIGHT AWAY. 
  IF HE APPROVES, I WILL REMOVE THIS DISTRIBUTION RIGHT AWAY 
  AND WORK OUT (AFTER BEEN GIVEN MAINTAINERSHIP STATUS) 
  A DISTRIBUTION WITHOUT THIS NOTE AND WHICH INDEXES CORRECTLY.

B<Term::Size> is a Perl module which provides a straightforward way to
retrieve the terminal size.

Both functions take an optional filehandle argument, which defaults to
C<*STDIN{IO}>.  They both return a list of two values, which are the
current width and height, respectively, of the terminal associated with
the specified filehandle.

C<Term::Size::chars> returns the size in units of characters, whereas
C<Term::Size::pixels> uses units of pixels.

In a scalar context, both functions return the first element of the
list, that is, the terminal width.

The functions may be imported.

If you need to pass a filehandle to either of the C<Term::Size>
functions, beware that the C<*STDOUT{IO}> syntax is only supported in
Perl 5.004 and later.  If you have an earlier version of Perl, or are
interested in backwards compatibility, use C<*STDOUT> instead.

=head1 EXAMPLES

1. Refuse to run in a too narrow window.

    use Term::Size;

    die "Need 80 column screen" if Term::Size::chars *STDOUT{IO} < 80;

2. Track window size changes.

    use Term::Size 'chars';

    my $changed = 1;

    while (1) {
            local $SIG{'WINCH'} = sub { $changed = 1 };

            if ($changed) {
                    ($cols, $rows) = chars;
                    # Redraw, or whatever.
                    $changed = 0;
            }
    }

=head1 RETURN VALUES

Both functions return C<undef> if there is an error.

If the terminal size information is not available, the functions
will normally return C<(0, 0)>, but this depends on your system.  On
character only terminals, C<pixels> will normally return C<(0, 0)>.

=head1 BUGS

It only works on Unix systems.

=head1 AUTHOR

Tim Goodwin, <tim@uunet.pipex.com>, 1997-04-23.

Candidate for maintainership:
Adriano Ferreira, <ferreira@cpan.org>, 2006-05-19.

=cut

bootstrap Term::Size::Unix $VERSION;

1;

__END__

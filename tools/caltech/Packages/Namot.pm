package Packages::Namot;

BEGIN {
    unshift @INC, "/ul/tpascal/.libs/Namot2";
}

require Exporter;
use p5namot;

p5namot::Cmd("set hush ERROR off");
p5namot::Cmd("set hush INFO off");
p5namot::Cmd("set hush REQUESTED off");
p5namot::Cmd("set hush WARNING off");

our (@ISA, @EXPORT, $VERSION, @EXPORT_OK);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(DoNamotCmd); 

$VERSION = "1.00";

sub DoNamotCmd {
    my ($cmd_str) = $_[0];
    my ($i);

    p5namot::Cmd("set background black");
    
    for $i (1 .. $#{$cmd_str}) {
	p5namot::Cmd("$cmd_str->[$i]");
	p5namot::Cmd("render");
    }
    p5namot::Cmd("close");
}

1;


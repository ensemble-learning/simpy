package Packages::MovieGen;
BEGIN {
    push (@INC, "/ul/tpascal/.libs/Net/Ssh/lib/perl5/site_perl/5.6.1");
}
 
require Exporter;
use Net::SSH qw(ssh);
use File::Basename;

our (@ISA, @EXPORT, $VERSION, @EXPORT_OK);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(sshCmd CreateMovie);
$VERSION = "1.00";

sub sshCmd {
    my ($usrname, $host, $cmd) = @_;
    my ($result);

    $result = ssh("$usrname\@$host", "$cmd");
 
    $result ?
	$result = 0 :
        $result = 1;
    return $result;

}

sub CreateMovie {
    my ($indata) = @_;
    my ($result, $curr_file, $file_string, $cmd);

    $file_string = "";
    for($indata->{"start"} .. $indata->{"end"}) {
	$curr_file = $indata->{"filebase"} . $_ . "." . $indata->{"extension"};
	if (-e $curr_file) {
	    $file_string .= "$curr_file ";
	    $movie_file_string .= "/temp1/" . $indata->{"user"} . "/tempics/" . 
		basename($indata->{"filebase"}) . $_ . "." . 
		$indata->{"extension"} . " ";
	}
    }

    if ($indata->{"framerate"} < 2.5) {
	$indata->{"framerate"} = 2.5;
    }

    if ($file_string ne "") {
	$cmd = "mkdir -p /temp1/" . $indata->{"user"} . "/tempics";
	if (sshCmd($indata->{"user"}, $indata->{"host"}, $cmd)) {

	    system "scp $file_string " . $indata->{"host"} . ":/temp1/" . 
		$indata->{"user"} . "/tempics >& dumpfile.txt";
	    $cmd = "makemovie -o /temp1/" . $indata->{"user"};
	    $cmd .= "/tempics/" . $indata->{"savename"} . " -f qt -r ";
	    $cmd .=  $indata->{"framerate"} . " -c qt_video -s ";
	    $cmd .= $indata->{"moviesize"} . " $movie_file_string";
	    if (sshCmd($indata->{"user"}, $indata->{"host"}, $cmd)) {
		$cmd = "scp " . $indata->{"host"} . ":/temp1/" . $indata->{"user"};
		$cmd .= "/tempics/" . $indata->{"savename"} . " . >& dumpfile.txt";
		system $cmd;
		$cmd = "rm -f $movie_file_string /temp1/" . 
		    $indata->{"user"} . "/tempics/" . $indata->{"savename"};
		sshCmd($indata->{"user"}, $indata->{"host"}, $cmd);
		$result = 1;
	    } else {
		$result = 0;
	    }
	} else {
	    $result = 0;
	}
    } else {
	$result = 0;
    }
    
    return $result;
}

1;

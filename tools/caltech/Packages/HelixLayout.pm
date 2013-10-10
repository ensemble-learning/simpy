package Packages::HelixLayout;

require Exporter;

our(@EXPORT_OK, @ISA, @EXPORT, $VERSION, @helix, @layout_array);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(GetHelixInfo);
$VERSION = "1.0";



#Initilization of Class
sub spawn {
    my $invocant = shift;
    my $class = ref($invocant) || $invocant;
    my $self = { };
    bless($self, $class);
    return $self;
}

sub GetHelixInfo {

   return @helix;
}

sub DefineChainLayout {
# This subroutine will define the chain layout in a 4x4 array
# The offset will represent where the helices are arranged 5' 3' 3' 5'
# at the bottom or 3' 5' 5' 3'
# @chain_layout holds the initial chain layout with 1, 3 = 5', 2, 4 = 3'

    if ($_[0]) {
	$layout_array[0] = [ 1, 2, 4, 3 ];
        $layout_array[1] = [ 1, 2, 4, 3 ];
	$layout_array[2] = [ 1, 4, 2, 3 ];
	$layout_array[3] = [ 1, 4, 2, 3 ];
    }else { 
	$layout_array[0] = [ 1, 2, 4, 3 ];
	$layout_array[1] = [ 1, 4, 2, 3 ];
	$layout_array[2] = [ 1, 4, 2, 3 ];
	$layout_array[3] = [ 1, 2, 4, 3 ]; 
    }
}

sub  DetermineHelixLayout {
# helix[helix_no][strand_no][region_no]
# Flip the order of the bases for the second strand
# Added offset too handle missing crossovers

    my ($invocant) = shift;
    my ($i, $start, $end, $curr_region, $cross_counter, @Regions);
    my ($h1_s1, $h2_s2, $h1_s2, $h2_s1, $was_previous_cross, $j, $offset);

    my ($maj_groove, $min_groove, $is3Prime, $numAtEnd, $totalBases, @crossovers) = @_;

    @Regions = GetRegions($maj_groove, $min_groove, $totalBases, $numAtEnd);

    if ($#Regions > $#crossovers) {
	$offset = ($#Regions - $#crossovers) - 1;
    } else {
	$offset = 0;
    } 

    $curr_region = (($#crossovers + 1 - $offset) % 4);
    
    if ($offset == 2 && $curr_region == 2) {
	$curr_region = 0;
    }

#    print "OFF: $offset REGIONS: $#Regions CROSS: $#crossovers ", 
#    "CURR REGION: $curr_region IS3PRIME: $is3Prime\n";
    

    DefineChainLayout($is3Prime);
    ($h1_s1, $h1_s2, $h2_s1, $h2_s2) = ( @{ $layout_array[$curr_region] } );
#    print "START:: STRAND INFO: $h1_s1, $h1_s2, $h2_s1, $h2_s2\n\n";


    $was_previous_cross = 1;
    $cross_counter = 0;
    $start = $end = 1;

    for $i (0 .. $#Regions) {
#	If there is a crossover at this point
#	print "REGION: " . ($i + 1) . " START: " . $Regions[$i]->{"Start"};
#	print " END: " . $Regions[$i]->{"End"} . "\n";
	if ( IsCrossoverInRegion(\@crossovers, \%{ $Regions[$i] }) && $was_previous_cross) {
	    if (($i + 1) % 2 == 0 && $cross_counter > 0) {
		($h1_s1, $h2_s2) = ($h2_s2, $h1_s1);
	    } elsif($cross_counter > 0) {
		($h1_s2, $h2_s1) = ($h2_s1, $h1_s2);
	    }
	    $cross_counter++;
	    $was_previous_cross = 1;
	}elsif ($was_previous_cross) {
	    if (($i + 1) % 2 == 0 && $cross_counter > 0) {
		($h1_s1, $h2_s2) = ($h2_s2, $h1_s1);
	    } elsif($cross_counter > 0) {
		($h1_s2, $h2_s1) = ($h2_s1, $h1_s2);
	    }
	    $was_previous_cross = 0;
	}elsif (IsCrossoverInRegion(\@crossovers, \%{ $Regions[$i] }) ) {
	    $was_previous_cross = 1;
	}

#	print "STRAND INFO: h1_s1: $h1_s1, h1_s2: $h1_s2, h2_s1: $h2_s1, ",
#	"h2_s2: $h2_s2\n\n";

#	Get The Starting and Ending Info
	$start = $Regions[$i]->{"Start"} - 1;
	$end = $Regions[$i]->{"End"};

#	Assign Strands
	$helix[0][0][$i]{"Strand"} = $h1_s1;
	$helix[0][1][$i]{"Strand"} = $h1_s2;

	$helix[1][1][$i]{"Strand"} = $h2_s1;
	$helix[1][0][$i]{"Strand"} = $h2_s2;

#	Write Start - End Data for Helix1
	#Strand 1
	$helix[0][0][$i]{"StartUnit"} = $start + 1;
	$helix[0][0][$i]{"EndUnit"} = $end;
	
	#Strand 2
	$helix[0][1][$i]{"EndUnit"} = $totalBases - $start;
	$helix[0][1][$i]{"StartUnit"} = $totalBases - $end + 1;
	    
#	Write Data for Helix2
	#Strand 1
	$helix[1][0][$i]{"StartUnit"} = $start + 1;
	$helix[1][0][$i]{"EndUnit"} = $end;

	#Strand 2
	$helix[1][1][$i]{"EndUnit"} = $totalBases - $start;
	$helix[1][1][$i]{"StartUnit"} = $totalBases - $end + 1;
	    
	$curr_region++;
	
    }

}

sub GetRegions {

    my ($majG, $minG, $tBases, $endBases) = @_;
    my ($period) = $majG + $minG;
    my (@RGions, $base_count, $base_counter);

    $base_count = $tBases - $endBases;

    $RGions[0] = (
		  {
		      "Start" => 1,
		      "End"   => $endBases,
		  }
		 );

    $base_counter = $endBases;

    while ($base_count >= $period) {
	$RGions[$#RGions +1] = (
				{
				    "Start" => $base_counter + 1,
				    "End"   => $base_counter + $majG,
				}
				);
	$base_counter += $majG;
	$RGions[$#RGions +1] = (
				{
				    "Start" => $base_counter + 1,
				    "End"   => $base_counter + $minG,
				}
				);

	$base_counter += $minG;
	$base_count -= $period;
    }
#    print "EXIT GetRegions. base_counter: $base_counter, base_count: $base_count\n";
    if ($base_count >= $minG) {
        $RGions[$#RGions +1] = (
                                {
                                   "Start" => $base_counter + 1,
                                   "End"   => $tBases,
                                }
                                );
    } 
    $RGions[$#RGions]->{"End"} = $tBases;

    return @RGions;
}
			       
sub IsCrossoverInRegion {
    my ($crossovers, $Curr_Region) = @_;

    my ($start_base, $end_base) = ($Curr_Region->{"Start"}, $Curr_Region->{"End"});
    my ($curr_cross, $returnval);

    for $curr_cross (@$crossovers) {
	if ($curr_cross >= $start_base && $curr_cross <= $end_base) {
	    $returnval = 1;
	    last;
	} else {
	    $returnval = 0;
	}
    }

    if ($returnval) {
#	print "GOT CROSSOVER BETWEEN $start_base AND $end_base\n";
    } else {
#	print "NO CROSSOVER BETWEEN $start_base AND $end_base\n";
    }

    return $returnval;
}

1;

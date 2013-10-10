package Packages::ManipImages;

BEGIN {
    push (@INC, "/ul/tpascal/.libs/PerlMagick/lib/perl5/site_perl/5.8.0/i386-linux-thread-multi/");
}

use Image::Magick;
use 5.005_64;
require Exporter;

our(@ISA, @EXPORT, $VERSION, @EXPORT_OK);

@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(AnnotatePic ConvertPicture ColasceImages DetermineRGBVal);
$VERSION = "1.00";

sub ConvertPicture {
    my ($in_img, $out_img) = @_;

    if (-e $in_img && $out_img) {
	$out_file = Image::Magick->new();
	$result = $out_file->ReadImage($in_img);
	if (! "$result") {
	    $out_file->Write( $out_img );
	    $result = 1;
	} else {
	    $result = 0;
	}
    } else {
	$result = 0;
    }

    return $result;
}

sub AnnotatePic {
    my ($result, $img_height, $img_weight);
    my ($in_img, $in_text, $text_size, $l_offset, $w_offset, $l_loc, $w_loc) = @_;

    if (-e $in_img) {
	$out_file = Image::Magick->new();
	$result = $out_file->ReadImage($in_img);
	if (! $result) {
	    ($img_height, $img_width) = $out_file->Get('height', 'width');
	    if ($l_loc == 2) {
		$l_offset += ($img_height/2);
	    } elsif($l_loc == 3) {
		$l_offset = $img_height - $l_offset;
	    }
	    
	    if ($w_loc == 2) {
		$w_offset = ($img_width/2) + $w_offset;
	    } elsif ($w_loc == 3) {
		$w_offset = $img_width - $w_offset;
	    }
	    
	    $out_file->Annotate(
				text=>$in_text, 
				pointsize=>$text_size, 
				fill=>'white', 
				align=>'Center', 
				x=>$w_offset, 
				y=>$l_offset
				);
	    $out_file->Write($in_img);
	    $result = 1;
	} else {
	    $result = 0;
	}
    } else {
	$result = 0;
    }
    return $result;
}

sub ColasceImages {
# This will colasce two images with the following options:
#   $l_loc: VALUES 1 - 3, represents the positions right, middle, left
#   $w_loc: VALUES 1 - 3, represents the positions top, middile, bottom
#   $l_offset: VALUES int, represents the offset for the length
#   $w_offset: VALUES int, represents the offset for the width

    my ($img1, $img2, $save_nm, $l_offset, $w_offset, $l_loc, $w_loc) = @_;
    my ($MainImage, $ChildImage, $result);
    my ($img_height, $img_width, $height_place, $width_place);

    if (-e $img1 && -e $img2) {
	$MainImage = Image::Magick->new();
	$result = $MainImage->ReadImage($img1);
	if (! $result) {
	    ($img_height, $img_width) = $MainImage->Get('height', 'width');
	    $ChildImage = Image::Magick->new();
	    $result = $ChildImage->ReadImage($img2);
	    if (! $result) {
		if ($l_loc == 2) {
		    $l_offset += ($img_height/2);
		} elsif($l_loc == 3) {
		    $l_offset = $img_height - $l_offset;
		}
	    
		if ($w_loc == 2) {
		    $w_offset = ($img_width/2) + $w_offset;
		} elsif ($w_loc == 3) {
		    $w_offset = $img_width - $w_offset;
		}
	    
		$MainImage->Composite(
				      image=>$ChildImage, 
				      compose=>'over', 
				      x=>$w_offset, 
				      y=>$l_offset,
				      );
		$MainImage->Write($save_nm);
		$result = 1;
	    } else {
		$result = 0;
	    }
	} else {
	    $result = 0;
	}
    }
}

sub DetermineRGBVal {
    my ($curr_val, $max_val, $min_val) = @_;
    my ($return_val);
    
    # there are 255x255x255 different colors
    # $min_val = 000 000 000
    # $max_val = 255 000 000

    $return_val = ($max_val - $min_val);
    $return_val = ($curr_val - $min_val) / $return_val;
#    $return_val = sprintf( "%03d", int($return_val * 255) );

#    $return_val .= " 000 000";

    ($return_val < 1/3) ?
	return "white" :
	($return_val < 2/3) ?
	    return "green" :
	    return "red";
	
}

1;

package Math::FFT;

use strict;
use vars qw($VERSION @ISA);
require DynaLoader;

@ISA = qw(DynaLoader);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

$VERSION = '1.28';

bootstrap Math::FFT $VERSION;

# Preloaded methods go here.

sub new {
  my ($class, $data) = @_;
  die 'Must call constructor with an array reference for the data'
    unless ref($data) eq 'ARRAY';
  $data->[0] ||= 0;    # keep warnings happy
  my $n = @$data;
  my $nip = int(3 + sqrt($n));
  my $nw = int(2 + 5*$n/4);
  my $ip = pack("i$nip", ());
  my $w = pack("d$nw", ());
  bless {
	 type => '',
	 mean => '',
	 coeff => '',
	 n   => $n,
	 data => $data,
         ip => \$ip,
         w => \$w,
	}, $class;
}

# clone method to copy the ip and w arrays for data of equal size
sub clone {
  my ($self, $data) = @_;
  die 'Must call clone with an array reference for the data'
    unless ref($data) eq 'ARRAY';
  $data->[0] ||= 0;    # keep warnings happy
  my $n = @$data;
  die "Cannot clone data of unequal sizes" unless $n == $self->{n};
  my $class = ref($self);
  bless {
	 type => '',
	 coeff => '',
	 mean => '',
	 n   => $self->{n},
	 data => $data,
	 ip => $self->{ip},
	 w => $self->{w},
	}, $class;
}

# Complex Discrete Fourier Transform
sub cdft {
  my $self = shift;
  my $n = $self->{n};
  die "data size ($n) must be an integer power of 2" unless check_n($n);
  my $data = [ @{$self->{data}} ];
  _cdft($n, 1, $data, $self->{ip}, $self->{w});
  $self->{type} = 'cdft';
  $self->{coeff} = $data;
  return $data;
}

# Inverse Complex Discrete Fourier Transform
sub invcdft {
  my $self = shift;
  my $data;
  my $n = $self->{n};
  if (my $arg = shift) {
    die 'Must pass an array reference to invcdft'
      unless ref($arg) eq 'ARRAY';
    die "Size of data set must be $n"
      unless $n == @$arg;
    $data = [ @$arg ];
  }
  else {
    die 'Must invert data created with cdft'
      unless $self->{type} eq 'cdft';
    $data =  [ @{$self->{coeff}} ];
  }
  _cdft($n, -1, $data, $self->{ip}, $self->{w});
  $_ *= 2.0/$n for (@$data);
  return $data;
}

# Real Discrete Fourier Transform
sub rdft {
  my $self = shift;
  my $n = $self->{n};
  die "data size ($n) must be an integer power of 2" unless check_n($n);
  my $data = [ @{$self->{data}} ];
  _rdft($n, 1, $data, $self->{ip}, $self->{w});
  $self->{type} = 'rdft';
  $self->{coeff} = $data;
  return $data;
}

# Inverse Real Discrete Fourier Transform
sub invrdft {
  my $self = shift;
  my $data;
  my $n = $self->{n};
  if (my $arg = shift) {
    die 'Must pass an array reference to invrdft'
      unless ref($arg) eq 'ARRAY';
    die "Size of data set must be $n"
      unless $n == @$arg;
    $data = [ @$arg ];
  }
  else {
    die 'Must invert data created with rdft'
      unless $self->{type} eq 'rdft';
    $data =  [ @{$self->{coeff}} ];
  }
  _rdft($n, -1, $data, $self->{ip}, $self->{w});
  $_ *= 2.0/$n for (@$data);
  return $data;
}

# Discrete Cosine Transform
sub ddct {
  my $self = shift;
  my $n = $self->{n};
  die "data size ($n) must be an integer power of 2" unless check_n($n);
  my $data = [ @{$self->{data}} ];
  _ddct($n, -1, $data, $self->{ip}, $self->{w});
  $self->{type} = 'ddct';
  $self->{coeff} = $data;
  return $data;
}

# Inverse Discrete Cosine Transform
sub invddct {
  my $self = shift;
  my $data;
  my $n = $self->{n};
  if (my $arg = shift) {
    die 'Must pass an array reference to invddct'
      unless ref($arg) eq 'ARRAY';
    die "Size of data set must be $n"
      unless $n == @$arg;
    $data = [ @$arg ];
  }
  else {
    die 'Must invert data created with ddct'
      unless $self->{type} eq 'ddct';
    $data =  [ @{$self->{coeff}} ];
  }
  $data->[0] *= 0.5;
  _ddct($n, 1, $data, $self->{ip}, $self->{w});
  $_ *= 2.0/$n for (@$data);
  return $data;
}

# Discrete Sine Transform
sub ddst {
  my $self = shift;
  my $n = $self->{n};
  die "data size ($n) must be an integer power of 2" unless check_n($n);
  my $data = [ @{$self->{data}} ];
  _ddst($n, -1, $data, $self->{ip}, $self->{w});
  $self->{type} = 'ddst';
  $self->{coeff} = $data;
  return $data;
}

# Inverse Discrete Sine Transform
sub invddst {
  my $self = shift;
  my $data;
  my $n = $self->{n};
  if (my $arg = shift) {
    die 'Must pass an array reference to invddst'
      unless ref($arg) eq 'ARRAY';
    die "Size of data set must be $n"
      unless $n == @$arg;
    $data = [ @$arg ];
  }
  else {
    die 'Must invert data created with ddst'
      unless $self->{type} eq 'ddst';
    $data =  [ @{$self->{coeff}} ];
  }
  $data->[0] *= 0.5;
  _ddst($n, 1, $data, $self->{ip}, $self->{w});
  $_ *= 2.0/$n for (@$data);
  return $data;
}

# Cosine Transform of RDFT (Real Symmetric DFT)
sub dfct {
  my $self = shift;
  my $np1 = $self->{n};
  my $n = $np1 - 1;
  die "data size ($n) must be an integer power of 2" unless check_n($n);
  my $nt = int(2 + $n/2);
  my $t = [];
  my $data = [ @{$self->{data}} ];
  pdfct($nt, $n, $data, $t, $self->{ip}, $self->{w});
  $self->{type} = 'dfct';
  $self->{coeff} = $data;
  return $data;
}

# Inverse Cosine Transform of RDFT (Real Symmetric DFT)
sub invdfct {
  my $self = shift;
  my $data;
  my $np1 = $self->{n};
  my $n = $np1 - 1;
  if (my $arg = shift) {
    die 'Must pass an array reference to invdfct'
      unless ref($arg) eq 'ARRAY';
    die "Size of data set must be $n"
      unless $np1 == @$data;
    $data = [ @$arg ];
  }
  else {
    die 'Must invert data created with dfct'
      unless $self->{type} eq 'dfct';
    $data =  [ @{$self->{coeff}} ];
  }
  my $nt = int(2 + $n/2);
  my $t = [];
  $data->[0] *= 0.5;
  $data->[$n] *= 0.5;
  pdfct($nt, $n, $data, $t, $self->{ip}, $self->{w});
  $data->[0] *= 0.5;
  $data->[$n] *= 0.5;
  $_ *= 2.0/$n for (@$data);
  return $data;
}

# Sine Transform of RDFT (Real Anti-symmetric DFT)
sub dfst {
  my $self = shift;
  my $n = $self->{n};
  die "data size ($n) must be an integer power of 2" unless check_n($n);
  my $data = [ @{$self->{data}} ];
  my $nt = int(2 + $n/2);
  my $t = [];
  pdfst($nt, $n, $data, $t, $self->{ip}, $self->{w});
  $self->{type} = 'dfst';
  $self->{coeff} = $data;
  return $data;
}

# Inverse Sine Transform of RDFT (Real Anti-symmetric DFT)
sub invdfst {
  my $self = shift;
  my $n = $self->{n};
  my $data;
  if (my $arg = shift) {
    die 'Must pass an array reference to invdfst'
      unless ref($arg) eq 'ARRAY';
    die "Size of data set must be $n"
      unless $n == @$arg;
    $data = [ @$arg ];
  }
  else {
    die 'Must invert data created with dfst'
      unless $self->{type} eq 'dfst';
    $data =  [ @{$self->{coeff}} ];
  }
  my $nt = int(2 + $n/2);
  my $t = [];
  pdfst($nt, $n, $data, $t, $self->{ip}, $self->{w});
  $_ *= 2.0/$n for (@$data);
  return $data;
}

# check if $n is a power of 2
sub check_n {
  my $n = shift;
  my $y = log($n) / 0.693147180559945309417;
  return abs($y-int($y)) < 1e-6 ? 1 : 0;
}

sub correl {
  my ($self, $other) = @_;
  my $n = $self->{n};
  my $d1 = $self->{type} ? 
    ($self->{type} eq 'rdft' ? [ @{$self->{coeff}} ] :
    die 'correl must involve a real function' ) :
      $self->rdft &&  [ @{$self->{coeff}} ];
  my $d2 = [];
  if (ref($other) eq 'Math::FFT') {
    $d2 = $other->{type} ? 
      ($other->{type} eq 'rdft' ? [ @{$other->{coeff}}] :
       die 'correl must involve a real function' ) :
	 $other->rdft && [ @{$other->{coeff}}];
  }
  elsif (ref($other) eq 'ARRAY') {
    $d2 = [ @$other ];
    _rdft($n, 1, $d2, $self->{ip}, $self->{w});
  }
  else {
    die 'Must call correl with either a Math::FFT object or an array ref';
  }
  my $corr = [];
  _correl($n, $corr, $d1, $d2, $self->{ip}, $self->{w});
  return $corr;
}

sub convlv {
  my ($self, $r) = @_;
  die 'Must call convlv with an array reference for the response data'
    unless ref($r) eq 'ARRAY';
  my $respn = [ @$r ];
  my $m = @$respn;
  die 'size of response data must be an odd integer' unless $m % 2 == 1;
  my $n = $self->{n};
  my $d1 = $self->{type} ? 
    ($self->{type} eq 'rdft' ? [ @{$self->{coeff}} ] :
    die 'correl must involve a real function' ) :
      $self->rdft &&  [ @{$self->{coeff}} ];
  for (my $i=1; $i<=($m-1)/2; $i++) {
    $respn->[$n-$i] = $respn->[$m-$i];
  }
  for (my $i=($m+3)/2; $i<=$n-($m-1)/2; $i++) {
    $respn->[$i-1] = 0.0;
  }
  my $convlv = [];
  _convlv($n, $convlv, $d1, $respn, $self->{ip}, $self->{w});
  return $convlv;
}

sub deconvlv {
  my ($self, $r) = @_;
  die 'Must call deconvlv with an array reference for the response data'
    unless ref($r) eq 'ARRAY';
  my $respn = [ @$r ];
  my $m = @$respn;
  die 'size of response data must be an odd integer' unless $m % 2 == 1;
  my $n = $self->{n};
  my $d1 = $self->{type} ? 
    ($self->{type} eq 'rdft' ? [ @{$self->{coeff}} ] :
    die 'correl must involve a real function' ) :
      $self->rdft &&  [ @{$self->{coeff}} ];
  for (my $i=1; $i<=($m-1)/2; $i++) {
    $respn->[$n-$i] = $respn->[$m-$i];
  }
  for (my $i=($m+3)/2; $i<=$n-($m-1)/2; $i++) {
    $respn->[$i-1] = 0.0;
  }
  my $convlv = [];
  if (_deconvlv($n, $convlv, $d1, $respn, $self->{ip}, $self->{w}) != 0) {
    die "Singularity encountered for response in deconvlv";
  }
  return $convlv;
}

sub spctrm {
  my ($self, %args) = @_;
  my %accept = map {$_ => 1} qw(window segments number overlap);
  for (keys %args) {
    die "`$_' is not a valid argument to spctrm" if not $accept{$_};
  }
  my $win_fun = $args{window};
  if ($win_fun and ref($win_fun) ne 'CODE') {
    my %accept = map {$_ => 1} qw(hamm hann welch bartlett);
    die "`$win_fun' is not a known window function in spctrm"
      if not $accept{$win_fun};
  }
  die 'Please specify a value for "segments" in spctrm()'
    if ($args{number} and ! $args{segments}); 
  my $n = $self->{n};
  my $d;
  my $n2 = 0;
  my $spctrm = [];
  my $win_sub = {
		 'hamm' => sub {
		   my ($j, $n) = @_;
		   my $pi = 4.0*atan2(1,1);
		   return (1 - cos(2*$pi*$j/$n))/2;
		 },
		 'hann' => sub {
		   my ($j, $n) = @_;
		   my $pi = 4.0*atan2(1,1);
		   return (1 - cos(2*$pi*$j/$n))/2;
		 },
		 'welch' => sub {
		   my ($j, $n) = @_;
		   return 1 - 4*($j-$n/2)*($j-$n/2)/$n/$n;
		 },
		 'bartlett' => sub {
		   my ($j, $n) = @_;
		   return 1 - abs(2*($j-$n/2)/$n);
		 },
		};
  if (not $args{segments} or ($args{segments} == 1 and not $args{number})) {
    die "data size ($n) must an integer power of 2" unless check_n($n);
    if ($win_fun) {
      $d = [ @{$self->{data}}];
      $win_fun = $win_sub->{$win_fun} if ref($win_fun) ne 'CODE';
      for (my $j=0; $j<$n; $j++) {
	my $w = $win_fun->($j, $n);
	$d->[$j] *= $w;
	$n2 += $w * $w;
      }
      $n2 *= $n;
      _spctrm($n, $spctrm, $d, $self->{ip}, $self->{w}, $n2, 1);
    }
    else {
      $d = $self->{type} ? 
	($self->{type} eq 'rdft' ? $self->{coeff} :
	 die 'correl must involve a real function' ) :
	   $self->rdft && $self->{coeff};
      $n2 = $n*$n;
      _spctrm($n, $spctrm, $d, $self->{ip}, $self->{w}, $n2, 0);
    }
  }
  else {
    $d = [ @{$self->{data}}];
    my ($data, @w);
    my $k = $args{segments};
    my $m = $args{number};
    die 'Please specify a value for "number" in spctrm()'
      if ($k and ! $m); 
    die "number ($m) must an integer power of 2" unless check_n($m);
    my $m2 = $m+$m;
    my $overlap = $args{overlap};
    my $N = $overlap ? ($k+1)*$m : 2*$k*$m;
    die "Need $N data points (data only has $n)" if $N > $n;
    if ($win_fun) {
      $win_fun = $win_sub->{$win_fun} if ref($win_fun) ne 'CODE';
      for (my $j=0; $j<$m2; $j++) {
	$w[$j] = $win_fun->($j, $m2);
	$n2 += $w[$j]*$w[$j];
      }
    }
    else {
      $n2 = $m2;      
    }
    if ($overlap) {
      my @old =  splice(@$d, 0, $m);
      for (0..$k-1) {
	push @{$data->[$_]}, @old;
	my @new = splice(@$d, 0, $m);
	push @{$data->[$_]}, @new;
	@old = @new;
	if ($win_fun) {
	  my $j=0;
	  $data->[$_] = [ map {$w[$j++]*$_} @{$data->[$_]}];
	}
      }
    }
    else {
      for (0..$k-1) {
	push @{$data->[$_]}, splice(@$d, 0, $m2);
	if ($win_fun) {
	  my $j=0;
	  $data->[$_] = [ map {$w[$j++]*$_} @{$data->[$_]}];
	}
      }
    }
    my $tmp = [];
    my $nip = int(3 + sqrt($m2));
    my $nw = int(2 + 5*$m2/4);
    my $ip = pack("i$nip", ());
    my $w = pack("d$nw", ());
    _spctrm_bin($k, $m2, $spctrm, $data, \$ip, \$w, $n2, $tmp);
  }
  return $spctrm;
}

sub mean {
  my $self = shift;
  my $sum = 0;
  my ($n, $data);
  my $flag = 0;
  if ($data = shift) {
    die 'Must call with an array reference'
      unless ref($data) eq 'ARRAY';
    $n = @$data;
    $flag = 1;
  }
  else {
    $data = $self->{data};
    $n = $self->{n};
  }
  $sum += $_ for @$data;
  my $mean = $sum / $n;
  $self->{mean} = $mean unless $flag == 1;
  return $mean;
}

sub rms {
  my $self = shift;
  my $sum = 0;
  my ($n, $data);
  if ($data = shift) {
    die 'Must call with an array reference'
      unless ref($data) eq 'ARRAY';
    $n = @$data;
  }
  else {
    $data = $self->{data};
    $n = $self->{n};
  }
  $sum += $_*$_ for @$data;
  return sqrt($sum / $n);
}

sub stdev {
  my $self = shift;
  my ($n, $data, $mean);
  if ($data = shift) {
    die 'Must call with an array reference'
      unless ref($data) eq 'ARRAY';
    $n = @$data;
    $mean = $self->mean($data);
  }
  else {
    $data = $self->{data};
    $n = $self->{n};
    $mean = $self->{mean} || $self->mean;
  }
  die 'Cannot find the standard deviation with n = 1'
    if $n == 1;
  my $sum = 0;
  $sum += ($_ - $mean)*($_ - $mean) for @$data;
  return sqrt($sum / ($n-1));
}

sub range {
  my $self = shift;
  my ($n, $data);
  if ($data = shift) {
    die 'Must call with an array reference'
      unless ref($data) eq 'ARRAY';
    $n = @$data;
  }
  else {
    $data = $self->{data};
    $n = $self->{n};
  }
  my $min = $data->[0];
  my $max = $data->[0];
  for (@$data) {
    $min = $_ if $_ < $min;
    $max = $_ if $_ > $max;
  }
  return ($min, $max);
}

sub median {
  my $self = shift;
  my ($n, $data);
  if ($data = shift) {
    die 'Must call with an array reference'
      unless ref($data) eq 'ARRAY';
    $n = @$data;
  }
  else {
    $data = $self->{data};
    $n = $self->{n};
  }
  my @sorted = sort {$a <=> $b} @$data;
  return $n % 2 == 1 ?
    $sorted[($n-1)/2] : ($sorted[$n/2] + $sorted[$n/2-1])/2;
}

# Autoload methods go after =cut, and are processed by the autosplit program.


1;

__END__

=head1 NAME

Math::FFT - Perl module to calculate Fast Fourier Transforms

=head1 SYNOPSIS

  use Math::FFT;
  my $PI = 3.1415926539;
  my $N = 64;
  my ($series, $other_series);
  for (my $k=0; $k<$N; $k++) {
      $series->[$k] = sin(4*$k*$PI/$N) + cos(6*$k*$PI/$N);
  }
  my $fft = new Math::FFT($series);
  my $coeff = $fft->rdft();
  my $spectrum = $fft->spctrm;
  my $original_data = $fft->invrdft($coeff);

  for (my $k=0; $k<$N; $k++) {
      $other_series->[$k] = sin(16*$k*$PI/$N) + cos(8*$k*$PI/$N);
  }
  my $other_fft = $fft->clone($other_series);
  my $other_coeff = $other_fft->rdft();
  my $correlation = $fft->correl($other_fft);

=head1 DESCRIPTION

This module implements some algorithms for calculating
Fast Fourier Transforms for one-dimensional data sets of size 2^n.
The data, assumed to arise from a constant sampling rate, is
represented by an array reference C<$data> (as described in the
methods below), which is then used to create a C<Math::FFT> object as

  my $fft = new Math::FFT($data);

The methods available include the following.

=head2 FFT METHODS

=over

=item C<$coeff = $fft-E<gt>cdft();>

This calculates the complex discrete Fourier transform
for a data set C<x[j]>. Here, C<$data> is a reference to an 
array C<data[0...2*n-1]> holding the data

  data[2*j] = Re(x[j]),
  data[2*j+1] = Im(x[j]), 0<=j<n

An array reference C<$coeff> is returned consisting of

  coeff[2*k] = Re(X[k]),
  coeff[2*k+1] = Im(X[k]), 0<=k<n

where

   X[k] = sum_j=0^n-1 x[j]*exp(2*pi*i*j*k/n), 0<=k<n

=item C<$orig_data = $fft-E<gt>invcdft([$coeff]);>

Calculates the inverse complex discrete Fourier transform
on a data set C<x[j]>. If C<$coeff> is not given, it will be set 
equal to an earlier call to C<$fft-E<gt>cdft()>. C<$coeff> is 
a reference to an array C<coeff[0...2*n-1]> holding the data

  coeff[2*j] = Re(x[j]),
  coeff[2*j+1] = Im(x[j]), 0<=j<n

An array reference C<$orig_data> is returned consisting of

  orig_data[2*k] = Re(X[k]),
  orig_data[2*k+1] = Im(X[k]), 0<=k<n

where, excluding the scale,

   X[k] = sum_j=0^n-1 x[j]*exp(-2*pi*i*j*k/n), 0<=k<n

A scaling C<$orig_data-E<gt>[$i] *= 2.0/$n> is then done so that
C<$orig_data> coincides with the original C<$data>.

=item C<$coeff = $fft-E<gt>rdft();>

This calculates the real discrete Fourier transform
for a data set C<x[j]>. On input, $data is a reference to an
array C<data[0...n-1]> holding the data. An array reference
C<$coeff> is returned consisting of

  coeff[2*k] = R[k], 0<=k<n/2
  coeff[2*k+1] = I[k], 0<k<n/2
  coeff[1] = R[n/2]

where

  R[k] = sum_j=0^n-1 data[j]*cos(2*pi*j*k/n), 0<=k<=n/2
  I[k] = sum_j=0^n-1 data[j]*sin(2*pi*j*k/n), 0<k<n/2

=item C<$orig_data = $fft-E<gt>invrdft([$coeff]);>

Calculates the inverse real discrete Fourier transform
on a data set C<coeff[j]>. If C<$coeff> is not given, it will be set 
equal to an earlier call to C<$fft-E<gt>rdft()>. C<$coeff> 
is a reference to an array C<coeff[0...n-1]> holding the data

  coeff[2*j] = R[j], 0<=j<n/2
  coeff[2*j+1] = I[j], 0<j<n/2
  coeff[1] = R[n/2]

An array reference C<$orig_data> is returned where, excluding the scale,

  orig_data[k] = (R[0] + R[n/2]*cos(pi*k))/2 + 
    sum_j=1^n/2-1 R[j]*cos(2*pi*j*k/n) + 
      sum_j=1^n/2-1 I[j]*sin(2*pi*j*k/n), 0<=k<n

A scaling C<$orig_data-E<gt>[$i] *= 2.0/$n> is then done so that
C<$orig_data> coincides with the original C<$data>.

=item C<$coeff = $fft-E<gt>ddct();>

Computes the discrete cosine tranform on a data set
C<data[0...n-1]> contained in an array reference C<$data>. An
array reference C<$coeff> is returned consisting of

  coeff[k] = C[k], 0<=k<n

where

  C[k] = sum_j=0^n-1 data[j]*cos(pi*(j+1/2)*k/n), 0<=k<n 

=item C<$orig_data = $fft-E<gt>invddct([$coeff]);>

Computes the inverse discrete cosine tranform on a data set
C<coeff[0...n-1]> contained in an array reference C<$coeff>. 
If C<$coeff> is not given, it will be set equal to an earlier 
call to C<$fft-E<gt>ddct()>. An array reference C<$orig_data> 
is returned consisting of

  orig_data[k] = C[k], 0<=k<n

where, excluding the scale,

  C[k] = sum_j=0^n-1 coeff[j]*cos(pi*j*(k+1/2)/n), 0<=k<n

A scaling C<$orig_data-E<gt>[$i] *= 2.0/$n> is then done so that
C<$orig_data> coincides with the original C<$data>.

=item C<$coeff = $fft-E<gt>ddst();>

Computes the discrete sine transform of a data set 
C<data[0...n-1]> contained in an array reference C<$data>. An
array reference C<$coeff> is returned consisting of

 coeff[k] = S[k], 0<k<n
 coeff[0] = S[n]

where

 S[k] = sum_j=0^n-1 data[j]*sin(pi*(j+1/2)*k/n), 0<k<=n

=item C<$orig_data = $fft-E<gt>invddst($coeff);>

Computes the inverse discrete sine transform of a data set 
C<coeff[0...n-1]> contained in an array reference C<$coeff>, arranged as 

 coeff[j] = A[j], 0<j<n
 coeff[0] = A[n]

If C<$coeff> is not given, it will be set equal to an earlier 
call to C<$fft-E<gt>ddst()>. An array reference C<$orig_data> 
is returned consisting of

 orig_data[k] = S[k], 0<=k<n

where, excluding a scale,

 S[k] =  sum_j=1^n A[j]*sin(pi*j*(k+1/2)/n), 0<=k<n

The scaling C<$a-E<gt>[$i] *= 2.0/$n> is then done so that
C<$orig_data> coincides with the original C<$data>.

=item C<$coeff = $fft-E<gt>dfct();>

Computes the real symmetric discrete Fourier transform of a
data set C<data[0...n]> contained in the array reference C<$data>. An
array reference C<$coeff> is returned consisting of 

  coeff[k] = C[k], 0<=k<=n

where

  C[k] = sum_j=0^n data[j]*cos(pi*j*k/n), 0<=k<=n

=item C<$orig_data = $fft-E<gt>invdfct($coeff);>

Computes the inverse real symmetric discrete Fourier transform of a
data set C<coeff[0...n]> contained in the array reference C<$coeff>. 
If C<$coeff> is not given, it will be set equal to an earlier 
call to C<$fft-E<gt>dfct()>. An array reference C<$orig_data> 
is returned consisting of

  orig_data[k] = C[k], 0<=k<=n

where, excluding the scale,

  C[k] = sum_j=0^n coeff[j]*cos(pi*j*k/n), 0<=k<=n

A scaling C<$coeff-E<gt>[0] *= 0.5>, C<$coeff-E<gt>[$n] *= 0.5>, and 
C<$orig_data-E<gt>[$i] *= 2.0/$n> is then done so that
C<$orig_data> coincides with the original C<$data>.

=item C<$coeff = $fft-E<gt>dfst();>

Computes the real anti-symmetric discrete Fourier transform of a
data set C<data[0...n-1]> contained in the array reference C<$data>. An
array reference C<$coeff> is returned consisting of 

  coeff[k] = C[k], 0<k<n

where

  C[k] = sum_j=0^n data[j]*sin(pi*j*k/n), 0<k<n

(C<coeff[0]> is used for a work area)

=item C<$orig_data = $fft-E<gt>invdfst($coeff);>

Computes the inverse real anti-symmetric discrete Fourier transform of a
data set C<coeff[0...n-1]> contained in the array reference C<$coeff>.
If C<$coeff> is not given, it will be set equal to an earlier 
call to C<$fft-E<gt>dfst()>. An array reference C<$orig_data> is 
returned consisting of

  orig_data[k] = C[k], 0<k<n

where, excluding the scale,

  C[k] = sum_j=0^n coeff[j]*sin(pi*j*k/n), 0<k<n

A scaling C<$orig_data-E<gt>[$i] *= 2.0/$n> is then done so that
C<$orig_data> coincides with the original C<$data>.

=back

=head2 CLONING

The algorithm used in the transforms makes use of arrays for a work 
area and for a cos/sin lookup table dependent only on the size of 
the data set. These arrays are initialized when the C<Math::FFT> object 
is created and then are populated when a transform method is first 
invoked. After this, they persist for the lifetime of the object.

This aspect is exploited in a C<cloning> method; if a C<Math::FFT>
object is created for a data set C<$data1> of size C<N>:

  $fft1 = new Math::FFT($data1);

then a new C<Math::FFT> object can be created for a second data 
set C<$data2> of the I<same> size C<N> by

   $fft2 = $fft1->clone($data2);

The C<$fft2> object will copy the reuseable work area and
lookup table calculated from C<$fft1>.

=head2 APPLICATIONS

This module includes some common applications - correlation,
convolution and deconvolution, and power spectrum - that
arise with real data sets. The conventions used here
follow that of I<Numerical Recipes in C>, by Press, Teukolsky,
Vetterling, and Flannery, in which further details of the
algorithms are given. Note in particular the treatment of end
effects by zero padding, which is assumed to be done by the
user, if required.

=over

=item Correlation

The correlation between two functions is defined as

             /
   Corr(t) = | ds g(s+t) h(s) 
             /

This may be calculated, for two array references C<$data1>
and C<$data2> of the same size C<$n>, as either

   $fft1 = new Math::FFT($data1);
   $fft2 = new Math::FFT($data2);
   $corr = $fft1->correl($fft2);

or as

   $fft1 = new Math::FFT($data1);
   $corr = $fft1->correl($data2);

The array reference C<$corr> is returned in wrap-around 
order - correlations at increasingly positive lags are in 
C<$corr-E<gt>[0]> (zero lag) on up to C<$corr-E<gt>[$n/2-1]>, 
while correlations at increasingly negative lags are in 
C<$corr-E<gt>[$n-1]> on down to C<$corr-E<gt>[$n/2]>. The sign 
convention used is such that if C<$data1> lags C<$data2> (that 
is, is shifted to the right), then C<$corr> will show a peak 
at positive lags.

=item Convolution

The convolution of two functions is defined as

               /
   Convlv(t) = | ds g(s) h(t-s) 
               /

This is similar to calculating the correlation between the
two functions, but typically the functions here have a quite
different physical interpretation - one is a signal which 
persists indefinitely in time, and the other is a response 
function of limited duration. The convolution may be calculated, 
for two array references C<$data> and C<$respn>, as

   $fft = new Math::FFT($data);
   $convlv = $fft->convlv($respn);

with the returned C<$convlv> being an array reference. The method 
assumes that the response function C<$respn> has an I<odd> number 
of elements C<$m> less than or equal to the number of elements C<$n> 
of C<$data>. C<$respn> is assumed to be stored in wrap-around order - 
the first half contains the response at positive times, while the 
second half, counting down from C<$respn-E<gt>[$m-1]>, contains the
response at negative times.

=item Deconvolution

Deconvolution undoes the effects of convoluting a signal
with a known response function. In other words, in the relation

               /
   Convlv(t) = | ds g(s) h(t-s) 
               /

deconvolution reconstructs the original signal, given the convolution
and the response function. The method is implemented, for two array 
references C<$data> and C<$respn>, as

   $fft = new Math::FFT($data);
   $deconvlv = $fft->deconvlv($respn);

As a result, if the convolution of a data set C<$data> with
a response function C<$respn> is calculated as

   $fft1 = new Math::FFT($data);
   $convlv = $fft1->convlv($respn);

then the deconvolution

   $fft2 = new Math::FFT($convlv);
   $deconvlv = $fft2->deconvlv($respn);

will give an array reference C<$deconvlv> containing the
same elements as the original data C<$data>.

=item Power Spectrum

If the FFT of a real function of C<N> elements is calculated, 
the C<N/2+1> elements of the power spectrum are defined, in terms 
of the (complex) Fourier coefficients C<C[k]>, as

   P[0] = |C[0]|^2 / N^2
   P[k] = 2 |C[k]|^2 / N^2   (k = 1, 2 ,..., N/2-1)
   P[N/2] = |C[N/2]|^2 / N^2

Often for these purposes the data is partitioned into C<K>
segments, each containing C<2M> elements. The power spectrum
for each segment is calculated, and the net power spectrum
is the average of all of these segmented spectra.

Partitioning may be done in one of two ways: I<non-overlapping> and
I<overlapping>. Non-overlapping is useful when the data set
is gathered in real time, where the number of data points
can be varied at will. Overlapping is useful where there is
a fixed number of data points. In non-overlapping, the first
<2M> elements constitute segment 1, the next C<2M> elements
are segment 2, and so on up to segment C<K>, for a total of
C<2KM> sampled points. In overlapping, the first and second
C<M> elements are segment 1, the second and third C<M> elements
are segment 2, and so on, for a total of C<(K+1)M> sampled points.

A problem that may arise in this procedure is I<leakage>: the
power spectrum calculated for one bin contains contributions
from nearby bins. To lessen this effect I<data windowing> is
often used: multiply the original data C<d[j]> by a window
function C<w[j]>, where j = 0, 1, ..., N-1. Some popular choices 
of such functions are

              | j - N/2 |
  w[j] = 1 -  | ------- |     ... Bartlett   
              |   N/2   |


              / j - N/2 \ 2
  w[j] = 1 -  | ------- |     ... Welch  
              \   N/2   /


           1   /                    \
  w[j] =  ---  |1 - cos(2 pi j / N) |     ... Hann  
           2   \                    /


The C<spctrm> method, used as

    $fft = Math::FFT->new($data);
    $spectrum = $fft->spctrm(%options);

returns an array reference C<$spectrum> representing the power 
spectrum for a data set represented by an array reference C<$data>.
The options available are

=over

=item C<window =E<gt> window_name>

This specifies the window function; if not given, no such
function is used. Accepted values (see above) are C<"bartlett">, 
C<"welch">, C<"hann">, and C<\&my_window>, where C<my_window> is a 
user specified subroutine which must be of the form, for example,

   sub my_window {
      my ($j, $n) = @_;
      return 1 - abs(2*($j-$n/2)/$n);
   }

which implements the Bartlett window.

=item C<overlap =E<gt> 1>

This specifies whether overlapping should be done; if true (1),
overlapping will be used, whereas if false (0), or not
specified, no overlapping is used.

=item C<segments =E<gt> n>

This specifies that the data will be partitioned into C<n>
segments. If not specified, no segmentation will be done.

=item C<number =E<gt> m>

This specifies that C<2m> data points will be used for 
each segment, and must be a power of 2. The power 
spectrum returned will consist of C<m+1> elements.

=back

=back

=head2 STATISTICAL FUNCTIONS

For convenience, a number of common statistical functions are 
included for analyzing real data. After creating the object as

  my $fft = new Math::FFT($data);

for a data set represented by the array reference C<$data>
of size C<N>, these methods may be called as follows.

=over

=item C<$mean = $fft-E<gt>mean([$data]);>

This returns the mean

  1/N * sum_j=0^N-1 data[j]

If an array reference C<$data> is not given, the data set used 
in creating C<$fft> will be used.

=item C<$stdev = $fft-E<gt>stdev([$data]);>

This returns the standard deviation

  sqrt{ 1/(N-1) * sum_j=0^N-1 (data[j] - mean)**2 }

If an array reference C<$data> is not given, the data set used 
in creating C<$fft> will be used.

=item C<$rms = $fft-E<gt>rms([$data]);>

This returns the root mean square

  sqrt{ 1/N * sum_j=0^N-1 (data[j])**2 }

If an array reference C<$data> is not given, the data set used 
in creating C<$fft> will be used.

=item C<($min, $max) = $fft-E<gt>range([$data]);>

This returns the minimum and maximum values of the data set.
If an array reference C<$data> is not given, the data set used 
in creating C<$fft> will be used.

=item C<$median = $fft-E<gt>median([$data]);>

This returns the median of a data set. The median is defined,
for the I<sorted> data set, as either the middle element, if the
number of elements is odd, or as the interpolated value of
the the two values on either side of the middle, if the number
of elements is even. If an array reference C<$data> is not given, 
the data set used in creating C<$fft> will be used.

=back

=head1 BUGS

Please report any to Randy Kobes <randy@theoryx5.uwinnipeg.ca>

=head1 SEE ALSO

L<Math::Pari> and L<PDL>

=head1 COPYRIGHT

The algorithm used in this module to calculate the Fourier
transforms is based on the C routine of fft4g.c available
at http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html, which is
copyrighted 1996-99 by Takuya OOURA. The file arrays.c included 
here to handle passing arrays to and from C comes from the PGPLOT 
module of Karl Glazebrook <kgb@aaoepp.aao.gov.au>. The perl code 
of Math::FFT is copyright 2000,2005 by Randy Kobes <r.kobes@uwinnipeg.ca>,
and is distributed under the same terms as Perl itself.

=cut

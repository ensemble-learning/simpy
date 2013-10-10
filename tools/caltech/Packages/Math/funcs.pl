$DEBUG = 0;

######### help funcs

sub ok ($$) {
    my($number, $result) = @_ ;

    print "ok $number\n"     if $result ;
    print "not ok $number\n" if !$result ;
}

sub random_matrix ($)
{
    my ($size) = @_;
    my $M = Math::MatrixReal->new($size, $size);
    for (my $i=1; $i<=$size; $i++)
    {
        for (my $j=1; $j<=$size; $j++)
        {
            $M->assign($i,$j,rand());
        }
    }
    return $M;
}

sub ok_matrix ($$$)
{
    my ($no, $M1, $M2) = @_;
    my $tmp = $M1->shadow();
    $tmp->subtract($M1,$M2);
    my $v = $tmp->norm_one();
    ok($no, ($v < 1e-8));
    print " ($no: |Delta| = $v)\n" if $DEBUG;
}

sub ok_matrix_orthogonal ($$)
{
    my ($no, $M) = @_;
    my $tmp = $M->shadow();
    $tmp->one();
    my $transp = $M->shadow();
    $transp->transpose($M);
    $tmp->subtract($M->multiply($transp), $tmp);
    my $v = $tmp->norm_one();
    ok($no, ($v < 1e-8));
    print " ($no: |M * ~M - I| = $v)\n" if $DEBUG;
}

sub ok_eigenvectors ($$$$)
{
    my ($no, $M, $L, $V) = @_;
    # Now check that all of them correspond to eigenvalue * eigenvector
    my ($rows, $columns) = $M->dim();
    unless ($rows == $columns) {
        ok("$no", 0);
        return;
    }
    # Computes the result of all eigenvectors...
    my $test = $M * $V;
    my $test2 = $V->clone();
    for (my $i = 1; $i <= $columns; $i++)
    {
        my $lambda = $L->element($i,1);
        for (my $j = 1; $j <= $rows; $j++)
        { # Compute new vector via lambda * x
            $test2->assign($j, $i, $lambda * $test2->element($j, $i));
        }
      }
    ok_matrix("$no",$test,$test2);
    return;
}
sub similar($$$) {
    my ($x,$y) = @_;
    my $eps = shift || 1e-8;
    abs($x - $y) < $eps ? return 1 : return 0;
}
1;


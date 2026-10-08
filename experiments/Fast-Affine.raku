
use Math::NIntegrate::Utilities;
use Math::NIntegrate::VariableTransformer;
use Math::NIntegrate::VariableTransformer::Affine;
use Math::NIntegrate::Region;

use Math::SparseMatrix;
use Math::SparseMatrix::Native;
use Math::SparseMatrix::Native::Utilities;

use Data::TypeSystem;

# Same/similar timings with both types of transforming
my $type = 'dot'; # dot | separate | Affine

# Emulating MonteCarloRule for a 2D integral.
my $nrow = 2;
my $ncol = 1_000;
my $density = 1;
my $nnz = ($nrow * $ncol * $density).Int;
my $seed = 3432;
my $n = 1_000;
my $working-precision = Num;

my $matrix1 = Math::SparseMatrix::Native::CSRStruct.new.random(:$nrow, :$ncol, :$nnz, :$seed);
say (:$matrix1);

my $tstart = now;
given $type {
    when 'separate' {
        for ^$n {
            my $a1 = numerical(10, $working-precision);
            my $b1 = numerical(15, $working-precision);
            my $a2 = numerical(3, $working-precision);
            my $b2 = numerical(5, $working-precision);

            my $length1 = $b1 - $a1;
            my $col1 = $matrix1.row-at(0);
            $col1.multiply($length1).add($a1);

            my $length2 = $b2 - $a2;
            my $col2 = $matrix1.row-at(1);
            $col2.multiply($length2).add($a2);
        }
    }

    when 'dot' {
        for ^$n {
            my $a1 = numerical(10, $working-precision);
            my $b1 = numerical(15, $working-precision);
            my $a2 = numerical(3, $working-precision);
            my $b2 = numerical(5, $working-precision);

            my $eye = Math::SparseMatrix::Native::CSRStruct.new(dense-matrix => [[$b1 - $a1, 0], [0, $b2 - $a2]]);

            my $scaled1 = $eye.dot($matrix1);
            my $offset1 = Math::SparseMatrix::Native::CSRStruct.new(dense-matrix => [$a1 xx $ncol, $a2 xx $ncol]>>.Array.Array);

            my $res1 = $scaled1.add($offset1);
            #say (:$res1);
        }
    }

    when 'Affine' {

        my $region = Math::NIntegrate::Region.new(min => [10, 3], max => [15, 5], dimension => 2);

        my @min-original-bounds = 10, 3;
        my @max-original-bounds = 15, 5;
        my @min-transform-bounds = 0 xx $region.dimension;
        my @max-transform-bounds = 1 xx $region.dimension;

        my %args =
                region => $region,
                working-precision => Rat,
                ;

        my $vt = Math::NIntegrate::VariableTransformer::Affine.new(
                :@min-original-bounds, :@max-original-bounds, :@min-transform-bounds, :@max-transform-bounds,
                |%args
                );

        my @points = |$matrix1.transpose.Array;
        say deduce-type(@points);
        for ^$n {
            for @points -> @point {
                my $p = $vt.transform(:@point, jacobian => 1, :!functional-bounds);
            }
        }

    }
}
my $tend = now;

say "Total time : { $tend - $tstart }";
say "Mean time  : { ($tend - $tstart) / $n }"

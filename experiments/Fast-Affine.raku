#!/usr/bin/env raku
use v6.d;

my %*SUB-MAIN-OPTS = :named-anywhere;

use Math::NIntegrate::Utilities;
use Math::NIntegrate::VariableTransformer;
use Math::NIntegrate::VariableTransformer::Affine;
use Math::NIntegrate::VariableTransformer::AffineEnBloc;
use Math::NIntegrate::Region;

use Math::SparseMatrix;
use Math::SparseMatrix::Native;
use Math::SparseMatrix::Native::Utilities;

use Data::TypeSystem;


sub MAIN(
        Str:D $type = 'dot', #= One of: "dot", "separate", "Affine", "AffineEnBloc"
        Int:D :$n = 100,        #= Number of benchmarking computation runs,
        Int:D :$ncol = 1_000,   #= Number of points (columns of the representation matrix)
        Int:D :$nrow = 2,       #= Dimension of the points (rows of the representation matrix)
        Int:D :$seed = 3234,    #= Random seed
        Bool:D :$echo = False,  #= Should intermediate results be printed out or not?
         ) {

    # Emulating MonteCarloRule for a 2D integral.
    my $density = 1;
    my $nnz = ($nrow * $ncol * $density).Int;
    my $working-precision = Num;

    my $matrix1 = Math::SparseMatrix::Native::CSRStruct.new.random(:$nrow, :$ncol, :$nnz, :$seed);
    say (:$matrix1) if $echo;

    my ($a1, $b1) = 10, 15;
    my ($a2, $b2) = 3, 5;

    my @result;
    my $tstart = now;
    given $type {
        when 'separate' {
            # Same/similar timings with both "dot" and "separate" if no row-binding for "separate".
            # This only works for 2D
            my $res1;
            for ^$n {
                $a1 = numerical($a1, $working-precision);
                $b1 = numerical($b1, $working-precision);
                $a2 = numerical($a2, $working-precision);
                $b2 = numerical($b2, $working-precision);

                my $length1 = $b1 - $a1;
                my $col1 = $matrix1.row-at(0);
                $col1.multiply($length1).add($a1);

                my $length2 = $b2 - $a2;
                my $col2 = $matrix1.row-at(1);
                $col2.multiply($length2).add($a2);

                $res1 = $col1.clone.row-bind($col2)
            }
            say ($res1) if $echo;
            @result = |$res1.transpose.Array;
        }

        when 'dot' {
            my $res1;
            for ^$n {
                $a1 = numerical($a1, $working-precision);
                $b1 = numerical($b1, $working-precision);
                $a2 = numerical($a2, $working-precision);
                $b2 = numerical($b2, $working-precision);

                my $eye = Math::SparseMatrix::Native::CSRStruct.new(dense-matrix => [[$b1 - $a1, 0], [0, $b2 - $a2]]);

                my $scaled1 = $eye.dot($matrix1);
                my $offset1 = Math::SparseMatrix::Native::CSRStruct.new(dense-matrix => [$a1 xx $ncol, $a2 xx $ncol]>>.Array.Array);

                $res1 = $scaled1.add($offset1);
            }
            say ($res1) if $echo;
            # Pretty print the sparse matrix
            #if $echo {
            #   Math::SparseMatrix.new(matrix =>  Math::SparseMatrix::NativeAdapter.new($res1)).print
            #}
            @result = |$res1.transpose.Array;
        }

        when 'Affine' {
            my $region = Math::NIntegrate::Region.new(min => [$a1, $a2], max => [$b1, $b2], dimension => 2);

            my @min-transform-bounds = $a1, $a2;
            my @max-transform-bounds = $b1, $b2;
            my @min-original-bounds = 0 xx $region.dimension;
            my @max-original-bounds = 1 xx $region.dimension;

            my %args =
                    region => $region,
                    working-precision => Rat,
                    ;

            my $vt = Math::NIntegrate::VariableTransformer::Affine.new(
                    :@min-original-bounds, :@max-original-bounds, :@min-transform-bounds, :@max-transform-bounds,
                    |%args
                    );

            my @points = |$matrix1.transpose.Array;
            say deduce-type(@points) if $echo;
            for ^$n {
                @result = [];
                for @points -> @point {
                    my $p = $vt.transform(:@point, jacobian => 1, :!functional-bounds);
                    @result.push($p<point>)
                }
            }
        }

        when 'AffineEnBloc' {
            my $region = Math::NIntegrate::Region.new(min => [$a1, $a2], max => [$b1, $b2], dimension => 2);

            my @min-transform-bounds = $a1, $a2;
            my @max-transform-bounds = $b1, $b2;
            my @min-original-bounds = 0 xx $region.dimension;
            my @max-original-bounds = 1 xx $region.dimension;

            my %args =
                    region => $region,
                    working-precision => Rat,
                    ;

            my $vt = Math::NIntegrate::VariableTransformer::AffineEnBloc.new(
                    :@min-original-bounds, :@max-original-bounds, :@min-transform-bounds, :@max-transform-bounds,
                    |%args
                    );

            my @points = |$matrix1.transpose.Array;
            say deduce-type(@points) if $echo;
            for ^$n {
                @result = |$vt.transform(:@points, jacobian => 1, :!functional-bounds)<points>;
            }
        }
    }
    my $tend = now;

    say "Total time : { $tend - $tstart }";
    say "Mean time  : { ($tend - $tstart) / $n }";

    if $echo { .say for @result }
}
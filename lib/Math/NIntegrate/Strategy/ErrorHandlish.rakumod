use v6.d;

#| Error message and no-result role
role Math::NIntegrate::Strategy::ErrorHandlish {

    my %.no-result = integral => Whatever, error => Whatever;

    has $.msgProbOrig = 'Problems integrating the original set of regions.';
    has $.msgProbSplit = 'Problems integrating the split regions.';
    has $.msgProbIMT = 'Problems integrating the region after IMT application.';
    has $.msgProbRegion = 'Problems integrating the region at level';
    has $.msgSuppress = 'Suppressing further messages of this type.';

    method max-iterative-refinements-reached(@min, @max, $max-recursion) {
        my @range-info = |(@min>>.Numeric Z @max>>.Numeric).map({ "[{$_.head}, {$_.tail}]" });
        "Failed to converge to prescribed accuracy after {$max-recursion} " ~
                "recursive refinement{ $max-recursion == 1 ?? '' !! 's'} " ~
                "in the region {@range-info.join(' x ')}."
    }

    method max-recursion-reached(@args, $max-recursion) {
        "Failed to converge to prescribed accuracy after $max-recursion recursive bisections in near {@args.raku}."
    }
}

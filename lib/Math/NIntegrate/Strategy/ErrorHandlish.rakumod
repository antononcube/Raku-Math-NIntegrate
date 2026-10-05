use v6.d;

#| Error message and no-result role
role Math::NIntegrate::Strategy::ErrorHandlish {

    my %.no-result = integral => Whatever, error => Whatever;

    has $.msgProbOrig = 'Problems integrating the original set of regions.';
    has $.msgProbSplit = 'Problems integrating the split regions.';
    has $.msgProbIMT = 'Problems integrating the region after IMT application.';
    has $.msgProbRegion = 'Problems integrating the region at level';
    has $.msgSuppress = 'Suppressing further messages of this type.'
}

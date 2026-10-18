use v6.d;

use Math::NIntegrate::Strategy::MonteCarlo;

# This class does not need to exist since Builder handles the creation of a
# ::Strategy::MonteCarlo object with a MonteCarloRule with a low-discrepancy point generator.

# Inheriting from- or delegating to MonteCarlo ?
class Math::NIntegrate::Strategy::QuasiMonteCarlo
        is Math::NIntegrate::Strategy::MonteCarlo {

    submethod TWEAK(:$random-seed, *%args) {

        note 'The value of $random-seed is ignored by the QuasiMonteCarlo strategy.'
        with $random-seed;
    }
}

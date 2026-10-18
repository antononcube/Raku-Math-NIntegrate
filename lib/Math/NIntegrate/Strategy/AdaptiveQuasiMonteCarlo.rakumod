use v6.d;

use Math::NIntegrate::Strategy::AdaptiveMonteCarlo;

# Very similar to QuasiMonteCarlo
class Math::NIntegrate::Strategy::AdaptiveQuasiMonteCarlo
        is Math::NIntegrate::Strategy::AdaptiveMonteCarlo {

    submethod TWEAK(:$random-seed, *%args) {

        note 'The value of $random-seed is ignored by the AdaptiveQuasiMonteCarlo strategy.'
        with $random-seed;
    }
}

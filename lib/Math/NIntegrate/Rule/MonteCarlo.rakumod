use v6.d;

use Math::NIntegrate::Utilities;
use Math::NIntegrate::Rule::General;

class Math::NIntegrate::Rule::MonteCarlo
        is Math::NIntegrate::Rule::General {

    has $.points is required;
    has &.point-generator = WhateverCode;
    has &.axis-selector = WhateverCode;

    submethod BUILD(UInt:D :$!points, :&!point-generator = WhateverCode, :$working-precision = Num) {
        # In general, working precision of Monte-Carlo methods does not matter --
        # the precision- and accuracy goals are low.
        # But since all the rule classes take working-precision as an argument it is also present here.
        without &!point-generator {
            # The most universal point generator has very fine granularity: one point per invocation.
            &!point-generator = -> UInt:D $n, UInt:D $axis, UInt:D $dim, UInt:D $points-per-step {1.rand}
        }
    }

    submethod TWEAK(:%args) {
        # Why is this needed?
        # It passifies Math::NIntegrate::Rule::General workings.
        # For now Math::NIntegrate::Region takes Math::NIntegrate::Rule::General not Math::NIntegrate::Rule.
        # Hence Math::NIntegrate::Rule::MonteCarlo is a Math::NIntegrate::Rule::General.
        self.abscissas = |(&!point-generator xx $!points);
        self.weights = |((1/$!points) xx $!points);
        self.error-weights = |((1/$!points) xx $!points);
    }

    multi method new($points, $working-precision = Num, :&point-generator = WhateverCode) {
        self.bless(:$points, :&point-generator, :$working-precision)
    }

    multi method new($points, :prec(:$working-precision) = Num, :&point-generator = WhateverCode) {
        self.bless(:$points, :&point-generator, :$working-precision)
    }

    multi method new(:$points, :$working-precision = Num, :&point-generator = WhateverCode) {
        self.bless(:$points, :$working-precision)
    }

    method integrate($region) {

        my ($sum, $sqsum, $n);
        $sum = $region.reuse-values<sum> // 0;
        $sqsum = $region.reuse-values<sqsum> // 0;
        $n = $region.reuse-values<n> // 0;

        self.abscissas = [];
        my $dim = $region.dimension;
        for ^$!points {
            $n++;
            my @point = (^$dim).map({ &!point-generator($n, $_, $dim, $!points) });

            # To be used by the axis selector
            self.abscissas.push(@point);

            my $value = $region.eval-integrand(@point);
            $sum += $value;
            $sqsum += $value ** 2;
        }

        $region.reuse-values = %(:$sum, :$sqsum, :$n);

        self.integral = $sum / $n;
        self.error = sqrt( ($sqsum / $n - ($sum / $n) ** 2) / $n );
        self.largest-error-axis = &!axis-selector ?? &!axis-selector(self.abscissas) !! 0;

        return self;
    }
}

use v6.d;

use Math::NIntegrate::Utilities;
use Math::NIntegrate::Rule::General;

class Math::NIntegrate::Rule::MonteCarlo
        is Math::NIntegrate::Rule::General {

    has $.points is required;
    has &.point-generator = WhateverCode;
    has &.axis-selector = WhateverCode;

    #======================================================
    # Creators
    #======================================================

    submethod BUILD(UInt:D :$!points, UInt:D :$dimension!, :&!point-generator = WhateverCode, :$working-precision = Num) {
        # In general, working precision of Monte-Carlo methods does not matter --
        # the precision- and accuracy goals are low.
        # But since all the rule classes take working-precision as an argument it is also present here.
        without &!point-generator {
            # The most universal point generator has very fine granularity: one point coordinate per invocation.
            &!point-generator = -> UInt:D $n, UInt:D $axis, UInt:D $dim, UInt:D $points-per-step {1.rand}
        }

        self.dimension = $dimension;

        self.abscissas = [];
        for ^$!points -> $n {
            my @point = (^self.dimension).map({ &!point-generator($n, $_, self.dimension, $!points) });
            self.abscissas.push(@point);
        }

        self.weights = |((1/$!points) xx $!points);
        self.error-weights = |((1/$!points) xx $!points);
    }

    multi method new($points, $working-precision = Num, UInt:D :$dimension = 1, :&point-generator = WhateverCode) {
        self.bless(:$points, :$dimension, :&point-generator, :$working-precision)
    }

    multi method new($points, UInt:D :$dimension = 1, :prec(:$working-precision) = Num, :&point-generator = WhateverCode) {
        self.bless(:$points, :$dimension, :$working-precision, :&point-generator)
    }

    multi method new(:$points, UInt:D :$dimension, :$working-precision = Num, :&point-generator = WhateverCode) {
        self.bless(:$points, :$dimension, :$working-precision, :&point-generator)
    }

    #| Copy attributes from an object
    method copy(Math::NIntegrate::Rule::MonteCarlo:D $from, Bool:D :deep(:deep-copy(:$clone)) = False) {
        self.Math::NIntegrate::Rule::General::copy($from, :$clone);
        $!points = $from.points;
        &!point-generator = $from.point-generator;
        &!axis-selector = $from.axis-selector;

        return self
    }

    #| Clone the object
    method clone(-->Math::NIntegrate::Rule::MonteCarlo:D) {
        my $obj = Math::NIntegrate::Rule::MonteCarlo.new(:$!points, dimension => self.dimension, :&!point-generator);
        say (:$obj);
        $obj.copy(self, :clone)
    }

    #======================================================
    # Integration
    #======================================================

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

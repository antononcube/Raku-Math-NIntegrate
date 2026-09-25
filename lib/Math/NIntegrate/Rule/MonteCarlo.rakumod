use v6.d;

use Math::NIntegrate::Utilities;
use Math::NIntegrate::Rule::General;

class Math::NIntegrate::Rule::MonteCarlo
        is Math::NIntegrate::Rule::General {

    has $.points is required;
    has &.point-generator = WhateverCode;

    submethod BUILD(UInt:D :$!points, :&point-generator = WhateverCode, :$working-precision = Num) {
        my %res = self.make-weights($!points, $working-precision);
        self.abscissas = |%res<abscissas>;
        self.weights = |%res<weights>;
        self.error-weights = |%res<error-weights>;
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

    method make-weights(UInt:D $points, $working-precision = Num) {!!!}
}

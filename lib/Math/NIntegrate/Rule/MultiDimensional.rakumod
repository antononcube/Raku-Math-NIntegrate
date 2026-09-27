use v6.d;

use Math::NIntegrate::Rule::General;
use Math::NIntegrate::Utilities;

class Math::NIntegrate::Rule::MultiDimensional
        is Math::NIntegrate::Rule::General {

    has $.generators;
    has @.rule-points;
    has @.scales;
    has @.norms;

    #======================================================
    # Creators
    #======================================================

    submethod TWEAK(:$!generators, :$dimension!, :$working-precision = Num) {
        # The dimension attribute is inherited from Rule::General, so assign it before
        # creating dimension-dependent DCUHRE rule data.
        self.dimension = $dimension;

        my %data = do given $!generators {
            when $_ == 7 {
                self.d07hre(self.dimension) unless %!data.elems;
            }
            when $_ == 9 {
                self.d09hre(self.dimension) unless %!data.elems;
            }
            default {
                die "No multidimensional rule is available with $!generators generators.";
            }
        }
        self.abscissas = %data<generators>;
        self.weights = %data<weights>;
        self.error-weights = %data<error-weights>;
        @!rule-points = %data<rule-points>;

        # Fill-in the scales and norms -- see DEINHR
        for ^3 -> $k {
            for ^$!generators -> $i {
                my @we;

                if !is-zero(self.weights[$k + 1][$i]) {
                    @!scales[$k][$i] = - self.weights[$k + 2][$i] / self.weights[$k + 1][$i]
                } else {
                    @!scales[$k][$i] = 100
                }

                for ^$!generators -> $j {
                    @we[$j] = self.weights[$k + 2][$j] + @!scales[$k][$i] * @!scales[$k + 1][$j]
                }

                @!norms[$k][$i] = 0;

                for ^$!generators -> $j {
                    @!norms[$k][$i] += @!rule-points[$j] * @we[$j].abs
                }

                @!norms[$k][$i] = 2 ** self.dimension / @!norms[$k][$i]
            }
        }
    }

    multi method new($dimension, :$generators = 7, :prec(:$working-precision) = Num) {
        self.bless(:$generators, :$dimension, :$working-precision)
    }

    multi method new(:$generators = 7, :dim(:$dimension)!, :prec(:$working-precision) = Num) {
        self.bless(:$generators, :$dimension, :$working-precision)
    }

    #| Copy from object
    method copy(Math::NIntegrate::Rule::MultiDimensional:D $from,  Bool:D :deep(:deep-copy(:$clone)) = False) {
        # Delegate to parent class
        self.Math::NIntegrate::Rule::General::copy($from, :$clone);

        $!generators = $from.generators;
        # The rule data does not mutate, hence, cloning it is not needed.
        # But it is small enough, so, doing for consistency.
        @!rule-points = $clone ?? $from.rule-points !! $from.rule-points.clone;
        @!scales = $clone ?? $from.scales !! $from.scales.clone;
        @!norms = $clone ?? $from.norms !! $from.norms.clone;

        return self
    }

    method clone(-->Math::NIntegrate::Rule::MultiDimensional) {
        Math::NIntegrate::Rule::MultiDimensional.new(dimension => $.dimension, generators => $!generators).copy(self, :clone)
    }

    #======================================================
    # Data methods
    #======================================================

    #| Initialise the DCUHRE degree-seven fully symmetric rule.
    #|
    #| The returned arrays retain the layout of the original D07HRE routine:
    #| `weights[rule][generator]` and `generators[axis][generator]`. Rule zero
    #| is the degree-seven integration rule; rules one through four are its
    #| embedded null rules. `rule-points` gives each generator orbit's size.
    method d07hre(
            UInt:D $dimension,
            --> Map:D
                  ) {
        die 'D07HRE requires a positive dimension.' unless $dimension > 1;

        my UInt:D $weight-length = 6;

        my @weights = (^5).map({ [0 xx $weight-length] }).Array;
        my @generators = (^$dimension).map({ [0 xx $weight-length] }).Array;
        my @rule-points = (2 * $dimension) xx $weight-length;

        my Num $two-to-dimension = 2 ** $dimension;
        @rule-points[$weight-length - 1] = $two-to-dimension;
        @rule-points[$weight-length - 2] = 2 * $dimension * ($dimension - 1);
        @rule-points[0] = 1;

        # Squared generator parameters.
        my Num $lambda0 = 0.4707e0;
        my Num $lambda-p = 0.5625e0;
        my Num $lambda1 = 4 / (15 - 5 / $lambda0);
        my Num $ratio = (1 - $lambda1 / $lambda0) / 27;
        my Num $lambda2 = (5 - 7 * $lambda1 - 35 * $ratio) / (7 - 35 * $lambda1 / 3 - 35 * $ratio / $lambda0);

        # Degree-seven rule weights.
        @weights[0][5] = 1e0 / (3e0 * $lambda0) ** 3 / $two-to-dimension;
        @weights[0][4] = (1e0 - 5e0 * $lambda0 / 3e0) / (60e0 * ($lambda1 - $lambda0) * $lambda1 ** 2);
        @weights[0][2] = (1e0 - 5e0 * $lambda2 / 3e0
                - 5e0 * $two-to-dimension * @weights[0][5] * $lambda0 * ($lambda0 - $lambda2))
                / (10e0 * $lambda1 * ($lambda1 - $lambda2))
                - 2e0 * ($dimension - 1) * @weights[0][4];
        @weights[0][1] = (1e0 - 5e0 * $lambda1 / 3e0
                - 5e0 * $two-to-dimension * @weights[0][5] * $lambda0 * ($lambda0 - $lambda1))
                / (10e0 * $lambda2 * ($lambda2 - $lambda1));

        # The lower-degree rules are subsequently converted into null-rule
        # weights by subtracting the degree-seven weights below.
        @weights[1][5] = 1e0 / (36e0 * $lambda0 ** 3) / $two-to-dimension;
        @weights[1][4] = (1e0 - 9e0 * $two-to-dimension * @weights[1][5] * $lambda0 ** 2)
                / (36e0 * $lambda1 ** 2);
        @weights[1][2] = (1e0 - 5e0 * $lambda2 / 3e0
                - 5e0 * $two-to-dimension * @weights[1][5] * $lambda0 * ($lambda0 - $lambda2))
                / (10e0 * $lambda1 * ($lambda1 - $lambda2))
                - 2e0 * ($dimension - 1) * @weights[1][4];
        @weights[1][1] = (1e0 - 5e0 * $lambda1 / 3e0
                - 5e0 * $two-to-dimension * @weights[1][5] * $lambda0 * ($lambda0 - $lambda1))
                / (10e0 * $lambda2 * ($lambda2 - $lambda1));

        @weights[2][5] = 5e0 / (108e0 * $lambda0 ** 3) / $two-to-dimension;
        @weights[2][4] = (1e0 - 9e0 * $two-to-dimension * @weights[2][5] * $lambda0 ** 2)
                / (36e0 * $lambda1 ** 2);
        @weights[2][2] = (1e0 - 5e0 * $lambda-p / 3e0
                - 5e0 * $two-to-dimension * @weights[2][5] * $lambda0 * ($lambda0 - $lambda-p))
                / (10e0 * $lambda1 * ($lambda1 - $lambda-p))
                - 2e0 * ($dimension - 1) * @weights[2][4];
        @weights[2][3] = (1e0 - 5e0 * $lambda1 / 3e0
                - 5e0 * $two-to-dimension * @weights[2][5] * $lambda0 * ($lambda0 - $lambda1))
                / (10e0 * $lambda-p * ($lambda-p - $lambda1));

        @weights[3][5] = 1e0 / (54e0 * $lambda0 ** 3) / $two-to-dimension;
        @weights[3][4] = (1e0 - 18e0 * $two-to-dimension * @weights[3][5] * $lambda0 ** 2)
                / (72e0 * $lambda1 ** 2);
        @weights[3][2] = (1e0 - 10e0 * $lambda2 / 3e0
                - 10e0 * $two-to-dimension * @weights[3][5] * $lambda0 * ($lambda0 - $lambda2))
                / (20e0 * $lambda1 * ($lambda1 - $lambda2))
                - 2e0 * ($dimension - 1) * @weights[3][4];
        @weights[3][1] = (1e0 - 10e0 * $lambda1 / 3e0
                - 10e0 * $two-to-dimension * @weights[3][5] * $lambda0 * ($lambda0 - $lambda1))
                / (20e0 * $lambda2 * ($lambda2 - $lambda1));

        # Generator values are the positive square roots of their parameters.
        $lambda0 = $lambda0.sqrt;
        $lambda1 = $lambda1.sqrt;
        $lambda2 = $lambda2.sqrt;
        $lambda-p = $lambda-p.sqrt;
        for ^$dimension -> $axis {
            @generators[$axis][5] = $lambda0;
        }
        @generators[0][4] = $lambda1;
        @generators[1][4] = $lambda1;
        @generators[0][1] = $lambda2;
        @generators[0][2] = $lambda1;
        @generators[0][3] = $lambda-p;

        # Convert embedded-rule weights to null-rule weights, then scale the
        # integration-rule weights for the [-1, 1]^dimension hypercube.
        @weights[0][0] = $two-to-dimension;
        for 1 ..^ 5 -> $rule {
            for 1 ..^ $weight-length -> $generator {
                @weights[$rule][$generator] -= @weights[0][$generator];
                @weights[$rule][0] -= @rule-points[$generator] * @weights[$rule][$generator];
            }
        }
        for 1 ..^ $weight-length -> $generator {
            @weights[0][$generator] *= $two-to-dimension;
            @weights[0][0] -= @rule-points[$generator] * @weights[0][$generator];
        }

        my @error-coefficients = 5, 5, 1, 5, 0.5, 0.25;

        # Make the generator points to be in [-1/2, 1/2]^dim
        @generators = @generators.map({ $_ <<*>> 0.5 });

        return {
            :@weights,
            :@generators,
            :@error-coefficients,
            :@rule-points,
        };
    }

    #| Initialise the DCUHRE degree-nine fully symmetric rule.
    #|
    #| This is a direct translation of D09HRE.  As for `d07hre`, rule zero
    #| is the integration rule and rules one through four are null rules.
    method d09hre(
            UInt:D $dimension,
            --> Map:D
                  ) {
        die 'D09HRE requires a positive dimension.' unless $dimension > 0;

        my UInt:D $weight-length = 9;


        my @weights = (^5).map({ [0e0 xx $weight-length] }).Array;
        my @generators = (^$dimension).map({ [0e0 xx $weight-length] }).Array;
        my @rule-points = (2e0 * $dimension) xx $weight-length;

        my Num $two-to-dimension = 2e0 ** $dimension;
        @rule-points[8] = $two-to-dimension;
        @rule-points[7] = 4e0 * $dimension * ($dimension - 1) * ($dimension - 2) / 3e0
        if $dimension > 2;
        @rule-points[6] = 4e0 * $dimension * ($dimension - 1);
        @rule-points[5] = 2e0 * $dimension * ($dimension - 1);
        @rule-points[0] = 1e0;

        # Squared generator parameters.
        my Num $lambda0 = 0.4707e0;
        my Num $lambda1 = 4e0 / (15e0 - 5e0 / $lambda0);
        my Num $ratio = (1e0 - $lambda1 / $lambda0) / 27e0;
        my Num $lambda2 = (5e0 - 7e0 * $lambda1 - 35e0 * $ratio)
                / (7e0 - 35e0 * $lambda1 / 3e0 - 35e0 * $ratio / $lambda0);
        $ratio *= (1e0 - $lambda2 / $lambda0) / 3e0;
        my Num $lambda3 = (7e0 - 9e0 * ($lambda2 + $lambda1)
                + 63e0 * $lambda2 * $lambda1 / 5e0 - 63e0 * $ratio)
                / (9e0 - 63e0 * ($lambda2 + $lambda1) / 5e0
                        + 21e0 * $lambda2 * $lambda1 - 63e0 * $ratio / $lambda0);
        my Num $lambda-p = 0.0625e0;

        # Degree-nine rule weights.
        @weights[0][8] = 1e0 / (3e0 * $lambda0) ** 4 / $two-to-dimension;
        @weights[0][7] = (1e0 - 1e0 / (3e0 * $lambda0)) / (6e0 * $lambda1) ** 3
        if $dimension > 2;
        @weights[0][6] = (1e0 - 7e0 * ($lambda0 + $lambda1) / 5e0
                + 7e0 * $lambda0 * $lambda1 / 3e0)
                / (84e0 * $lambda1 * $lambda2 * ($lambda2 - $lambda0)
                        * ($lambda2 - $lambda1));
        @weights[0][5] = (1e0 - 7e0 * ($lambda0 + $lambda2) / 5e0
                + 7e0 * $lambda0 * $lambda2 / 3e0)
                / (84e0 * $lambda1 ** 2 * ($lambda1 - $lambda0)
                        * ($lambda1 - $lambda2))
                - @weights[0][6] * $lambda2 / $lambda1
                - 2e0 * ($dimension - 2) * @weights[0][7];
        @weights[0][3] = (1e0 - 9e0 * (($lambda0 + $lambda1 + $lambda2) / 7e0
                - ($lambda0 * $lambda1 + $lambda0 * $lambda2 + $lambda1 * $lambda2) / 5e0)
                - 3e0 * $lambda0 * $lambda1 * $lambda2)
                / (18e0 * $lambda3 * ($lambda3 - $lambda0) * ($lambda3 - $lambda1)
                        * ($lambda3 - $lambda2));
        @weights[0][2] = (1e0 - 9e0 * (($lambda0 + $lambda1 + $lambda3) / 7e0
                - ($lambda0 * $lambda1 + $lambda0 * $lambda3 + $lambda1 * $lambda3) / 5e0)
                - 3e0 * $lambda0 * $lambda1 * $lambda3)
                / (18e0 * $lambda2 * ($lambda2 - $lambda0) * ($lambda2 - $lambda1)
                        * ($lambda2 - $lambda3)) - 2e0 * ($dimension - 1) * @weights[0][6];
        @weights[0][1] = (1e0 - 9e0 * (($lambda0 + $lambda2 + $lambda3) / 7e0
                - ($lambda0 * $lambda2 + $lambda0 * $lambda3 + $lambda2 * $lambda3) / 5e0)
                - 3e0 * $lambda0 * $lambda2 * $lambda3)
                / (18e0 * $lambda1 * ($lambda1 - $lambda0) * ($lambda1 - $lambda2)
                        * ($lambda1 - $lambda3))
                - 2e0 * ($dimension - 1) * (@weights[0][6] + @weights[0][5]
                        + ($dimension - 2) * @weights[0][7]);

        # Two degree-seven, one degree-five, and one degree-three rules.
        @weights[1][8] = 1e0 / (108e0 * $lambda0 ** 4) / $two-to-dimension;
        @weights[2][8] = 5e0 / (324e0 * $lambda0 ** 4) / $two-to-dimension;
        @weights[3][8] = 2e0 / (81e0 * $lambda0 ** 4) / $two-to-dimension;
        for 1 .. 3 -> $rule {
            @weights[$rule][7] = ((($rule == 3 ?? 2e0 !! 1e0)
                    - 27e0 * $two-to-dimension * @weights[$rule][8] * $lambda0 ** 3)
                    / (6e0 * $lambda1) ** 3) if $dimension > 2;
        }

        for 1, 2, 3 -> $rule {
            my $factor = $rule == 3 ?? 2e0 !! 1e0;
            @weights[$rule][6] = ($factor - 5e0 * $lambda1 / 3e0
                    - 15e0 * $two-to-dimension * @weights[$rule][8] * $lambda0
                            ** ($rule == 3 ?? 1 !! 2) * ($lambda0 - $lambda1))
                    / (60e0 * $lambda1 * $lambda2 * ($lambda2 - $lambda1));
            @weights[$rule][5] = (1e0 - 9e0 * (8e0 * $lambda1 * $lambda2
                    * @weights[$rule][6] + $two-to-dimension * @weights[$rule][8]
            * $lambda0 ** 2)) / (36e0 * $lambda1 ** 2)
                    - 2e0 * @weights[$rule][7] * ($dimension - 2);
        }

        @weights[1][3] = (1e0 - 7e0 * (($lambda1 + $lambda2) / 5e0
                - $lambda1 * $lambda2 / 3e0 + $two-to-dimension * @weights[1][8]
                * $lambda0 * ($lambda0 - $lambda1) * ($lambda0 - $lambda2)))
                / (14e0 * $lambda3 * ($lambda3 - $lambda1) * ($lambda3 - $lambda2));
        @weights[1][2] = (1e0 - 7e0 * (($lambda1 + $lambda3) / 5e0
                - $lambda1 * $lambda3 / 3e0 + $two-to-dimension * @weights[1][8]
                * $lambda0 * ($lambda0 - $lambda1) * ($lambda0 - $lambda3)))
                / (14e0 * $lambda2 * ($lambda2 - $lambda1) * ($lambda2 - $lambda3))
                - 2e0 * ($dimension - 1) * @weights[1][6];
        @weights[1][1] = (1e0 - 7e0 * (($lambda2 + $lambda3) / 5e0
                - $lambda2 * $lambda3 / 3e0 + $two-to-dimension * @weights[1][8]
                * $lambda0 * ($lambda0 - $lambda2) * ($lambda0 - $lambda3)))
                / (14e0 * $lambda1 * ($lambda1 - $lambda2) * ($lambda1 - $lambda3))
                - 2e0 * ($dimension - 1) * (@weights[1][6] + @weights[1][5]
                        + ($dimension - 2) * @weights[1][7]);

        @weights[2][4] = (1e0 - 7e0 * (($lambda1 + $lambda2) / 5e0
                - $lambda1 * $lambda2 / 3e0 + $two-to-dimension * @weights[2][8]
                * $lambda0 * ($lambda0 - $lambda1) * ($lambda0 - $lambda2)))
                / (14e0 * $lambda-p * ($lambda-p - $lambda1) * ($lambda-p - $lambda2));
        @weights[2][2] = (1e0 - 7e0 * (($lambda1 + $lambda-p) / 5e0
                - $lambda1 * $lambda-p / 3e0 + $two-to-dimension * @weights[2][8]
                * $lambda0 * ($lambda0 - $lambda1) * ($lambda0 - $lambda-p)))
                / (14e0 * $lambda2 * ($lambda2 - $lambda1) * ($lambda2 - $lambda-p))
                - 2e0 * ($dimension - 1) * @weights[2][6];
        @weights[2][1] = (1e0 - 7e0 * (($lambda2 + $lambda-p) / 5e0
                - $lambda2 * $lambda-p / 3e0 + $two-to-dimension * @weights[2][8]
                * $lambda0 * ($lambda0 - $lambda2) * ($lambda0 - $lambda-p)))
                / (14e0 * $lambda1 * ($lambda1 - $lambda2) * ($lambda1 - $lambda-p))
                - 2e0 * ($dimension - 1) * (@weights[2][6] + @weights[2][5]
                        + ($dimension - 2) * @weights[2][7]);

        @weights[3][3] = (2e0 - 7e0 * (($lambda1 + $lambda2) / 5e0
                - $lambda1 * $lambda2 / 3e0 + $two-to-dimension * @weights[3][8]
                * $lambda0 * ($lambda0 - $lambda1) * ($lambda0 - $lambda2)))
                / (14e0 * $lambda3 * ($lambda3 - $lambda1) * ($lambda3 - $lambda2));
        @weights[3][2] = (2e0 - 7e0 * (($lambda1 + $lambda3) / 5e0
                - $lambda1 * $lambda3 / 3e0 + $two-to-dimension * @weights[3][8]
                * $lambda0 * ($lambda0 - $lambda1) * ($lambda0 - $lambda3)))
                / (14e0 * $lambda2 * ($lambda2 - $lambda1) * ($lambda2 - $lambda3))
                - 2e0 * ($dimension - 1) * @weights[3][6];
        @weights[3][1] = (2e0 - 7e0 * (($lambda2 + $lambda3) / 5e0
                - $lambda2 * $lambda3 / 3e0 + $two-to-dimension * @weights[3][8]
                * $lambda0 * ($lambda0 - $lambda2) * ($lambda0 - $lambda3)))
                / (14e0 * $lambda1 * ($lambda1 - $lambda2) * ($lambda1 - $lambda3))
                - 2e0 * ($dimension - 1) * (@weights[3][6] + @weights[3][5]
                        + ($dimension - 2) * @weights[3][7]);
        @weights[4][1] = 1e0 / (6e0 * $lambda1);

        # Generator values are the positive square roots of their parameters.
        $lambda0 = $lambda0.sqrt;
        $lambda1 = $lambda1.sqrt;
        $lambda2 = $lambda2.sqrt;
        $lambda3 = $lambda3.sqrt;
        $lambda-p = $lambda-p.sqrt;
        for ^$dimension -> $axis {
            @generators[$axis][8] = $lambda0;
        }
        if $dimension > 2 {
            @generators[0][7] = $lambda1;
            @generators[1][7] = $lambda1;
            @generators[2][7] = $lambda1;
        }
        @generators[0][6] = $lambda1;
        @generators[1][6] = $lambda2 if $dimension > 1;
        @generators[0][5] = $lambda1;
        @generators[1][5] = $lambda1 if $dimension > 1;
        @generators[0][4] = $lambda-p;
        @generators[0][3] = $lambda3;
        @generators[0][2] = $lambda2;
        @generators[0][1] = $lambda1;

        # Convert embedded-rule weights to null-rule weights, then scale the
        # integration-rule weights for the [-1, 1]^dimension hypercube.
        @weights[0][0] = $two-to-dimension;
        for 1 ..^ 5 -> $rule {
            for 1 ..^ $weight-length -> $generator {
                @weights[$rule][$generator] -= @weights[0][$generator];
                @weights[$rule][0] -= @rule-points[$generator] * @weights[$rule][$generator];
            }
        }
        for 1 ..^ $weight-length -> $generator {
            @weights[0][$generator] *= $two-to-dimension;
            @weights[0][0] -= @rule-points[$generator] * @weights[0][$generator];
        }

        my @error-coefficients = 5, 5, 1, 5, 0.5, 0.25;

        # Make the generator points to be in [-1/2, 1/2]^dim
        @generators = @generators.map({ $_ <<*>> 0.5 });

        return {
            :@weights,
            :@generators,
            :@error-coefficients,
            :@rule-points,
        };
    }

    #======================================================
    # Integration helpers
    #======================================================

    method defshr() {!!!}

    method derlhr() {!!!}

    #======================================================
    # Integration
    #======================================================

    #| Apply the selected fully symmetric rule to a region.
    method integrate($region --> Math::NIntegrate::Rule::MultiDimensional:D) {
        die 'The rule and region dimensions must match.'
        unless $region.dimension == self.dimension;


        return self;
    }
}
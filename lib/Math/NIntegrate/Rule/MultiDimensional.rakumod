use v6.d;

use Math::NIntegrate::Rule::General;

class Math::NIntegrate::Rule::MultiDimensional
        is Math::NIntegrate::Rule::General {

    has $.generators;
    has %.data;

    submethod TWEAK(:$!generators, :$dimension!, :$working-precision = Num) {
        # `dimension` is inherited from Rule::General, so assign it before
        # creating dimension-dependent DCUHRE rule data.
        self.dimension = $dimension;
        given $!generators {
            when $_ == 7 {
                %!data = self.d07hre(self.dimension) unless %!data.elems;
            }
            default {
                die "No multidimensional rule is available with $!generators generators.";
            }
        }
    }

    multi method new($dimension, :$generators = 7, :prec(:$working-precision) = Num) {
        self.bless(:$generators, :$dimension, :$working-precision)
    }

    multi method new(:$generators = 7, :dim(:$dimension)!, :prec(:$working-precision) = Num) {
        self.bless(:$generators, :$dimension, :$working-precision)
    }

    #| Initialise the DCUHRE degree-seven fully symmetric rule.
    #|
    #| The returned arrays retain the layout of the original D07HRE routine:
    #| `weights[rule][generator]` and `generators[axis][generator]`. Rule zero
    #| is the degree-seven integration rule; rules one through four are its
    #| embedded null rules. `rule-points` gives each generator orbit's size.
    method d07hre(
            UInt:D $dimension,
            UInt:D :$weight-length = 6
            --> Map:D
                  ) {
        die 'D07HRE requires a positive dimension.' unless $dimension > 0;
        die 'D07HRE requires a weight length of 6.' unless $weight-length == 6;

        my @weights = (^5).map({ [0e0 xx $weight-length] }).Array;
        my @generators = (^$dimension).map({ [0e0 xx $weight-length] }).Array;
        my @rule-points = (2e0 * $dimension) xx $weight-length;

        my Num $two-to-dimension = 2e0 ** $dimension;
        @rule-points[$weight-length - 1] = $two-to-dimension;
        @rule-points[$weight-length - 2] = 2e0 * $dimension * ($dimension - 1);
        @rule-points[0] = 1e0;

        # Squared generator parameters.
        my Num $lambda0 = 0.4707e0;
        my Num $lambda-p = 0.5625e0;
        my Num $lambda1 = 4e0 / (15e0 - 5e0 / $lambda0);
        my Num $ratio = (1e0 - $lambda1 / $lambda0) / 27e0;
        my Num $lambda2 = (5e0 - 7e0 * $lambda1 - 35e0 * $ratio)
                / (7e0 - 35e0 * $lambda1 / 3e0 - 35e0 * $ratio / $lambda0);

        # Degree-seven rule weights.
        @weights[0][5] = 1e0 / (3e0 * $lambda0) ** 3 / $two-to-dimension;
        @weights[0][4] = (1e0 - 5e0 * $lambda0 / 3e0)
                / (60e0 * ($lambda1 - $lambda0) * $lambda1 ** 2);
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
        @generators[1][4] = $lambda1 if $dimension > 1;
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

        my @error-coefficients = 5e0, 5e0, 1e0, 5e0, 0.5e0, 0.25e0;

        return {
            :@weights,
            :@generators,
            :@error-coefficients,
            :@rule-points,
        };
    }

    #| Apply the selected fully symmetric rule to a region.
    method integrate($region --> Math::NIntegrate::Rule::MultiDimensional:D) {
        die 'The rule and region dimensions must match.'
        unless $region.dimension == self.dimension;

        my @weights = |%!data<weights>;
        my @generator-columns = ^@weights[0].elems .map: -> $generator-index {
            %!data<generators>».[$generator-index];
        };
        my $reference-volume = 2e0 ** self.dimension;
        my $rule-scale = 1e0 / $reference-volume;
        my @rule-values = 0e0 xx @weights.elems;

        self.abscissas = [];

        for @generator-columns.kv -> $generator-index, @generator {
            my %seen;
            my $symmetric-sum = 0e0;

            # Generate all sign changes of each permutation.  A generator can
            # contain zero or repeated coordinates, so points are de-duplicated.
            for @generator.permutations -> @permutation {
                for ^$reference-volume.Int -> $sign-mask {
                    my @reference-point = @permutation.kv.map: -> $axis, $value {
                        my $sign = $sign-mask +& (1 +< $axis) ?? -1e0 !! 1e0;
                        $value == 0e0 ?? 0e0 !! $sign * $value;
                    };
                    my $key = @reference-point.join(',');
                    next if %seen{$key}:exists;
                    %seen{$key} = True;

                    # Rules in this package use unit-cube abscissas; the
                    # region maps them to its actual integration bounds.
                    my @point = @reference-point.map({ (1e0 + $_) / 2e0 });
                    self.abscissas.push(@point);
                    my $value = $region.eval-integrand(@point);
                    next unless $value ~~ Numeric:D
                            && !($value.isNaN || $value ~~ Inf | -Inf);
                    $symmetric-sum += $value;
                }
            }

            for ^@weights.elems -> $rule-index {
                @rule-values[$rule-index] += @weights[$rule-index][$generator-index]
                        * $symmetric-sum;
            }
        }

        self.integral = $rule-scale * @rule-values[0];

        # Construct the three normalized combinations of successive null
        # rules, as in DINHRE/DRLHRE, and apply the local error heuristic.
        my @normalized-null-values;
        for ^3 -> $null-index {
            my $largest = 0e0;
            for ^@generator-columns.elems -> $scale-index {
                my $denominator = @weights[$null-index + 1][$scale-index];
                my $scale = $denominator == 0e0
                        ?? 100e0
                        !! -@weights[$null-index + 2][$scale-index] / $denominator;
                my $one-norm = [+] (^@generator-columns.elems).map: -> $generator-index {
                    %!data<rule-points>[$generator-index]
                            * (@weights[$null-index + 2][$generator-index]
                            + $scale * @weights[$null-index + 1][$generator-index]).abs;
                };
                my $candidate = ($rule-scale
                        * (@rule-values[$null-index + 2]
                                + $scale * @rule-values[$null-index + 1])).abs
                        * $reference-volume / $one-norm;
                $largest = max($largest, $candidate);
            }
            @normalized-null-values.push($largest);
        }

        my @error-coefficients = |%!data<error-coefficients>;
        self.error = @error-coefficients[0] * @normalized-null-values[0]
                <= @normalized-null-values[1]
                && @error-coefficients[1] * @normalized-null-values[1]
                <= @normalized-null-values[2]
                ?? @error-coefficients[2] * @normalized-null-values[0]
                !! @error-coefficients[3] * @normalized-null-values.max;
        self.largest-error-axis = 0;

        return self;
    }
}
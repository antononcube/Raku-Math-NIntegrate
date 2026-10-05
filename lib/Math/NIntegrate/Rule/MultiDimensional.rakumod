use v6.d;

use Math::NIntegrate::Rule::General;
use Math::NIntegrate::Utilities;
use Math::NIntegrate::Codes;

class Math::NIntegrate::Rule::MultiDimensional
        is Math::NIntegrate::Rule::General {

    has $.generators;
    has @.rule-points;
    has @.scales;
    has @.norms;

    #======================================================
    # Creators
    #======================================================

    method get-wt-length() {
        do given $!generators {
            when $_ == 9 { self.dimension == 2 ?? 8 !! 9 }
            when $_ == 7 { 6 }
            default {
                # Should we give error;
                warn 'INCORRECT_ARGUMENTS: potential wrong argument.';
                self.rule-points.elems
            }
        }
    }

    submethod TWEAK(:$!generators, :$dimension!, :$working-precision = Num) {
        # The dimension attribute is inherited from Rule::General, so assign it before
        # creating dimension-dependent DCUHRE rule data.
        self.dimension = $dimension;

        my %data = do given $!generators {
            when $_ == 7 {
                self.d07hre(self.dimension) unless self.abscissas.elems;
            }
            when $_ == 9 {
                self.d09hre(self.dimension) unless self.abscissas.elems;
            }
            default {
                die "No multidimensional rule is available with $!generators generators.";
            }
        }

        self.abscissas = |%data<generators>;
        self.weights = |%data<weights>;
        self.error-weights = |%data<error-coefficients>;
        @!rule-points = |%data<rule-points>;

        # Fill-in the scales and norms -- see DEINHR
        #my $wtLength = self.rule-points.elems;
        my $wtLength = self.get-wt-length;

        my @we = 0 xx 14;
        for ^3 -> $k {
            for ^$wtLength -> $i {

                if !is-zero(self.weights[$k + 1][$i]) {
                    @!scales[$k][$i] = - self.weights[$k + 2][$i] / self.weights[$k + 1][$i]
                } else {
                    @!scales[$k][$i] = 100
                }

                for ^$wtLength -> $j {
                    @we[$j] = self.weights[$k + 2][$j] + @!scales[$k][$i] * self.weights[$k + 1][$j]
                }

                @!norms[$k][$i] = 0;

                for ^$wtLength -> $j {
                    @!norms[$k][$i] += @!rule-points[$j] * @we[$j].abs
                }

                @!norms[$k][$i] = 2 ** self.dimension / @!norms[$k][$i]
            }
        }
    }

    multi method new($dimension, :$generators = 9, :prec(:$working-precision) = Num) {
        self.bless(:$generators, :$dimension, :$working-precision)
    }

    multi method new(:$generators = 9, :dim(:$dimension)!, :prec(:$working-precision) = Num) {
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

        my $two-to-dimension = 2 ** $dimension;
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

        # Store offsets whose fully symmetric orbits lie in [-1/2, 1/2]^dim.
        @generators = @generators.map({ $_ <<*>> 0.5 });

        return {
            :@weights,
            :@generators,
            :@error-coefficients,
            :@rule-points,
        };
    }

    #| D09HRE
    method d09hre(UInt:D $ndim, :$working-precision = Num --> Map:D) {
        my $wtleng = $ndim == 2 ?? 8 !! 9;

        my @w      = do for ^5 { [0 xx $wtleng] };
        my @g      = do for ^$ndim { [0 xx $wtleng] };
        my @errcof = 0 xx 6;
        my @rulpts = 0 xx $wtleng;

        my $ratio;
        my $lam0;
        my $lam1;
        my $lam2;
        my $lam3;
        my $lamp;
        my $twondm;

        # Initialize generators, weights and RULPTS
        for 1..$wtleng -> $j {
            for 1..$ndim -> $i {
                @g[$i - 1][$j - 1] = 0;
            }
            for 1..5 -> $i {
                @w[$i - 1][$j - 1] = 0;
            }
            @rulpts[$j - 1] = 2 * $ndim;
        }
        $twondm = 2 ** $ndim;
        @rulpts[$wtleng - 1] = $twondm;
        if $ndim > 2 {
            @rulpts[7] = (4 * $ndim * ($ndim - 1) * ($ndim - 2)) div 3;
        }
        @rulpts[6] = 4 * $ndim * ($ndim - 1);
        @rulpts[5] = 2 * $ndim * ($ndim - 1);
        @rulpts[0] = 1;

        # Compute squared generator parameters
        $lam0 = 0.4707;
        $lam1 = 4 / (15 - 5 / $lam0);
        $ratio = (1 - $lam1 / $lam0) / 27;
        $lam2 = (5 - 7 * $lam1 - 35 * $ratio) / (7 - 35 * $lam1 / 3 - 35 * $ratio / $lam0);
        $ratio = $ratio * (1 - $lam2 / $lam0) / 3;
        $lam3 = (7 - 9 * ($lam2 + $lam1) + 63 * $lam2 * $lam1 / 5 - 63 * $ratio) /
                (9 - 63 * ($lam2 + $lam1) / 5 + 21 * $lam2 * $lam1 - 63 * $ratio / $lam0);
        $lamp = 0.0625;

        # Compute degree 9 rule weights
        @w[0][$wtleng - 1] = 1 / (3 * $lam0) ** 4 / $twondm;
        if $ndim > 2 {
            @w[0][7] = (1 - 1 / (3 * $lam0)) / (6 * $lam1) ** 3;
        }
        @w[0][6] = (1 - 7 * ($lam0 + $lam1) / 5 + 7 * $lam0 * $lam1 / 3) /
                (84 * $lam1 * $lam2 * ($lam2 - $lam0) * ($lam2 - $lam1));
        @w[0][5] = (1 - 7 * ($lam0 + $lam2) / 5 + 7 * $lam0 * $lam2 / 3) /
                (84 * $lam1 * $lam1 * ($lam1 - $lam0) * ($lam1 - $lam2)) -
                @w[0][6] * $lam2 / $lam1 - 2 * ($ndim - 2) * @w[0][7];
        @w[0][3] = (1 - 9 * (($lam0 + $lam1 + $lam2) / 7 - ($lam0 * $lam1 + $lam0 * $lam2 +
                $lam1 * $lam2) / 5) - 3 * $lam0 * $lam1 * $lam2) /
                (18 * $lam3 * ($lam3 - $lam0) * ($lam3 - $lam1) * ($lam3 - $lam2));
        @w[0][2] = (1 - 9 * (($lam0 + $lam1 + $lam3) / 7 - ($lam0 * $lam1 + $lam0 * $lam3 +
                $lam1 * $lam3) / 5) - 3 * $lam0 * $lam1 * $lam3) /
                (18 * $lam2 * ($lam2 - $lam0) * ($lam2 - $lam1) * ($lam2 - $lam3)) -
                2 * ($ndim - 1) * @w[0][6];
        @w[0][1] = (1 - 9 * (($lam0 + $lam2 + $lam3) / 7 - ($lam0 * $lam2 + $lam0 * $lam3 +
                $lam2 * $lam3) / 5) - 3 * $lam0 * $lam2 * $lam3) /
                (18 * $lam1 * ($lam1 - $lam0) * ($lam1 - $lam2) * ($lam1 - $lam3)) -
                2 * ($ndim - 1) * (@w[0][6] + @w[0][5] + ($ndim - 2) * @w[0][7]);

        # Compute weights for 2 degree 7, 1 degree 5 and 1 degree 3 rules
        @w[1][$wtleng - 1] = 1 / (108 * $lam0 ** 4) / $twondm;
        if $ndim > 2 {
            @w[1][7] = (1 - 27 * $twondm * @w[1][8] * $lam0 ** 3) / (6 * $lam1) ** 3;
        }
        @w[1][6] = (1 - 5 * $lam1 / 3 - 15 * $twondm * @w[1][$wtleng - 1] * $lam0 ** 2 * ($lam0 - $lam1)) /
                (60 * $lam1 * $lam2 * ($lam2 - $lam1));
        @w[1][5] = (1 - 9 * (8 * $lam1 * $lam2 * @w[1][6] + $twondm * @w[1][$wtleng - 1] * $lam0 ** 2)) /
                (36 * $lam1 * $lam1) - 2 * @w[1][7] * ($ndim - 2);
        @w[1][3] = (1 - 7 * (($lam1 + $lam2) / 5 - $lam1 * $lam2 / 3 + $twondm * @w[1][$wtleng - 1] *
                $lam0 * ($lam0 - $lam1) * ($lam0 - $lam2))) /
                (14 * $lam3 * ($lam3 - $lam1) * ($lam3 - $lam2));
        @w[1][2] = (1 - 7 * (($lam1 + $lam3) / 5 - $lam1 * $lam3 / 3 + $twondm * @w[1][$wtleng - 1] *
                $lam0 * ($lam0 - $lam1) * ($lam0 - $lam3))) /
                (14 * $lam2 * ($lam2 - $lam1) * ($lam2 - $lam3)) - 2 * ($ndim - 1) * @w[1][6];
        @w[1][1] = (1 - 7 * (($lam2 + $lam3) / 5 - $lam2 * $lam3 / 3 + $twondm * @w[1][$wtleng - 1] *
                $lam0 * ($lam0 - $lam2) * ($lam0 - $lam3))) /
                (14 * $lam1 * ($lam1 - $lam2) * ($lam1 - $lam3)) -
                2 * ($ndim - 1) * (@w[1][6] + @w[1][5] + ($ndim - 2) * @w[1][7]);
        @w[2][$wtleng - 1] = 5 / (324 * $lam0 ** 4) / $twondm;
        if $ndim > 2 {
            @w[2][7] = (1 - 27 * $twondm * @w[2][8] * $lam0 ** 3) / (6 * $lam1) ** 3;
        }
        @w[2][6] = (1 - 5 * $lam1 / 3 - 15 * $twondm * @w[2][$wtleng - 1] * $lam0 ** 2 * ($lam0 - $lam1)) /
                (60 * $lam1 * $lam2 * ($lam2 - $lam1));
        @w[2][5] = (1 - 9 * (8 * $lam1 * $lam2 * @w[2][6] + $twondm * @w[2][$wtleng - 1] * $lam0 ** 2)) /
                (36 * $lam1 * $lam1) - 2 * @w[2][7] * ($ndim - 2);
        @w[2][4] = (1 - 7 * (($lam1 + $lam2) / 5 - $lam1 * $lam2 / 3 + $twondm * @w[2][$wtleng - 1] *
                $lam0 * ($lam0 - $lam1) * ($lam0 - $lam2))) /
                (14 * $lamp * ($lamp - $lam1) * ($lamp - $lam2));
        @w[2][2] = (1 - 7 * (($lam1 + $lamp) / 5 - $lam1 * $lamp / 3 + $twondm * @w[2][$wtleng - 1] *
                $lam0 * ($lam0 - $lam1) * ($lam0 - $lamp))) /
                (14 * $lam2 * ($lam2 - $lam1) * ($lam2 - $lamp)) - 2 * ($ndim - 1) * @w[2][6];
        @w[2][1] = (1 - 7 * (($lam2 + $lamp) / 5 - $lam2 * $lamp / 3 + $twondm * @w[2][$wtleng - 1] *
                $lam0 * ($lam0 - $lam2) * ($lam0 - $lamp))) /
                (14 * $lam1 * ($lam1 - $lam2) * ($lam1 - $lamp)) -
                2 * ($ndim - 1) * (@w[2][6] + @w[2][5] + ($ndim - 2) * @w[2][7]);
        @w[3][$wtleng - 1] = 2 / (81 * $lam0 ** 4) / $twondm;
        if $ndim > 2 {
            @w[3][7] = (2 - 27 * $twondm * @w[3][8] * $lam0 ** 3) / (6 * $lam1) ** 3;
        }
        @w[3][6] = (2 - 15 * $lam1 / 9 - 15 * $twondm * @w[3][$wtleng - 1] * $lam0 * ($lam0 - $lam1)) /
                (60 * $lam1 * $lam2 * ($lam2 - $lam1));
        @w[3][5] = (1 - 9 * (8 * $lam1 * $lam2 * @w[3][6] + $twondm * @w[3][$wtleng - 1] * $lam0 ** 2)) /
                (36 * $lam1 * $lam1) - 2 * @w[3][7] * ($ndim - 2);
        @w[3][3] = (2 - 7 * (($lam1 + $lam2) / 5 - $lam1 * $lam2 / 3 + $twondm * @w[3][$wtleng - 1] *
                $lam0 * ($lam0 - $lam1) * ($lam0 - $lam2))) /
                (14 * $lam3 * ($lam3 - $lam1) * ($lam3 - $lam2));
        @w[3][2] = (2 - 7 * (($lam1 + $lam3) / 5 - $lam1 * $lam3 / 3 + $twondm * @w[3][$wtleng - 1] *
                $lam0 * ($lam0 - $lam1) * ($lam0 - $lam3))) /
                (14 * $lam2 * ($lam2 - $lam1) * ($lam2 - $lam3)) - 2 * ($ndim - 1) * @w[3][6];
        @w[3][1] = (2 - 7 * (($lam2 + $lam3) / 5 - $lam2 * $lam3 / 3 + $twondm * @w[3][$wtleng - 1] *
                $lam0 * ($lam0 - $lam2) * ($lam0 - $lam3))) /
                (14 * $lam1 * ($lam1 - $lam2) * ($lam1 - $lam3)) -
                2 * ($ndim - 1) * (@w[3][6] + @w[3][5] + ($ndim - 2) * @w[3][7]);
        @w[4][1] = 1 / (6 * $lam1);

        # Set generator values
        $lam0 = sqrt($lam0);
        $lam1 = sqrt($lam1);
        $lam2 = sqrt($lam2);
        $lam3 = sqrt($lam3);
        $lamp = sqrt($lamp);
        for 1..$ndim -> $i {
            @g[$i - 1][$wtleng - 1] = $lam0;
        }
        if $ndim > 2 {
            @g[0][7] = $lam1;
            @g[1][7] = $lam1;
            @g[2][7] = $lam1;
        }
        @g[0][6] = $lam1;
        @g[1][6] = $lam2;
        @g[0][5] = $lam1;
        @g[1][5] = $lam1;
        @g[0][4] = $lamp;
        @g[0][3] = $lam3;
        @g[0][2] = $lam2;
        @g[0][1] = $lam1;

        # Compute final weight values.
        # The null rule weights are computed from differences between
        # the degree 9 rule weights and lower degree rule weights.
        @w[0][0] = $twondm;
        for 2..5 -> $j {
            for 2..$wtleng -> $i {
                @w[$j - 1][$i - 1] = @w[$j - 1][$i - 1] - @w[0][$i - 1];
                @w[$j - 1][0] = @w[$j - 1][0] - @rulpts[$i - 1] * @w[$j - 1][$i - 1];
            }
        }
        for 2..$wtleng -> $i {
            @w[0][$i - 1] = $twondm * @w[0][$i - 1];
            @w[0][0] = @w[0][0] - @rulpts[$i - 1] * @w[0][$i - 1];
        }

        # Set error coefficients
        @errcof[0] = 5;
        @errcof[1] = 5;
        @errcof[2] = 1;
        @errcof[3] = 5;
        @errcof[4] = 0.5;
        @errcof[5] = 0.25;

        # Store offsets whose fully symmetric orbits lie in [-1/2, 1/2]^dim.
        @g = @g.map({ $_ <<*>> 0.5 });

        my @weights = @w.map({ $_.map(*.&numerical(:$working-precision)) });
        my @generators = @g.map({ $_.map(*.&numerical(:$working-precision)) });
        my @error-coefficients = @errcof.map({ numerical($_, $working-precision) });
        my @rule-points = @rulpts;

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
    # The memoization did not produce noticeable time changes for the 2D integral
    # nintegrate({ 1 / ($^x + $^y).sqrt }, <x 0 1>, <y 0 1>)
    # has %!permutations;

    #| Generate the distinct points in the fully symmetric orbit of a
    #| non-negative generator.
    sub distinct-permutations(@generators --> Array:D) {
        my @result;

        # Memoization
        # my $key = @generators.join('::');
        #if %!permutations{$key}:exists {
        #    #say "retrieved for $key";
        #    my @res = |%!permutations{$key};
        #    return @res;
        #}

        sub visit(@remaining, @prefix) {
            unless @remaining.elems {
                my @nonzero-axes = ^@prefix.elems
                        .grep({ @prefix[$_] != 0 });

                # A zero coordinate has only one distinct sign.
                for ^(1 +< @nonzero-axes.elems) -> $sign-mask {
                    my @point = @prefix.Array;
                    for @nonzero-axes.kv -> $bit, $axis {
                        @point[$axis] = -@point[$axis]
                        if $sign-mask +& (1 +< $bit);
                    }
                    @result.push(@point);
                }
                return;
            }

            # Selecting each distinct remaining value once generates the
            # multiset permutations without first producing duplicates.
            my @used-values;
            for @remaining.kv -> $index, $value {
                next if @used-values.first(* == $value, :k).defined;
                @used-values.push($value);

                my @rest = @remaining.Array;
                @rest.splice($index, 1);
                visit(@rest, [|@prefix, $value]);
            }
        }

        visit(@generators.Array, []);
        #%!permutations{$key} = @result;
        return @result;
    }

    #| Compute the fully symmetric sum for one generator orbit.
    method defshr(@generators is copy, $region, Numeric:D :$offset = 0) {
        die 'DEFSHR generator dimension does not match the rule dimension.'
        unless @generators.elems == self.dimension;
        die 'DEFSHR generators must be non-negative and non-increasing.'
        unless [&&] @generators.map(* >= 0)
                && [&&] (1 ..^ self.dimension).map({
                    @generators[$_ - 1] >= @generators[$_]
                });

        my Numeric $fulsms = 0e0;
        my Numeric $fulabs = 0e0;
        my Numeric $funvls = 0e0;

        for distinct-permutations(@generators) -> @g {
            my @point = @g.map(* + $offset);
            $funvls = $region.eval-integrand(@point);
            die 'DEFSHR obtained a non-numeric integrand value.'
            unless $funvls ~~ Numeric:D;
            $fulsms += $funvls;
            $fulabs += $funvls.abs;
        }

        return { :$fulsms, :$fulabs, :$funvls };
    }

    #| Compute a DCUHRE basic rule, its error estimate, and its split axis.
    #| The optional roundoff guard suppresses null estimates which contain no
    #| signal above floating-point accumulation error.
    method derlhr(
            $region,
            Bool:D :$roundoff-guard = False
            --> Map:D) {
        die 'The rule and region dimensions must match.'
        unless $region.dimension == self.dimension;

        my UInt:D $wtleng = self.get-wt-length;

        # Correspondences with the original DERLHR signature.
        my @g = self.abscissas;
        my @w = self.weights;
        my @errcof = self.error-weights;
        my @scales = @!scales;
        my @norms = @!norms;

        die 'DERLHR rule data has inconsistent dimensions.'
        unless @g.elems == self.dimension
                && [&&] @g.map(*.elems == $wtleng)
                        && @w.elems >= 5
                        && [&&] @w.map(*.elems == $wtleng)
                                && @errcof.elems >= 4
                                && @scales.elems >= 3
                                && @norms.elems >= 3;

        # Region evaluates points on the unit cube. Stored DCUHRE generators
        # are already offsets in [-1/2, 1/2] from the unit-cube center.
        my Numeric:D $center-coordinate = 0.5e0;
        my Numeric:D $region-volume = 0.5e0 ** self.dimension;
        my @center = $center-coordinate xx self.dimension;
        my @x = @center.Array;

        my UInt:D $division-axis = 0;
        my Numeric $center-value = $region.eval-integrand(@x);
        die 'DERLHR obtained a non-numeric integrand value at the center.'
        unless $center-value ~~ Numeric:D;

        my Numeric $basval = @w[0][0] * $center-value;
        my Numeric $absolute-rule-sum = (@w[0][0] * $center-value).abs;
        my @null = 0e0 xx 8;
        for ^4 -> $null-rule {
            @null[$null-rule] = @w[$null-rule + 1][0] * $center-value;
        }

        my Numeric $difference-maximum = 0e0;
        my @diff = 0e0 xx self.dimension;
        my @order = ^self.dimension;
        my Numeric:D $ratio = (@g[0][2] / @g[0][1]) ** 2;

        # Compute fourth differences and accumulate generator columns 1 and 2.
        for ^self.dimension -> $axis {
            @x = @center.Array;
            @x[$axis] = $center-coordinate - @g[0][1];
            my Numeric $near-minus = $region.eval-integrand(@x);
            @x[$axis] = $center-coordinate + @g[0][1];
            my Numeric $near-plus = $region.eval-integrand(@x);
            @x[$axis] = $center-coordinate - @g[0][2];
            my Numeric $far-minus = $region.eval-integrand(@x);
            @x[$axis] = $center-coordinate + @g[0][2];
            my Numeric $far-plus = $region.eval-integrand(@x);

            die 'DERLHR obtained a non-numeric axial integrand value.'
            unless ($near-minus, $near-plus, $far-minus, $far-plus).all
                    ~~ Numeric:D;

            my Numeric:D $near-sum = $near-minus + $near-plus;
            my Numeric:D $far-sum = $far-minus + $far-plus;
            my Numeric:D $fourth-difference =
                    (2e0 * (1e0 - $ratio) * $center-value
                            - $far-sum + $ratio * $near-sum).abs;
            my Numeric $difference-sum = 0e0;

            # Match DERLHR's roundoff guard.
            $difference-sum += $fourth-difference
            if $center-value.abs + $fourth-difference / 4e0
                    > $center-value.abs;

            for ^4 -> $null-rule {
                @null[$null-rule] += @w[$null-rule + 1][1] * $near-sum
                        + @w[$null-rule + 1][2] * $far-sum;
            }
            $basval += @w[0][1] * $near-sum + @w[0][2] * $far-sum;
            $absolute-rule-sum += @w[0][1].abs
                    * ($near-minus.abs + $near-plus.abs)
                    + @w[0][2].abs
                    * ($far-minus.abs + $far-plus.abs);

            if $difference-sum > $difference-maximum {
                $difference-maximum = $difference-sum;
                $division-axis = $axis;
            }
            @diff[$axis] = $difference-sum;
        }

        # Finish the basic and null rules using fully symmetric sums.
        for 3 ..^ $wtleng -> $generator-index {
            my @generator = @g.map(*.[$generator-index]);
            my %sums = self.defshr(
                    @generator,
                    $region,
                    offset => $center-coordinate);

            for ^4 -> $null-rule {
                @null[$null-rule] +=
                        @w[$null-rule + 1][$generator-index] * %sums<fulsms>;
            }
            $basval += @w[0][$generator-index] * %sums<fulsms>;
            $absolute-rule-sum += @w[0][$generator-index].abs * %sums<fulabs>;
        }

        # Find the greatest normalized estimate in each plane spanned by two
        # successive null rules.
        for ^3 -> $null-rule {
            my Numeric $search = 0e0;
            for ^$wtleng -> $generator-index {
                $search = max(
                        $search,
                        (@null[$null-rule + 1]
                                + @scales[$null-rule][$generator-index]
                                * @null[$null-rule]).abs
                                * @norms[$null-rule][$generator-index]);
            }
            @null[$null-rule] = $search;
        }

        my Numeric $roundoff-threshold = 0e0;
        if $roundoff-guard {
            $roundoff-threshold = 50e0 * $MACHINE_EPSILON * max($absolute-rule-sum, $basval.abs, $center-value.abs);

            for ^3 -> $null-rule {
                @null[$null-rule] = 0e0
                if @null[$null-rule] <= $roundoff-threshold;
            }
        }

        my Numeric $rgnerr =
                @errcof[0] * @null[0] <= @null[1] && @errcof[1] * @null[1] <= @null[2]
                ?? @errcof[2] * @null[0]
                !! @errcof[3] * max(@null[0], @null[1], @null[2]);

        $rgnerr = max($rgnerr, $roundoff-threshold)
        if $roundoff-guard;

        $basval *= $region-volume;
        $rgnerr *= $region-volume;
        my Numeric:D $greate = $rgnerr;

        return {
            :$basval,
            :$rgnerr,
            direct => $division-axis,
            :$greate,
            :@diff,
            :@order,
        };
    }

    #======================================================
    # Integration
    #======================================================

    #| Apply the selected fully symmetric rule to a region.
    method integrate(
            $region,
            Bool:D :$roundoff-guard = False
            --> Math::NIntegrate::Rule::MultiDimensional:D) {
        my %result = self.derlhr($region, :$roundoff-guard);
        self.integral = %result<basval>;
        self.error = %result<rgnerr>;
        self.largest-error-axis = %result<direct>;

        return self;
    }
}
use v6.d;

use Math::NIntegrate::Rule;

# Leaf in the Composite design pattern
class Math::NIntegrate::Rule::General
        is Math::NIntegrate::Rule {
    # Integration rule data
    has @.abscissas is rw;
    has @.weights is rw;
    has @.error-weights is rw;
    has UInt:D $.dimension is rw = 1;

    # Evaluation data
    has Numeric $.integral is rw;
    has Numeric $.error is rw;
    has UInt $.largest-error-axis is rw;

    #======================================================
    # Creators
    #======================================================

    submethod BUILD(:@!abscissas = Empty, :@!weights = Empty, :@!error-weights = Empty) {}

    multi method new(@abscissas, @weights, @error-weights) {
        self.bless(:@abscissas, :@weights, :@error-weights);
    }

    multi method new(:@abscissas, :@weights, :@error-weights) {
        self.bless(:@abscissas, :@weights, :@error-weights);
    }

    #| Copy attributes from an object
    method copy(Math::NIntegrate::Rule::General:D $from, Bool:D :deep(:deep-copy(:$clone)) = False) {
        # @!transforms are functions/subs, so, $from.transforms>>.clone.Array seems inappropriate.
        @!abscissas = $clone ?? $from.abscissas.clone !! $from.abscissas;
        @!weights = $clone ?? $from.weights.clone !! $from.weights;
        @!error-weights = $clone ?? $from.error-weights.clone !! $from.error-weights;

        return self
    }

    #| Clone the object
    method clone(-->Math::NIntegrate::Rule::General:D) {
        Math::NIntegrate::Rule::General.new(Empty, Empty, Empty).copy(self, :clone)
    }

    #======================================================
    # Integration
    #======================================================

    method integrate($region) {

        return self!integrate-en-bloc($region) if $region.is-en-bloc-ready;

        my $integralLocal = 0;
        my $errorLocal = 0;
        my @values;

        # Clear the cache of Jacobian values before integrating.
        # This is done here in order to localize the memoization per rule evaluation.
        # Hence, manage cache size.
        Math::NIntegrate::Utilities::empty-jacobians-cache;

        for ^@!abscissas.elems -> $i {
            my @point = @!abscissas[$i] ~~ Numeric:D ?? @!abscissas[$i] !! |@!abscissas[$i];
            my $value = $region.eval-integrand(@point);

            # Ignoring non-numerical results
            if $value ~~ Numeric:D && !($value.isNaN || $value ~~ Inf | -Inf) {
                # Warnings for NaN and Inf should be given
                @values.push($value);
                $integralLocal += @!weights[$i] * $value;
                $errorLocal += @!error-weights[$i] * $value;
            }
        }

        $!integral = $integralLocal;
        $!error = $errorLocal.abs;
        $!largest-error-axis = 0;

        return self;
    }


    #======================================================
    # Integration En Bloc
    #======================================================

    # Initially I considered to have a separate class, ::Rule::GeneralEnBloc.
    # See the comment next to method ::Rule::MonteCarlo!integrate-en-bloc.

    method !integrate-en-bloc($region) {

        my $integralLocal = 0;
        my $errorLocal = 0;
        my @values;

        # Also has to be done in this method
        Math::NIntegrate::Utilities::empty-jacobians-cache;

        my $dim = $region.dimension;

        # Assign the 1D abscissas as a list -- this probably should be done beforehand
        my @points = @!abscissas.map({ [$_, ]});

        # Transform abscissas
        my %rescaled = |$region.variable-transformer.vtAffine.transform-en-bloc(:@points, jacobian => 1, context => $region.variable-transformer);
        @points = |%rescaled<points>;

        # This has to be done here -- region splitting reverses the range boundaries for regions with end points,
        # hence, the variable transformer factors would reflect that.
        # Since there is no variable transformer below when $region.eval-integrand is called that Jacobian
        # adjustment factor has to be computed here.
        my $factor = [*] |$region.variable-transformer.jacobian-factors;
        %rescaled<jacobian> *= $factor;

        # Prevent variable transformation
        my $vt = $region.variable-transformer;
        $region.variable-transformer = Nil;

        # Integrand evaluation without variable transformation
        @values = @points.kv.map(-> $i, @a {
            my $value = $region.eval-integrand(@a);
            if $value ~~ Numeric:D && !($value.isNaN || $value ~~ Inf | -Inf) {
                # Warnings for NaN and Inf should be given
                $value *= %rescaled<jacobian>;
                @values.push($value);
                $integralLocal += @!weights[$i] * $value;
                $errorLocal += @!error-weights[$i] * $value;
            }
        });

        # Recover the variable transformer
        $region.variable-transformer = $vt;

        $!integral = $integralLocal;
        $!error = $errorLocal.abs;
        $!largest-error-axis = 0;

        return self;
    }

    #======================================================
    # Representation
    #======================================================
    method Str(::?CLASS:D:-->Str) {
        self.gist
    }

    multi method Hash(::?CLASS:D:-->Map:D) {
        my %h =
                name => self.^name,
                :$!dimension,
                :@!abscissas,
                :@!weights,
                :@!error-weights,
                ;
        return %h
    }

    multi method gist(::?CLASS:D:-->Str) {
        self.^name ~ self.Hash.map({ $_.key => $_.value ~~ Positional:D ?? $_.value.elems !! $_.value }).grep(*.key ne 'name').sort(*.key).List.raku;
    }
}
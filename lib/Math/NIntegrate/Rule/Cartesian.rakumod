use v6.d;

use Math::NIntegrate::Rule::General;

# Composite design pattern application over 1D integration rules.
class Math::NIntegrate::Rule::Cartesian
        is Math::NIntegrate::Rule::General {

    # I do not see why the dimension can be Whatever,
    # i.e. determined by the number of components.
    # At moment of this object creation the dimension of the integral should be known.

    has @.components;

    my $msg1Donly = 'Only one dimensional rules are allowed to be Cartesian rule components.';

    #======================================================
    # Creators
    #======================================================
    submethod TWEAK(:@components, :$dimension) {
        self.dimension = $dimension
    }

    method new(:@components, :dim(:$dimension) is copy = Whatever) {

        die 'The argument $dimension is expected to be an integer greater than 1 or Whatever.'
        unless $dimension.isa(Whatever) || $dimension > 1;

        if $dimension.isa(Whatever) {
            die 'If $dimension is Whatever the Cartesian rule creation with less than 2 components is not allowed.'
            unless @components.elems > 1;
            $dimension = @components.elems
        }

        if @components.elems > 0 {
            # Check that all components are 1D integration rules

            die 'Only integration rules can be Cartesia rule components.'
            unless @components.all ~~ Math::NIntegrate::Rule::General;

            die $msg1Donly unless @components>>.dimension.max == 1
        }

        self.bless(:@components, :$dimension)
    }

    #======================================================
    # Components management methods
    #======================================================
    #| Add a rule
    method add(Math::NIntegrate::Rule:D $obj) {
        # Check that the rule is 1D
        die $msg1Donly if $obj.dimension > 1;
        @!components.push($obj);
        return self
    }

    #| Remove a rule object
    multi method remove(Math::NIntegrate::Rule:D $obj) {
        @!components .= grep({ $_ !=== $obj });
        return self
    }

    #| Remove a rule at given index
    multi method remove(Int:D $index) {
        @!components = @!components.splice($index, 1);
        return self
    }

    #======================================================
    # Composing methods
    #======================================================

    #| Copy from object
    method copy(Math::NIntegrate::Rule::Cartesian:D $from,  Bool:D :deep(:deep-copy(:$clone)) = False) {
        # Delegate to parent class
        self.Math::NIntegrate::Rule::General::copy($from, :$clone);

        # The elements of the components array do not need to be cloned
        # because they are just sources of abscissas and weights -- the integration does not happen with them.
        # I.e. their evaluation data attributes (integral, error, last-error-axis) are not used.
        # But to be safer they are cloned if :clone.
        @!components = $clone ?? $from.components !! $from.components>>.clone;

        return self
    }

    method clone(-->Math::NIntegrate::Rule::Cartesian) {
        Math::NIntegrate::Rule::Cartesian.new(dimension => $.dimension).copy(self, :clone)
    }

    #======================================================
    # Populate abscissas and weights
    #======================================================

    #| Make the rule data by Cartesian product of the components abscissas and weights
    method make-rule-data() {

        return if self.abscissas;

        # If the dimension of the composite rule is greater than the number of components
        # the components array is extended to correspond to the dimension.
        # Note that components are not cloned, since they are just sources of data.
        if self.dimension > @!components.elems {
            @!components = flat @!components xx ceiling(self.dimension / @!components.elems);
            @!components .= head(self.dimension)
        }

        self.abscissas = cross(|@!components.map(*.abscissas));
        self.weights = cross(|@!components.map(*.weights)).map({ [*] |$_ });
        self.error-weights = cross(|@!components.map(*.error-weights)).map({ [*] |$_ });

        return self
    }

    #======================================================
    # Integration
    #======================================================

    method integrate($region) {

        self.make-rule-data();

        self.Math::NIntegrate::Rule::General::integrate($region);

        return self;
    }
}

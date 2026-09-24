use v6.d;

# Preprocessing
use Math::NIntegrate::Spec;

# Core classes
use Math::NIntegrate::Region;
use Math::NIntegrate::NumericalFunction;
use Math::NIntegrate::VariableTransformer::Composite;

# Strategies
use Math::NIntegrate::Strategy;
use Math::NIntegrate::Strategy::GlobalAdaptive;

# Rules
use Math::NIntegrate::Rule;
use Math::NIntegrate::Rule::Cartesian;
use Math::NIntegrate::Rule::ClenshawCurtis;
use Math::NIntegrate::Rule::GaussKronrod;
use Math::NIntegrate::Rule::Trapezoidal;

class Math::NIntegrate::Builder {

    has Math::NIntegrate::Strategy $.strategy;

    #------------------------------------------------------
    #| Make variable transformer composite
    method make-variable-transformer-composite(
            :$region
            --> Math::NIntegrate::VariableTransformer::Composite:D) {

        my $obj = Math::NIntegrate::VariableTransformer::Composite.new(:$region);

        # Should the stack have at least one variable transformer?
        # The Infinity transform can be seen as finite-range proxy of the region.
        # The Affine transform is for mapping or the finite boundaries to integration rules [0, 1] abscissas ranges.

        my @min-original-bounds = $region.min;
        my @max-original-bounds = $region.max;
        my $vtInf = Math::NIntegrate::VariableTransformer::Infinity.new(:@min-original-bounds, :@max-original-bounds);
        $obj.add($vtInf);
        return $obj
    }

    #------------------------------------------------------
    #| Make region
    method make-region(
            Math::NIntegrate::NumericalFunction:D $numerical-function,
            @min,
            @max,
            :$rule = Nil,
            :$working-precision = Num
            --> Math::NIntegrate::Region) {

        fail 'DIMENSIONS_DO_NOT_MATCH: region boundaries dimensions do not match.'
        if @min.elems != @max.elems;

        # First a region object is made
        my $region = Math::NIntegrate::Region.new(:$numerical-function, :@min, :@max, :$rule, dimension => @min.elems);

        # Variable transformer (Composite) is made with $region
        my $var-trans = self.make-variable-transformer-composite(:$region);

        # Variable transformer object is set in $region
        $region.variable-transformer = $var-trans;

        return $region
    }

    #------------------------------------------------------
    #| Make numerical function
    method make-numerical-function(
            &function,
            :$working-precision = Num,
            --> Math::NIntegrate::NumericalFunction:D) {

        # Numerical function
        my $nf = Math::NIntegrate::NumericalFunction.new(:&function, $working-precision);

        return $nf
    }

    #-------------------------------------------------------------------------
    #| Make an integration rule object according to method-spec and dimension.
    method make-rule(
            UInt:D :dim(:$dimension)!,
            :$method = Whatever,
            :$working-precision = Num
            --> Math::NIntegrate::Rule:D
                     ) {
        fail 'INCORRECT_ARGUMENTS: the dimension for integration rules is expected to a positive integer.'
        unless $dimension > 0;

        # Not a rigorous check -- just prevent developer mistakes.
        fail 'INCORRECT_ARGUMENTS: the method argument is expected to a Map:D, Whatever, or WhateverCode.'
        unless $method.isa(Whatever) || $method.isa(WhateverCode) || $method ~~ Map:D;

        # A more rigorous rule check.
        with $method {
            fail 'INCORRECT_ARGUMENTS: if $method is Map:D then is expected to have the keys "type", "name", and other options.'
            unless ($method.keys (&) <type name>).elems == 2;

            fail 'INCORRECT_ARGUMENTS: if $method is Map:D its key "type" is expected to have the value "rule".'
            unless $method<type> eq 'rule';
        }

        # Default 1D rule
        my %default-spec = type => 'rule', name => 'ClenshawCurtisRule', points => 5;

        # Create rule by spec
        my $rule = do given $method {
            when $_.isa(Whatever) || $_.isa(WhateverCode) {
                return $dimension == 1
                        ?? self.make-rule(method => %default-spec, :$dimension, :$working-precision)
                        !! self.make-rule(method => {type => 'rule', name => 'CartesianRule', method => %default-spec}, :$dimension, :$working-precision)
            }

            when $_<name> eq 'TrapezoidalRule' {
                Math::NIntegrate::Rule::Trapezoidal.new(points => $_<points> // 10, :$working-precision);
            }

            when $_<name> eq 'ClenshawCurtisRule' {
                Math::NIntegrate::Rule::ClenshawCurtis.new(points => $_<points> // 5, :$working-precision);
            }

            when $_<name> eq 'GaussKronrodRule' {
                die 'The Gauss-Kronrod rule is not implemented yet.'
                # Math::NIntegrate::Rule::GaussKronrod.new(gauss-points => $_<gauss-points> // $_<points> // 5, :$working-precision);
            }

            when $_<name> eq 'LobattoKronrodRule' {
                die 'The Lobatto-Kronrod rule is not implemented yet.'
                # Math::NIntegrate::Rule::LobattoKronrod.new(gauss-points => $_<gauss-points> // $_<points> // 5, :$working-precision);
            }

            when $_<name> eq 'CartesianRule' {
                die 'The Cartesian rule is not implemented yet.'
                # Math::NIntegrate::Rule::Cartesian.new(gauss-points => $_<gauss-points> // $_<points> // 5, :$working-precision);
            }

            default {
                die "Unknown integration rule: $_<name>."
            }
        }

        return $rule
    }

    #------------------------------------------------------------------------------------
    #| Make an integration strategy object according to method-spec and dimension.
    method make-strategy(
            :$method = Whatever,
            :dim(:$dimension) is copy = Whatever,
            :$singularity-depth is copy = Whatever,
            :$max-recursion is copy = Whatever,
            :$min-recursion is copy = Whatever,
            :$max-points is copy = Whatever,
            :$working-precision = Num,
            :$precision-goal = Whatever,
            :$accuracy-goal = Whatever,
            *%args
            --> Math::NIntegrate::Strategy:D
                         ) {

        # In the current design the strategy object is created before the regions.
        # Hence, only dimension is potentially needed for some of the strategies.
        # The most used strategies: GlobalAdaptive, LocalAdaptive, and AdaptiveMonteCarlo
        # do not need to know the dimensions of the integrand and ranges.

        fail 'INCORRECT_ARGUMENTS: the dimension for integration strategies is expected to a positive integer.'
        unless $dimension > 0;

        # Not a rigorous check -- just prevent developer mistakes.
        fail 'INCORRECT_ARGUMENTS: the method argument is expected to a Map:D, Whatever, or WhateverCode.'
        unless $method.isa(Whatever) || $method.isa(WhateverCode) || $method ~~ Map:D;

        # A more rigorous rule check.
        with $method {
            fail 'INCORRECT_ARGUMENTS: if $method is Map:D then is expected to have the keys "type", "name", and other options.'
            unless ($method.keys (&) <type name>).elems == 2;

            fail 'INCORRECT_ARGUMENTS: if $method is Map:D its key "type" is expected to have the value "strategy".'
            unless $method<type> eq 'strategy';
        }

        # Reassign options
        $max-points = $method<max-points> // $max-points // Whatever;
        $max-recursion = $method<max-recursion> // $max-recursion // 12;
        $min-recursion = $method<min-recursion> // $min-recursion // 0;
        $singularity-depth = $method<singularity-depth> // $singularity-depth // 4;

        # Default strategy
        #my %default-spec = type => 'strategy', name => 'GlobalAdaptive', min-recursion => 0, max-recursion => 12, singularity-depth => 4, max-points => Whatever;
        my %default-spec = type => 'strategy', name => 'GlobalAdaptive', method => Whatever;

        # Create strategy by spec
        my $strategy = do given $method {
            when $_.isa(Whatever) || $_.isa(WhateverCode) {
                return self.make-strategy(
                        method => %default-spec,
                        :$dimension,
                        :$singularity-depth,
                        :$max-points,
                        :$max-recursion,
                        :$min-recursion,
                        :$precision-goal,
                        :$accuracy-goal,
                        |%args
                        )
            }

            when $_<name> eq 'GlobalAdaptive' {
                Math::NIntegrate::Strategy::GlobalAdaptive.new(
                        :$min-recursion,
                        :$max-recursion,
                        :$singularity-depth,
                        :$max-points
                        );
            }

            when $_<name> eq 'LocalAdaptive' {
                die 'The LocalAdaptive strategy is not implemented yet.'
#                Math::NIntegrate::Strategy::LocalAdaptive.new(
#                        :$min-recursion,
#                        :$max-recursion,
#                        :$singularity-depth,
#                        :$max-points
#                        );
            }

            when $_<name> eq 'AdaptiveMonteCarlo' {
                die 'The AdaptiveMonteCarlo strategy is not implemented yet.'
                #                Math::NIntegrate::Strategy::LocalAdaptive.new(
                #                        :$min-recursion,
                #                        :$max-recursion,
                #                        :$singularity-depth,
                #                        :$max-points
                #                        );
            }

            default {
                die "Unknown integration strategy: $_<name>."
            }
        }

        return $strategy
    }

    #------------------------------------------------------
    #| Full integrator creation
    method make-integrator(
            Math::NIntegrate::NumericalFunction :$integrand!,
            :%ranges!,
            :$method = Whatever,
            :$singularity-depth = Whatever,
            :$max-recursion is copy = Whatever,
            :$min-recursion is copy = Whatever,
            :$max-points is copy = Whatever,
            :$working-precision = Num,
            :$precision-goal = Whatever,
            :$accuracy-goal = Whatever,
            *%args
            --> Math::NIntegrate::Strategy:D
                         ) {

        # The method option is parsed and validated at this point.

        # Default strategy options
        my %default = singularity-depth => 4, min-recursion => 4, max-recursion => 12;

        # Bounds
        # Ranges are already normalized at this point.
        my %bounds = min => [], max => [];
        %ranges.sort(*.value<index>).map({
            %bounds<min>.push($_.value<min>);
            %bounds<max>.push($_.value<max>)
        });

        # Integration rule
        my $rule = self.make-rule(dimension => %bounds<min>.elems, :$working-precision);

        # Make regions
        # The numerical function object for the integrand is already made
        my @regions = self.make-region($integrand, %bounds<min>, %bounds<max>, :$rule);

        # In the future:
        # - More than one region is obtained from the original ranges
        # - The integration rule object is set to all regions

        # Make the strategy
        $!strategy = self.make-strategy(dim => %ranges.elems, |%args);

        # Attach regions to strategies
        $!strategy.regions = |@regions;

        return $!strategy
    }
}

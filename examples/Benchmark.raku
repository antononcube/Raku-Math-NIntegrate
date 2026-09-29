#!/usr/bin/env raku
use v6.d;

# Created / shared by Andrew Shitov.
# Modified to take multiple methods by Anton Antonov.

use Math::NIntegrate;
use Math::NIntegrate::Builder;
use Math::NIntegrate::NumericalFunction;
use Math::NIntegrate::Rule::NewtonCotes;
use Math::NIntegrate::Rule::MonteCarlo;

# One integration method per run: MAIN(method, reps).
# Prints the results (for cross-engine comparison) and, on the last line,
# "TIME <ms>" for the timed workload only (module load excluded).

my @one-d =
        ('x²',        -> $x { $x ** 2 },       ['x', 0, 1]),
        ('1/√x',      -> $x { 1 / sqrt($x) },  ['x', 0, 2]),
        ('exp(-x²)',  -> $x { exp(-$x ** 2) }, ['x', 0, Inf]),
        ('sin²x',     -> $x { sin($x) ** 2 },  ['x', 0, 10]);

my @two-d =
        ('x+y²',      -> $x, $y { $x + $y ** 2 },       ['x', 0, 3], ['y', 0, 5]),
        ('1/√(x+y)',  -> $x, $y { 1 / sqrt($x + $y) },  ['x', 0, 1], ['y', 0, 1]);

my @three-d =
        ('exp(-r²)',  -> $x, $y, $z { exp(-($x*$x + $y*$y + $z*$z)) }, ['x', 0, 1], ['y', 0, 1], ['z', 0, 1]),;

sub via-builder($rule, &f, @min, @max, :$wp = Num) {
    my $b = Math::NIntegrate::Builder.new;
    my $nf = Math::NIntegrate::NumericalFunction.new(&f, $wp);
    my $region = $b.make-region($nf, @min, @max, :$rule);
    my $s = $b.make-strategy(dimension => @min.elems, min-recursion => 0, max-recursion => 12);
    $s.regions = $region;
    $s.algorithm(relative-tolerance => 1e-6, absolute-tolerance => 0)
}

my %work =
        'gauss-kronrod' => { @one-d.map({ .[0] => nintegrate(.[1], .[2], method => 'gauss-kronrod-rule'):pairs }) },
        'clenshaw-curtis' => { @one-d.map({ .[0] => nintegrate(.[1], .[2], method => 'clenshaw-curtis-rule'):pairs }) },
        'trapezoidal' => { @one-d.map({ .[0] => nintegrate(.[1], .[2], method => 'trapezoidal-rule'):pairs }) },
        'newton-cotes' => { @one-d.map({ .[0] => via-builder(Math::NIntegrate::Rule::NewtonCotes.new(5), .[1], [.[2][1],], [.[2][2],]) }) },
        'monte-carlo' => {
            # Random: a fixed amount of work (one rule application per region, 20 regions).
            (^20).map(-> $i {
                my $b = Math::NIntegrate::Builder.new;
                my $nf = Math::NIntegrate::NumericalFunction.new(-> $x { $x ** 2 });
                my $region = $b.make-region($nf, [0,], [2,], rule => Math::NIntegrate::Rule::MonteCarlo.new(1000));
                $region.apply-rule;
                "r$i" => %(integral => $region.integral, error => $region.error, region-count => 1)
            })
        },
        'multidimensional' => { (|@two-d, |@three-d).map({ .[0] => nintegrate(.[1], |.[2..*], method => 'multidimensional-rule'):pairs }) },
        'cartesian-gk' => { @two-d.map({ .[0] => nintegrate(.[1], |.[2..*], method => ('global-adaptive', method => ('gauss-kronrod-rule', points => 5))):pairs }) },
        'gauss-kronrod-rat' => { @one-d[0,1].map({ .[0] => nintegrate(.[1], .[2], method => 'gauss-kronrod-rule', working-precision => Rat, precision-goal => 8):pairs }) },
        'clenshaw-curtis-rat' => { @one-d[0,1].map({ .[0] => nintegrate(.[1], .[2], method => 'clenshaw-curtis-rule', working-precision => Rat, precision-goal => 8):pairs }) },
        'trapezoidal-rat' => { @one-d[0,1].map({ .[0] => nintegrate(.[1], .[2], method => 'trapezoidal-rule', working-precision => Rat, precision-goal => 8):pairs }) };

sub benchmark-run(Str:D $method, Int $reps = 1) {
    die "unknown method $method" unless %work{$method}:exists;
    my @res;
    my $t0 = now;
    @res = %work{$method}().eager for ^$reps;
    my $ms = (now - $t0) * 1000;
    for @res -> $p {
        my %r = $p.value;
        say sprintf('%-10s integral=%.10g error=%.3g regions=%s',
                $p.key, %r<integral>.Num, %r<error>.Num, %r<region-count>);
    }
    say sprintf('TIME %.1f', $ms);
    return {results => @res, time => $ms}
}

sub MAIN(Str:D $method = 'clenshaw-curtis, gauss-kronrod, newton-cotes, trapezoidal', Int $reps = 1) {
    my @methods = $method.split(/\s* ',' \s* | \s+/, :skip-empty)>>.trim;
    @methods.map({ benchmark-run($_, $reps) })
}

use v6.d;

unit module Math::NIntegrate::Rescale;

sub rescale-finite($x, $min, $max, $vmin, $vmax) {
    return (-($min * $vmax) + $max * $vmin + ($vmax - $vmin) * $x) / ($max - $min);
}

sub rescale-inf($x, $a, $b, $A, $B) {
    my $a2 = -Inf < $a < Inf ?? 1 !! $a;
    my $b2 = -Inf < $b < Inf ?? 1 !! $b;
    my $A2 = -Inf < $A < Inf ?? 1 !! $A;
    my $B2 = -Inf < $B < Inf ?? 1 !! $B;
    given ($a2, $b2, $A2, $B2) {
        when (1, Inf, 1, Inf) {
            return $b + ($x - $a) * ($B / abs($B)) / ($A / abs($A));
        }
        when (1, Inf, 1, 1) {
            return $B - (-$A + $B)/(1 + $x - $a);
        }
        when (-Inf, 1, 1, 1) {
            return $A - ($A - $B)/(1 - $x + $b);
        }
        when [1, Inf, Inf, 1] {
            return $b + $A / (($x - $a) * abs($A)) * ($B / abs($B));
        }
        when [Inf, 1, Inf, 1] {
            return rescale-inf($x, $a, Inf, $b, Inf);
        }
        when [Inf, 1, 1, Inf] {
            return rescale-inf($x, $a, Inf, Inf, $b);
        }
        when [1, Inf, -Inf, Inf] {
            return rescale-inf(rescale-inf($x, $a, Inf, 1, 0), 0, 1, -Inf, Inf);
        }
        when [Inf, 1, -Inf, Inf] {
            return rescale-inf(rescale-inf($x, Inf, $a, 0, 1), 0, 1, -Inf, Inf);
        }
        when [1, Inf, Inf, -Inf] {
            return rescale-inf(rescale-inf($x, $a, Inf, 1, 0), 0, 1, Inf, -Inf);
        }
        when [Inf, 1, Inf, -Inf] {
            return rescale-inf(rescale-inf($x, Inf, $a, 0, 1), 0, 1, Inf, -Inf);
        }
        when [-Inf, Inf, 1, Inf] {
            return rescale-inf(rescale-inf($x, -Inf, Inf, 0, 1), 0, 1, $b, Inf);
        }
        when [-Inf, Inf, Inf, 1] {
            return rescale-inf(rescale-inf($x, -Inf, Inf, 0, 1), 0, 1, Inf, $b);
        }
        when (Inf, -Inf, 1, 1) {
            return -rescale-inf($x, -Inf, Inf, $A, $B) if $A === Inf || $B === Inf;
        }
        when [-Inf, Inf, -Inf, Inf] {
            return $x;
        }
        when [-Inf, Inf, Inf, -Inf] {
            return -$x;
        }
        when [Inf, -Inf, -Inf, Inf] {
            return -$x;
        }
        when [Inf, -Inf, Inf, -Inf] {
            return $x;
        }
        default {
            note 'non-inf';
            rescale-finite($x, $a, $b, $A, $B)
        }
    }
}

#----------------------------------------------------------
proto sub rescale(Numeric:D $x, |) is export {*}

multi sub rescale(Numeric:D $x, $domain is copy = Whatever, $codomain is copy = Whatever) {
    my @res = rescale([$x,], $domain, $codomain);
    return @res.head;
}

multi sub rescale(@x, $domain is copy = Whatever, $codomain is copy = Whatever) {
    if $domain.isa(Whatever) { $domain = (@x.min, @x.max) }
    die "The second argument is expected to be a list of two numbers or Whatever."
    unless $domain ~~ (List:D | Array:D | Seq:D) && $domain.elems == 2 && $domain.all ~~ Numeric:D;

    if $codomain.isa(Whatever) { $codomain = (0, 1) }
    die "The third argument is expected to be a list of two numbers or Whatever."
    unless $codomain ~~ (List:D | Array:D | Seq:D) && $codomain.elems == 2 && $codomain.all ~~ Numeric:D;


    my ($min, $max) = |$domain;
    my ($vmin, $vmax) = |$codomain;

    if -Inf < $min < Inf && -Inf < $max < Inf && -Inf < $vmin < Inf && -Inf < $vmax < Inf {
        return @x.map({ rescale-finite($_, $min, $max, $vmin, $vmax) })
    } else {
        return @x.map({ rescale-inf($_, $min, $max, $vmin, $vmax) })
    }
}
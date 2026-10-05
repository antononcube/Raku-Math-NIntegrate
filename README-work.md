# Math::NIntegrate

[![MacOS](https://github.com/antononcube/Raku-Math-NIntegrate/actions/workflows/macos.yml/badge.svg)](https://github.com/antononcube/Raku-Math-NIntegrate/actions/workflows/macos.yml)
[![Linux](https://github.com/antononcube/Raku-Math-NIntegrate/actions/workflows/linux.yml/badge.svg)](https://github.com/antononcube/Raku-Math-NIntegrate/actions/workflows/linux.yml)
[![Win64](https://github.com/antononcube/Raku-Math-NIntegrate/actions/workflows/windows.yml/badge.svg)](https://github.com/antononcube/Raku-Math-NIntegrate/actions/workflows/windows.yml)
[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic%202.0-0298c3.svg)](https://opensource.org/licenses/Artistic-2.0)
[![https://raku.land/zef:antononcube/Math::NIntegrate](https://raku.land/zef:antononcube/Math::NIntegrate/badges/version)](https://raku.land/zef:antononcube/Math::NIntegrate)
[![https://raku.land/zef:antononcube/Math::NIntegrate](https://raku.land/zef:antononcube/Math::NIntegrate/badges/downloads)](https://raku.land/zef:antononcube/Math::NIntegrate)

## Introduction

This repository has the Raku code of a library for numerical integration.
The design and architecture of library's main function, `NIntegrate`, resembles that of 
[Wolfram Language `NIntegrate`](https://reference.wolfram.com/language/ref/NIntegrate.html),
[WRI1, WRI2].


----

## Installation

From Zef ecosystem:

```
zef install Math::NIntegrate
```

From GitHub:

```
zef install https://github.com/antononcube/Raku-Math-NIntegrate.git
```

----

## Motivation

While working for 
[Wolfram Research Inc.](https://en.wikipedia.org/wiki/Wolfram_Research),
during the years 2003-2006 I designed and implemented Mathematica's "new" `NIntegrate` and extensively documented it. 
(I wrote [WRI2].)

The Raku programming language has some built-in features that give the ability to:

- Program in an object-oriented manner

- Do [bignum arithmetic](https://en.wikipedia.org/wiki/Arbitrary-precision_arithmetic) 
  computations

- Do symbolic computations

These are also programming language features on which Wolfram Language's `NIntegrate` is based upon.

Hence, it seems natural to think that an implementation of a powerful numerical integration
framework in Raku that has unique features is achievable, worthy, and rewarding.

I plan to write a book that describes in detail software architecture designs and decisions
for the development of numerical integration frameworks. That book is going to use this library.

------

## Design (in brief)

The `NIntegrate` framework is based on Object-Oriented Programming (OOP)
[Design Patterns (DP)](https://en.wikipedia.org/wiki/Design_Patterns).
`NIntegrate` uses the software design patterns Strategy, Composite, Decorator, and others.
 
Here is a description using the OOP DP language: 

> - The basic objects of `NIntegrate` are integration regions. 
> - Each region has its own integration function and integration rule.
> - Conceptually, there are two main types of algorithms: integration strategies and integration rules.
> - The integration strategies use the 
> [Template method](https://en.wikipedia.org/wiki/Template_method_pattern) 
> for their "logic." 
> - The integration regions use 
> [Strategy](https://en.wikipedia.org/wiki/Strategy_pattern) 
> for the computation of integral and error estimates.
> - The integration regions can utilize singularity handler objects. 
> - Creations of integration rules generally use 
> [Builder](https://en.wikipedia.org/wiki/Builder_pattern).
> - Symbolic preprocessing is done through 
> [Decorator](https://en.wikipedia.org/wiki/Decorator_pattern).
> - User specifications are translated into creations of numerical integration algorithm objects 
> through
> [Interpreter](https://en.wikipedia.org/wiki/Interpreter_pattern).
> - The creation of the "final" integration algorithm object uses
> [Abstract factory](https://en.wikipedia.org/wiki/Abstract_factory_pattern).

------

## Usage examples

Basic usage examples:

```raku
use Math::NIntegrate;

nintegrate( -> $x { 1/sqrt($x) }, ['x', 0, 2] );
```

With the adverb "pairs" the result shows the integral and error estimates together with the number of subregions used:

```raku
nintegrate( -> $x { 1 / $x ** 2 }, <x 1 2>, working-precision => Rat ):pairs;
```

Compute a two-dimensional integral:

```raku
nintegrate( -> $x, $y { $x + $y² }, <x 0 2>, <y 0 12>);
```

Adaptive Monte-Carlo integration (not implemented yet):

```raku, eval=FALSE
nintegrate( -> $x, $y, $z { $x + $y² + 1/$z⁻¹ }, x => (0, 2), y => (0, 12), z => (1, 4), method => 'adaptive-monte-carlo' );
```

Integration with functional boundaries 
(using a 2D Cartesian integration rule specified with a 1D component):

```raku
nintegrate( -> *@_ { 1 }, ['x', 0, 1], ['y', 0, {$^x.sqrt}], method => 'clenshaw-curtis-rule'):pairs
```

Compute a two-dimensional integral with using a full specification of a Cartesian integration rule 
which based on a one-dimensional [Gauss-Kronrod rule](https://en.wikipedia.org/wiki/Gauss–Kronrod_quadrature_formula) with 5 Gauss points 
and a one-dimensional [Trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule) with 10 base-points:

```raku
nintegrate( { ($^x + $^y ** 2).sqrt }, <x 0 1>, <y 0 1>, 
        method => ('global-adaptive', 
                   method => ('CartesianRule',
                              method => (('gauss-kronrod-rule', points => 5), ('trapezoidal-rule', points => 10))
                   )
        )):pairs
```

Compare with the Wolfram Language results:

```shell
wolframscript -code 'Through[{Integrate, N@*Integrate, NIntegrate}[Sqrt[x+y^2], {x, 0, 1}, {y, 0, 1}]]'
```

Utilization through a DSL specification (not implemented yet):

```
integrate 1/x^2 over the range [1,2] with a local adaptive strategy and precision goal 5
```

------

## CLI

The package provides a Command Line Interface (CLI) script. Here is its usage message:

```shell
nintegrate --help
```

Here is an example invocation:

```shell
nintegrate '{$^x + $^y}' '<x 0 1>' '<y 0 10>' --pairs
```

------

## Status

The implementation status is tracked in the [Org-mode](https://orgmode.org) file [Work-plan.org](./org/Work-plan.org).

------

## References

[WRI1]
Wolfram Research (1988), 
[`NIntegrate`](https://reference.wolfram.com/language/ref/NIntegrate.html), 
Wolfram Language function, 
https://reference.wolfram.com/language/ref/NIntegrate.html (updated 2014).

[WRI2]
Wolfram Research (2006), 
[Advanced Numerical Integration in the Wolfram Language](https://reference.wolfram.com/language/tutorial/NIntegrateOverview.html),
Wolfram Monograph,
https://reference.wolfram.com/language/tutorial/NIntegrateOverview.html.

[AAr1] 
Anton Antonov,
[`NIntegrate` - The Missing Manual](https://github.com/antononcube/NIntegrateTheMissingManual-book),
(2019),
[GitHub/antononcube](https://github.com/antononcube).

------
Anton Antonov   
Windermere, Florida, USA  
2021-04-05   
2026-09-28
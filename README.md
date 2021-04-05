# Math::NIntegrate

## Introduction

This repository has the Raku code of a library for numerical integration.
The design and architecture of library's main function, `NIntegrate`, resembles that of 
[Mathematica's `NIntegrate`](https://reference.wolfram.com/language/ref/NIntegrate.html),
[WRI1, WRI2].

While working for 
[Wolfram Research Inc.](https://en.wikipedia.org/wiki/Wolfram_Research),
during the years 2003-2006 I designed and implemented Mathematica's "new" `NIntegrate` and extensively documented it. 
(I wrote [WRI2].)

The Raku programming language has some built-in features that give the ability to:

- Program in an object-oriented manner

- Do [bignum arithmetic](https://en.wikipedia.org/wiki/Arbitrary-precision_arithmetic) 
  computations

- Do symbolic computations

These are also programming language features on which Mathematica's `NIntegrate` is based upon.

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

## Usage example

Basic usage examples:

```raku
use Math::NIntegrate;

say NIntegrate( -> $x { 1/sqrt($x) }, x => (0, 2) )

say NIntegrate( -> $x { 1/$x^2 }, x => (1, 2), Method => 'LocalAdaptive', WorkingPrecision => 40 );

say NIntegrate( -> $x, $y { $x + $y^2 }, x => (0, 2), y => (0, 12), Method => ('GlobalAdaptive' Method => ('GaussKronrod', Points => 5) );

say NIntegrate( -> $x, $y, $z { $x + $y^2 + 1/$z }, x => (0, 2), y => (0, 12), z => (1, 4), Method => 'AdaptiveMonteCarlo' )
```

Utilization through a DSL specification:

```
integrate 1/x^2 over the range [1,2] with a local adaptive strategy
```

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
# Performance of integration algorithms

## Main challenge

The first versions "Math::NIntegrate" Raku's `nintegrate` is 100 times slower than Wolfram Language's `NIntegrate`.

The main challenge is that with the original Object-Oriented Design (OOP) of "Math::NIntegrate" 
too many variable transformations have to be done per integrand evaluation. 
Let us call that "too vertical evaluation". 

Hence, an obvious way to speed-up computations is to make those computations as "horizontal" as possible. 
For example, affine transformations can be en bloc, for the whole set of abscissas, instead of one by one.

Another way to compute "horizontally" is to push the variable transformations from Raku to C via NativeCall.
See ["Fast-Affine.raku](../../experiments/Fast-Affine.raku) -- by using native sparse matrix computations 
the affine transformations can be sped-up ≈10 times.

A third way is to use memoization for the results of routines that are frequently called with the same arguments.
See the memoization of `Math::NIntegrate::Utilities::length-calc-md` and corresponding "management call" in
`Math::NIntegrate::Rule::General.integrate`. 
Another memoization is for the distinct permutations for the fully symmetric rules.
See `Math::NIntegrate::Rule::MultiDimensional!distinct-permutations`.

---

## Profiling

Profiling with Rakudo is very cumbersome, so, I did profiling with [Raku++ (aka Rakupp)](https://github.com/ash/rakupp).

I profiled first with `::Rule::MonteCarlo`, because it is natural to use many (random) points 
per Monte Carlo rule evaluation and, therefore, the Monte Carlo rule integration has to be fast.
See ["Profiling-MonteCarlo.raku](../../examples/Profiling-MonteCarlo.raku).


---

## En bloc affine transformations (2026-10-09 ÷ 2026-10-10)

In version 0.0.7 en bloc affine transformations were implemented for constant, finite ranges -- see the method 
`Math::NIntegrate::VariableTransformation::Affine.transform-en-bloc`.

There is a corresponding en bloc integration method `integrate-en-bloc` in the classes
`Math::NIntegrate::Rule::MonteCarlo` and `Math::NIntegrate::Rule::General`.

Initially `Math::NIntegrate::Rule::General::Cartsian` had its own method `integrate-en-bloc`, but after
adding the attributes `@!last-abscissas` and `@!last-values` to `::Rule::General` the integration code of
`::Rule::Cartesian` was greatly simplified.

Those implementations gave:
- ≈ 5 times speed up on multidimensional crude Monte Carlo integration
- ≈ 4.5 times speed up for Cartesian rules
- ≈ 2 times speed up for one-dimensional integrals

In principle, for infinite ranges the same transformation pipeline can be applied 
for `::VariableTransform::Infinite`, but requires some additional verifications and restrictions.

Having the en bloc optimization for "just" constant, finite ranges, though, is probably good enough, 
since the symbolic strategy decorator `Math::NIntegrate::Strategy::UnitCubeRescaling` can be put both infinite and 
functional ranges into constant, finite ones.

---


## Additional notes

### Chat with Andrew Shitov on 2026-10-07

**Remark:** Andrew Shitov is the author of [Rakupp](https://github.com/ash/rakupp).

- AA: Hi! Rakupp is twice faster than Rakudo on “Math::NIntegrate” benchmarks with multidimensional integrals.
- AA: Also, memoization performance improvements I made based on Rakupp’s profiling have effect in Rakupp. But not in Rakudo.
- AA: “Math::NIntegrate” is still slow, though. One way to speed up its computation is the to use native sparse matrices.
- AA: (Hence, my “Math::SparseMatrix::Native” issue proclaiming, yesterday.)
- AA: BTW, I think “Sub::Memoized” did not work on Rakupp.
- AA: So, I did an ad hoc memoization. (Which I was inclined to do, anyway.)
- AS: << confirms failure to install "Sub::Memoized" >> 
- AA: That package adds a new trait. Cool and instructive. I am not sure is it of high priority.
- AA: I am getting curious to know why Rakudo is slower on my numerical integration benchmarks. I might try to use its profile today. (Which is way too cumbersome.)
- AS: << conjectures that there is no easy, single answer >>
- AA: Yeah, sure. But in some sense I am looking also for inspiration for what to optimize.
- AA: The problem with my current design is that there are too many nested variable transformations which are done per single integration point. Instead of doing those en bloc.
- AA: For certain simple cases — constant integration ranges without infinities — there is an obvious strategy how to speed up the computations.
- AA: At least in Rakupp that was obvious and gave results. :)
- AA: Another thing I noticed in Rakupp — the evaluation of the integrand function is / was faster than the variable transformation.
- AA: Integrands look like `{$^x + $^y ** 2}`.
- AA: Long way to go speed wise. For numerical integration WL is 100 and 50 times faster than Rakudo and Rakupp, respectively.
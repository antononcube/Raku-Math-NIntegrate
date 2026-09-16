# Methodology and Design 
## ***... for the Raku package "Math::NIntegrate"***

---

## Introduction

---

## Overall process

Each integration strategy of "Math::NIntegrate" creates and manipulates a collection of integration regions. 
Each integration region can have its own integrand and/or integration rule. 
"Math::NIntegrate"'s main integration strategy "GlobalAdaptive" keeps the integration regions in a heap according to their error. 
The sum of the integral estimates of all regions make the global integral estimate. 
The sum of the integral errors make the global error. 
If the global error is larger than the desired tolerance "GlobalAdaptive" splits the region with the largest error estimate into two regions and applies the integration rule. 
If too many splittings have been done then a singularity handler is applied over the last region split.

At each step of an integration strategy the option "IntegrationMonitor" obtains as an argument the list of integration regions used in that step. 
Below is a table that shows methods that can be applied to each integration region in that list.

| Property               | Description                                                                                                          |
|------------------------|----------------------------------------------------------------------------------------------------------------------|
| `"Axis"`               | Returns the axis with maximum error.                                                                                 |
| `"Boundaries"`         | Returns the boundaries of the region.                                                                                |
| `"Dimension"`          | Returns the dimension of the region.                                                                                 |
| `"Error"`              | Returns the error estimate of the last integral estimation.                                                          |
| `"GetRule"`            | Returns the attached integration rule.                                                                               |
| `"GetSamplingPoints"`  | Returns the sampling points used for the integration estimation.                                                     |
| `"GetValues"`          | Returns the integrand values over the sampling points used for the integration estimation.                           |
| `"Integral"`           | Returns the integral estimate of the last integral estimation.                                                       |
| `"Integrand"`          | Returns the integrand (numerical function) attached to the region.                                                   |
| `"NumericalFunction"`  | Returns the integrand (numerical function) attached to the region.                                                   |
| `"OriginalBoundaries"` | Returns the original boundaries of the region (when singularity handling is applied it differs from `"Boundaries"`). |
| `"Properties"`         | Returns all methods applicable.                                                                                      |
| `"WorkingPrecision"`   | Returns the precision with which the function `N` is used during evaluation.                                         |

--- 

## Object-oriented design

---

## References

### Article, blog posts

[AAmse1] Anton Antonov,
["Determining which rule NIntegrate selects automatically", answer](https://mathematica.stackexchange.com/a/96663),
(2015)
[MathematicaStackExchange](https://mathematica.stackexchange.com).

### Wolfram Language documentation

### Raku packages

[AAp1] Anton Antonov,
[Math::NIntegrate, Raku package](https://github.com/antononcube/Raku-Data-NIntegrate),
(2021-2026),
[GitHub/antononcube](https://github.com/antononcube).

[AAp1] Anton Antonov, 
[Data::Transformers, Raku package](https://github.com/antononcube/Raku-Data-Transformers),
(2026),
[GitHub/antononcube](https://github.com/antononcube).

[AAp2] Anton Antonov,
[LeftistHeap, Raku package](https://github.com/antononcube/Raku-LeftistHeap),
(2026),
[GitHub/antononcube](https://github.com/antononcube).

### Repositories

[AAr1] Anton Antonov,
["NIntegrate, the missing manual"](https://github.com/antononcube/NIntegrateTheMissingManual-book),
(2019-2026),
[GitHub/antononcube](https://github.com/antononcube).
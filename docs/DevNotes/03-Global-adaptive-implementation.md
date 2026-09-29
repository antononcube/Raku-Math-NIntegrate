## GlobalAdaptive algorithm implementation

- This the primary algorithm of the "Math::NIntegrate" framework.
- I implemented "LeftistHeap" a few weeks before my serious work "Math::NIntegrate" started.
  - For the purposes of "Graph" use in/with Rakupp, but also because of the need to have a heap in `NIntegrate`.
- Pretty straightforward implementation if region splitting and integration are implemented.
- Initially, the implementation was without the min-recursion partitioning and the without singularity handler application. 
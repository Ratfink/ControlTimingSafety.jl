# Probable Safety

In the paper "Statistical Hypothesis Testing of Controller Implementations
Under Timing Uncertainties," we developed a method for estimating the maximum
deviation of a control system that may experience deadline misses.  This
estimate is based on Jeffreys' Bayes factor hypothesis testing, and is
probabilistically guaranteed.  The method is implemented in the
[`estimate_deviation`](@ref) function.

```@docs
estimate_deviation
maximum_deviation_random
```

# ControlTimingSafety.jl

*Tools for verifying safety of control systems that may experience timing uncertainty*

This package provides tools for calculating safe overapproximations and
probabilistic estimations of the maximum deviation a control system may suffer
due to deadline misses.  If you find these tools useful in your research,
please consider [Citing](@ref) the papers that the algorithms originated in.

## Installation

From the Julia REPL, run:

```julia
using Pkg; Pkg.add(url="https://github.com/Ratfink/ControlTimingSafety.jl")
```

## Usage

Given a discrete-time state space model `sys` from the
[ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl) package,
and a feedback gain matrix `K`, first create an [`Automaton`](@ref) object using one of
the constructors for the strategies described in the paper:

```julia
hk = hold_kill(sys, K, max_misses)
hs = hold_skip_next(sys, K, max_misses)
zk = zero_kill(sys, K, max_misses)
zs = zero_skip_next(sys, K, max_misses)
```

To run the Bounded Runs Iteration algorithm for this automaton, first create
the initial set, then run the algorithm for a given per-iteration run
length `n` to produce the reachable set of length `H`. The number of iterations
is the smallest integer `t` such that `n`×`t`≥`H`.

```julia
augbounds = augment(automaton, bounds)
all_bounds = bounded_runs_iter(automaton, augbounds, n, H)
```

The deviation at each time step can then be computed thusly:

```julia
d = deviation(automaton, augbounds, all_bounds)
```

For a complete example, see the demo [Pluto](https://github.com/fonsp/Pluto.jl)
notebook, `notebooks/control_timing_safety.jl`, included in the package.

## Citing

If you use this code in your research, we ask that you consider citing the
relevant papers.  For the [`bounded_runs_iter`](@ref) function for
overapproximating the maximum deviation, please cite:

> Clara Hobbs, Bineet Ghosh, Shengjie Xu, Parasara Sridhar Duggirala, and
> Samarjit Chakraborty.
> "Safety Analysis of Embedded Controllers under Implementation Platform Timing
> Uncertainties."
> TCAD 2022.
> DOI: [10.1109/TCAD.2022.3198905](https://doi.org/10.1109/TCAD.2022.3198905)

For the [`estimate_deviation`](@ref) function, based on Jeffreys' Bayes factor
hypothesis testing, please cite:

> Bineet Ghosh, Clara Hobbs, Shengjie Xu, Parasara Sridhar Duggirala, James H.
> Anderson, P. S. Thiagarajan, and Samarjit Chakraborty.
> "Statistical Hypothesis Testing of Controller Implementations Under Timing
> Uncertainties."
> RTCSA 2022.
> DOI: [10.1109/RTCSA55878.2022.00008](https://doi.org/10.1109/RTCSA55878.2022.00008)

For the [`synthesize_constraints`](@ref) function to generate the list of safe weakly-hard constraints and [`schedule_xghtc`](@ref) function to synthesize schedules from weakly-hard constraints, please cite:
> Shengjie Xu, Bineet Ghosh, Clara Hobbs, P. S. Thiagarajan, and Samarjit 
> Chakraborty. 
> "Safety-Aware Flexible Schedule Synthesis for Cyber-Physical Systems Using
> Weakly-Hard Constraints." 
> ASPDAC 2023.
> DOI: [10.1145/3566097.3567848](https://doi.org/10.1145/3566097.3567848)

As we continue this line of work, we intend to incorporate new algorithms in
this package, and the list of papers will be updated accordingly (along with
specific the functions described in each).

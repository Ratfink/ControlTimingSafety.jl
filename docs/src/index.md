# ControlTimingSafety.jl

*Tools for verifying safety of control systems that may experience timing uncertainty*

This package, for the time being, consists primarily of an implementation of
the Bounded Runs Iteration algorithm from the paper "Safety Analysis of
Embedded Controllers under Implementation Platform Timing Uncertainties."  See
the [Citing](@ref) section below for more information.

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
length `n` and number of iterations `t`:

```julia
augbounds = augment(automaton, bounds)
all_bounds = bounded_runs_iter(automaton, augbounds, n, t)
```

The deviation at each time step can then be computed thusly:

```julia
d = deviation(automaton, augbounds, all_bounds)
```

For a complete example, see the demo [Pluto](https://github.com/fonsp/Pluto.jl)
notebook, `notebooks/control_timing_safety.jl`, included in the package.

## Citing

If you use this code in your research, we ask that you consider citing the
relevant paper.  Currently, this means the following:

> Clara Hobbs, Bineet Ghosh, Shengjie Xu, Parasara Sridhar Duggirala, and Samarjit Chakraborty.
> "Safety Analysis of Embedded Controllers under Implementation Platform Timing Uncertainties."
> TCAD 2022.

As we continue this line of work, we intend to incorporate new algorithms in
this package, and the list of papers will be updated accordingly (along with
specific the functions described in each).

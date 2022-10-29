# Automata

To represent the different control behaviors that can result from deadline
misses, we make use of finite-state *transducer automata*.  These take as input
a sequence of small integers representing deadline hits and misses, or
potentially other scheduler actions.

```@docs
Automaton
Automaton(::Automaton, ::Int64)
```

## Constructors for Miss Handling Strategies

As a user of ControlTimingSafety, you will likely not want to construct an
[`Automaton`](@ref) directly.  The package provides external constructors that
create an [`Automaton`](@ref) from a `StateSpace` model from
[ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl),
following different strategies of how to handle deadline misses.  Constraints
may be specified using `WeaklyHardConstraint` objects from
[RealTimeScheduling.jl](https://github.com/Ratfink/RealTimeScheduling.jl).

!!! note
    Currently, only `MissRow` constraints are supported, but more constraints
    will be added in the future!

```@docs
hold_kill
hold_skip_next
zero_kill
zero_skip_next
```

For convenience when building interfaces, e.g. notebooks, we also export a
sorted list of human-readable names for these constructors, and a map from
these names to the functions themselves.

```@docs
strat_names
strat_map
```

## Computing Evolution

Given an initial state `z0` and a sequence of scheduler actions `seq`, the
evolution of an [`Automaton`](@ref) can be computed.  Since the augmented state
space requires extra state variables outside the plant's state, the state `z0`
can be computed from a plant state `x0` using the [`augment`](@ref) function.

```@docs
augment
```

Once the initial augmented state is determined, we can then compute its
evolution.

```@docs
evol
evol_final
```

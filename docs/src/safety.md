# Checking Safety

The main safety-checking function of ControlTimingSafety is
[`bounded_runs_iter`](@ref), which implements the Bounded Runs Iteration
algorithm from "Safety Analysis of Embedded Controllers under Implementation
Platform Timing Uncertainties."  This implementation permits an initial point,
or an axis-aligned box initial set.  The helper function
[`bounded_runs`](@ref), which computes a single iteration, is also exported.

```@docs
bounded_runs_iter
bounded_runs
```

After computing a list of reachable sets, safety can be determined by computing
its distance from a nominal run.  This is done using the [`deviation`](@ref)
function.

```@docs
deviation
```

## Utility Functions

A few utility functions are exported as well.  These support the implementation
of the above functions, and are potentially useful to their callers as well.

```@docs
corners_from_bounds
merge_bounds
```

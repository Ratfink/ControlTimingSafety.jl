"""
    corners_from_bounds(bounds; cycle=false, dims=nothing)

Returns the corners of the n-dimensional interval represented by `bounds`.  If `cycle` is `true`, the first corner is repeated at the end, and the corners are given in Gray code order.  If `dims` is specified, only these dimensions are considered.
"""
function corners_from_bounds(bounds; cycle=false, dims=nothing)
    if dims == nothing
        dims = axes(bounds, 1)
    end
    ldims = length(dims)

    corners = cat(reshape([[c...] for c in Base.product(eachrow(bounds[dims,:])...)], 2^ldims)..., dims=2)
    if cycle
        gray(x) = x ⊻ (x >> 1)
        [corners[:,gray.(0:2^ldims-1) .+ 1] corners[:,1]]
    else
        corners
    end
end

_safefloatmin(x) = (isnan(x) || isinf(x)) ? +Inf : x
_safefloatmax(x) = (isnan(x) || isinf(x)) ? -Inf : x

"""
    merge_bounds(b)

Merges an array of bounding boxes `b` into one.
"""
function merge_bounds(b)
    mins = minimum(_safefloatmin, b[:,:,:,1], dims=1)
    maxs = maximum(_safefloatmax, b[:,:,:,2], dims=1)
    cat(dims=4, mins, maxs)[1,:,:,:]
end

"""
    bounded_runs(automaton, bounds, n)

Compute reachable sets for `n` time steps for the given `automaton`, starting from the initial set given by `bounds`.

`bounds` must be an `automaton.nz`×`2` matrix, whose first and second columns are the minimum and maximum along each dimension, respectively.

Returns `(bounds, locs)`, where `bounds` is an `automaton.A^n`×`n+1`×`automaton.nz`×`2` `Array{Float64}` giving the bounding box for each (run, time step), and `locs` is an `automaton.A^n` `Array{Int64}` of final locations for each run (or zero if there is no such run).
"""
function bounded_runs(automaton, bounds, n)
    corners = corners_from_bounds(bounds)

    # Stack
    z = Matrix{Float64}(undef, n+1, axes(corners)...)
    loc = Vector{Int64}(undef, n+1)
    act = Vector{Int64}(undef, n+1)

    # Bounding boxes for each final location, time step
    ret = Array{Float64}(undef, automaton.L, n+1, automaton.nz, 2) * NaN

    # Create the stack frame for time 0
    z[1,:,:] = corners
    loc[1] = automaton.l_int
    act[1] = 1
    # Initialize the stack pointer
    sp = 1
    # While we haven't popped all the way out
    while sp >= 1
        # If we've reached a leaf
        if sp == n+1
            # Calculate min and max for this final location at each time step
            ret[loc[sp],:,:,1] = minimum(_safefloatmin, cat(z, ret[loc[sp],:,:,1], dims=3), dims=3)
            ret[loc[sp],:,:,2] = maximum(_safefloatmax, cat(z, ret[loc[sp],:,:,2], dims=3), dims=3)
            sp -= 1
            # If we're out of actions from this step
        elseif act[sp] > automaton.A
            sp -= 1
            # If the transition is missing
        elseif ismissing(automaton.T[loc[sp], act[sp]])
            # Try the next transition
            act[sp] = act[sp] + 1
            # If the transition is present
        else
            z[sp+1,:] = automaton.Φ[automaton.μ[loc[sp], act[sp]]] * z[sp,:]
            loc[sp+1] = automaton.T[loc[sp], act[sp]]
            act[sp+1] = 1
            act[sp] = act[sp] + 1
            sp = sp + 1
        end
    end
    ret
end

"""
    bounded_runs_iter(automaton, bounds, n, t)

Iterate `bounded_runs(automaton, bounds, n)` for `t` iterations, returning the reachable set at each of the `n`×`t+1` time steps.
"""
function bounded_runs_iter(automaton, bounds, n, t)
    # Dimensions: time, augmented state, min/max
    all_bounds = Array{Float64}(undef, n*(t+1)+1, automaton.nz, 2)
    all_bounds[1,:,:] = bounds

    bounds = bounded_runs(automaton, bounds, n)
    all_bounds[1:n+1,:,:] = merge_bounds(bounds)[:,:,1:end]

    # Dimensions: initial location, final location, time, augmented state, min/max
    new_bounds = Array{Any}(undef, automaton.L, automaton.L, n+1, automaton.nz, 2)
    for i in 1:t
        # Simulate each box from previous iteration
        for i in 1:automaton.L
            a = Automaton_lint(automaton, i)
            new_bounds[i,:,:,:,:] = bounded_runs(a, bounds[i,end,:,:], n)
        end
        # Merge resulting boxes from these simulations
        for i in 1:automaton.L
            bounds[i,:,:,:] = merge_bounds(new_bounds[:,i,:,:,:])
        end
        # Save the bounds
        all_bounds[n*i+2:n*(i+1)+1,:,:] = merge_bounds(bounds)[2:end,:,:]
    end
    all_bounds
end

"""
    deviation(automaton, bounds, reachable; dims=[all], metric=Euclidean(), nominal=[2,2,...])

Compute the deviation from the `nominal` behavior (default: all `2`) that is possible for the given `automaton`, starting from the set of initial states `bounds`, within the `reachable` sets.  The deviation is computed using the specified `metric` from the [Distances.jl](https://www.juliapackages.com/p/distances) package.  If `dims` is specified, the deviation is computed for these dimensions only; otherwise, all dimensions are used.
"""
function deviation(automaton, bounds, reachable; dims=axes(bounds,1), metric=Euclidean(), nominal=repeat([2],size(reachable,1)-1))
    # Dimensions: state variables, points, time
    reachable_corners = cat([corners_from_bounds(reachable[t,:,:], dims=dims) for t in axes(reachable, 1)]..., dims=3)

    # Dimensions: state variables, points, time
    ev = Array{Float64}(undef, length(dims), 2^size(bounds,1), size(reachable, 1))
    corners = corners_from_bounds(bounds, dims=axes(bounds,1))
    for (i, c) in enumerate(eachcol(corners))
        e, _ = Evol(c, automaton, nominal)
        ev[:,i,:] = e'[dims,:]
    end

    # Compute Hausdorff distance at each time step
    H = Array{Float64}(undef, size(ev, 3))
    for t in axes(ev, 3)
        dist = pairwise(metric, reachable_corners[:,:,t], ev[:,:,t])
        H_row = maximum(minimum(dist, dims=1))
        H_col = maximum(minimum(dist, dims=2))
        H[t] = maximum((H_row, H_col))
    end
    H
end

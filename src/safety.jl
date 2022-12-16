"""
    corners_from_bounds(bounds::AbstractMatrix; cycle=false, dims=axes(bounds, 1))

Returns the corners of the n-dimensional interval represented by `bounds`.  If `cycle` is
`true`, the first corner is repeated at the end, and the corners are given in Gray code
order.  Only the dimensions from `dims` are considered.
"""
function corners_from_bounds(bounds::AbstractMatrix; cycle::Bool=false, dims=axes(bounds, 1))
    @boundscheck dims ⊆ axes(bounds, 1) || throw(ArgumentError("All entries of dims must be valid indices to the first dimension of bounds"))
    if size(bounds, 2) == 1
        return bounds
    end
    ldims = length(dims)

    corners = cat(reshape([[c...] for c in Base.product(eachrow(bounds[dims,:])...)],
			  2^ldims)..., dims=2)
    if cycle
        gray(x) = x ⊻ (x >> 1)
        [corners[:,gray.(0:2^ldims-1) .+ 1] corners[:,1]]
    else
        corners
    end
end

"""
    corners_from_bounds(bounds::AbstractVector; cycle=nothing, dims=nothing)

When applied to a vector, cast it to a one-column matrix.  `cycle` and `dims` are ignored.
"""
corners_from_bounds(bounds::AbstractVector; cycle=nothing, dims=nothing) = reshape(bounds, length(bounds), 1)

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
    bounded_runs(a::Automaton, z_0, n)

Compute reachable sets for `n` time steps for the given [`Automaton`](@ref) `a`, starting
from the initial set given by `z_0`.

`z_0` must be an `a.nz`-element vector, or an `a.nz`×`2` matrix whose first and second
columns are the minimum and maximum along each dimension, respectively.

Returns `(bounds, locs)`, where `bounds` is an `nactions(a)^n`×`n+1`×`a.nz`×`2`
`Array{Float64}` giving the bounding box for each (run, time step), and `locs` is an
`nactions(a)^n` `Array{Int64}` of final locations for each run (or zero if there is no such
run).

See also [`bounded_runs_iter`](@ref), which calls this function iteratively to efficiently
compute reachable sets for longer time horizons.  Typically one will call
[`deviation`](@ref) on the results of this function to determine deviation from a nominal
trajectory.
"""
function bounded_runs(a::Automaton, z_0::AbstractVecOrMat, n::Integer)
    corners = corners_from_bounds(z_0)

    # Stack
    # z gets one extra entry in the third dimension for cheap concatenation in leaf nodes
    z = Array{Float64}(undef, size(corners,1), size(corners,2)+1, n+1)
    loc = Vector{Int64}(undef, n+1)
    act = Vector{Int64}(undef, n+1)

    # Bounding boxes for each time step, final location
    ret = Array{Float64}(undef, a.nz, 2, n+1, nlocations(a))
    ret[:,1,:,:] .= Inf
    ret[:,2,:,:] .= -Inf

    # Create the stack frame for time 0
    z[:,begin:end-1,1] = corners
    loc[1] = a.l_int
    act[1] = 1
    # Initialize the stack pointer
    sp = 1
    # While we haven't popped all the way out
    while sp >= 1
        # If we've reached a leaf
        if sp == n+1
            # Calculate min and max for this final location at each time step
            z[:,end,:] = view(ret,:,1,:,loc[sp])
            minimum!(view(ret,:,1:1,:,loc[sp]), z, init=false)
            z[:,end,:] = view(ret,:,2,:,loc[sp])
            maximum!(view(ret,:,2:2,:,loc[sp]), z, init=false)
            sp -= 1
        # If we're out of actions from this step
        elseif act[sp] > nactions(a)
            sp -= 1
        # If the transition is missing
        elseif ismissing(a.T[loc[sp], act[sp]])
            # Try the next transition
            act[sp] += 1
        # If the transition is present
        else
            mul!(view(z,:,1:size(corners,2),sp+1), a.Φ[a.μ[loc[sp], act[sp]]], view(z,:,1:size(corners,2),sp))
            loc[sp+1] = a.T[loc[sp], act[sp]]
            act[sp+1] = 1
            act[sp] += 1
            sp += 1
        end
    end
    # TODO: the order of dimensions is largely an implementation detail, and since it was
    # non-optimal before, we should ultimately remove the need for this permutedims call
    permutedims(ret, [4, 3, 1, 2])
end

"""
    bounded_runs_iter(a, z_0, n, t)

Iterate [`bounded_runs`](@ref)`(a, z_0, n)` for `t` iterations, returning the reachable
set at each of the `n`×`t+1` time steps.

See also [`deviation`](@ref), which can be called with the result of this function to find
the deviation from a nominal trajectory.
"""
function bounded_runs_iter(a::Automaton, z_0::AbstractVecOrMat, n::Integer, t::Integer; safety_margin::Float64=Inf)
    # Dimensions: time, augmented state, min/max
    all_bounds = Array{Float64}(undef, n*(t+1)+1, a.nz, 2)
    if isa(z_0, AbstractVector)
        all_bounds[1,:,:] = [z_0 z_0]
    else
        all_bounds[1,:,:] = z_0
    end

    bounds = bounded_runs(a, z_0, n)
    all_bounds[1:n+1,:,:] = merge_bounds(bounds)[:,:,1:end]

    if isfinite(safety_margin)
        nom = evol(a, z_0, ones(n*(t+1)))
        d = deviation(a, z_0, all_bounds[1:n+1,:,:], nominal_trajectory=nom[1:n+1])
        if maximum(d) > safety_margin
            return all_bounds[1:n+1,:,:]
        end
    end

    # Dimensions: initial location, final location, time, augmented state, min/max
    new_bounds = Array{Float64}(undef, nlocations(a), nlocations(a), n+1, a.nz, 2)
    for i in 1:t
        # Simulate each box from previous iteration
        for i in a.L
            a = Automaton(a, i)
            new_bounds[i,:,:,:,:] = bounded_runs(a, bounds[i,end,:,:], n)
        end
        # Merge resulting boxes from these simulations
        for i in a.L
            bounds[i,:,:,:] = merge_bounds(new_bounds[:,i,:,:,:])
        end
        # Save the bounds
        all_bounds[n*i+2:n*(i+1)+1,:,:] = merge_bounds(bounds)[2:end,:,:]

        if isfinite(safety_margin)
            d = deviation(a, z_0, all_bounds[n*i+2:n*(i+1)+1,:,:], nominal_trajectory=nom[n*i+2:n*(i+1)+1])
            if maximum(d) > safety_margin
                return all_bounds[1:n*(i+1)+1,:,:]
            end
        end
    end
    all_bounds
end

"""
    deviation(a, z_0, reachable; metric=Euclidean(), nominal=[1,1,...])

Compute the deviation from the `nominal` behavior (default: all `1`) that is possible for
the given [`Automaton`](@ref) `a`, starting from the set of initial states `z_0`, within
the `reachable` sets.  The deviation is computed using the specified `metric` from the
[Distances.jl](https://www.juliapackages.com/p/distances) package, on the output of the
system (i.e. after multiplying by `a.C`).  `reachable` may be the three-dimensional output
of e.g. [`bounded_runs_iter`](@ref), or a matrix with dimensions (state, time), e.g. from
[`evol`](@ref).

If `nominal_trajectory` is given, this trajectory is used instead of computing the
trajectory from the `nominal` behavior.  This can improve efficiency when calling
`deviation` iteratively.

See also [`bounded_runs`](@ref) and [`bounded_runs_iter`](@ref), which can be used to
compute `reachable`.
"""
Base.@propagate_inbounds function deviation(a::Automaton, z_0::AbstractVecOrMat{Float64},
                   reachable::AbstractArray{Float64,3};
                   metric::PreMetric=Euclidean(),
                   nominal::AbstractVector{Int64}=ones(Int64,size(reachable,1)-1),
                   nominal_trajectory::Union{AbstractArray{Float64,3}, Nothing}=nothing)
    @boundscheck length(nominal) == size(reachable, 1) - 1 || throw(DimensionMismatch("nominal must have length size(reachable, 1) - 1"))

    # Dimensions: state variables, points, time
    reachable_corners = cat([corners_from_bounds(a.C * reachable[t,:,:]) for t in axes(reachable, 1)]..., dims=3)

    # Dimensions: state variables, points, time
    if z_0 isa AbstractVector{Float64}
        # XXX Not the most memory-efficient solution, but keeps us from having
        # to maintain two nearly-identical methods of the function.
        z_0 = [z_0 z_0]
    end
    if nominal_trajectory === nothing
        nominal_trajectory = Array{Float64}(undef, size(a.C, 1), 2^size(z_0,1), size(reachable, 1))
        corners = corners_from_bounds(z_0)
        for (i, c) in enumerate(eachcol(corners))
            e = evol(a, c, nominal)
            nominal_trajectory[:,i,:] = a.C * e'
        end
    end

    # Compute Hausdorff distance at each time step
    H = Array{Float64}(undef, size(nominal_trajectory, 3))
    for t in axes(nominal_trajectory, 3)
        dist = pairwise(metric, reachable_corners[:,:,t], nominal_trajectory[:,:,t])
        H_row = maximum(minimum(dist, dims=1))
        H_col = maximum(minimum(dist, dims=2))
        H[t] = max(H_row, H_col)
    end
    H
end

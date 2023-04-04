"""
    maximum_deviation_random(a, sampler, z_0, samples; estimate, metric, nominal,
                             nominal_trajectory)

Calculate `samples` random behaviors using `sampler`, and the corresponding trajectories of
the [`Automaton`](@ref) `a` from the initial state `z_0`.  Return the maximum deviation
from a nominal trajectory.  The parameters `metric`, `nominal`, and `nominal_trajectory` are
as in [`deviation`](@ref).  If `estimate` is specified, stop early if this deviation
estimate is exceeded, returning the exceeding deviation.
"""
Base.@propagate_inbounds function maximum_deviation_random(a::Automaton,
        sampler::RealTimeScheduling.SamplerWeaklyHard, samples::Int64,
        z_0::AbstractVecOrMat{Float64};
        estimate::Union{Float64, Nothing}=nothing, metric::PreMetric=Euclidean(),
        nominal::AbstractVector{Int64}=ones(Int64,sampler.H),
        nominal_trajectory::Union{AbstractArray{Float64}, Nothing}=nothing)
    @boundscheck length(nominal) == sampler.H || throw(DimensionMismatch("nominal must have length sampler.H"))
    @boundscheck nominal_trajectory === nothing || size(nominal_trajectory) == (size(z_0,1), sampler.H+1) || throw(DimensionMismatch("nominal_trajectory must have size (size(z_0,1), sampler.H+1)"))
    corners = unique(corners_from_bounds(z_0), dims=2)
    # Compute the nominal trajectory if not provided
    if nominal_trajectory === nothing
        nominal_trajectory = Array{Float64}(undef, size(z_0,1), sampler.H+1, size(corners,2))
        for i in axes(corners, 2)
            nominal_trajectory[:,:,i] = evol(a, corners[:,i], nominal)'
        end
    end
    Cnom = Array{Float64}(undef, size(a.C,1), sampler.H+1, size(corners,2))
    for i in axes(corners, 2)
        Cnom[:,:,i] = a.C * nominal_trajectory[:,:,i]
    end

    # Pre-allocate memory for things
    s = falses(sampler.H)
    beh = zeros(Int64, sampler.H)
    if z_0 isa AbstractVector
        z = zeros(sampler.H + 1, a.nz, 1)
        Cz = zeros(size(a.C, 1), sampler.H + 1, 1)
    else
        z = zeros(sampler.H + 1, a.nz, size(corners, 2))
        Cz = zeros(size(a.C, 1), sampler.H + 1, size(z, 3))
    end
    z[1,:,:] = corners
    dist = zeros(Float64, length(nominal)+1)

    H = zeros(sampler.H+1)
    dist = Array{Float64}(undef, size(corners,2), size(corners,2))
    r_row = zeros(1, size(corners,2))
    r_col = zeros(size(corners,2), 1)
    for _ in 1:samples
        # Generate a random behavior
        rand!(s, sampler)
        beh .= 2 .- s
        for i in axes(z, 3)
            # Compute its evolution
            evol!(a, view(z,:,:,i), beh)
            # Compute the output
            mul!(view(Cz,:,:,i), a.C, view(z,:,:,i)')
        end
        # Compute the distance between the output and nominal trajectory
        for t in eachindex(H)
            pairwise!(dist, metric, Cz[:,t,:], Cnom[:,t,:])
            H_row = maximum(minimum!(r_row, dist))
            H_col = maximum(minimum!(r_col, dist))
            H[t] = max(H_row, H_col, H[t])
        end
        # Exit early if the estimate is exceeded
        if estimate !== nothing && any(H .> estimate)
            return maximum(H)
        end
    end
    maximum(H)
end

"""
    estimate_deviation(a, sampler, z_0, c, B; nominal, nominal_trajectory,
                       estimate, estimate_samples, ϵ)

Estimate the deviation possible in the [`Automaton`](@ref) `a`, using `sampler` to generate
random behaviors.  The initial state is given as `z_0`.  `c` specifies the confidence with
which the returned value is asserted to be an upper bound, and `B` is the Bayes factor used
in hypothesis testing.  The parameters `metric`, `nominal`, and `nominal_trajectory` are as
in [`deviation`](@ref).

If `estimate` is given, this value is used as the initial estimate for the first hypothesis
formulated.  Otherwise, `estimate_samples` random runs (default: 50) are taken, and the
greatest deviation seen is used for the first hypothesis.  A pad of `ϵ` is added to each
subsequent hypothesis.

For more information, see Bineet Ghosh et al., "Statistical Hypothesis Testing of Controller
Implementations Under Timing Uncertainties," RTCSA 2022. 
DOI: [10.1109/RTCSA55878.2022.00008](https://doi.org/10.1109/RTCSA55878.2022.00008)
"""
Base.@propagate_inbounds function estimate_deviation(a::Automaton,
        sampler::RealTimeScheduling.SamplerWeaklyHard, z_0, c::Float64, B::Float64;
        metric::PreMetric=Euclidean(), nominal::AbstractVector{Int64}=ones(Int64,sampler.H),
        nominal_trajectory::Union{AbstractArray{Float64}, Nothing}=nothing,
        estimate::Union{Float64, Nothing}=nothing, estimate_samples::Int64=50, ϵ::Float64=0.)
    @boundscheck length(nominal) == sampler.H || throw(DimensionMismatch("nominal must have length sampler.H"))
    @boundscheck nominal_trajectory === nothing || size(nominal_trajectory) == (size(z_0,1), sampler.H+1) || throw(DimensionMismatch("nominal_trajectory must have size (size(z_0,1), sampler.H+1)"))

    # Compute the number of samples to draw for the required confidence and Bayes factor
    K = ceil(Int64, -log(c, c/(1-c)*B+1))
    if estimate === nothing
        # Take our initial estimate from a small random sample
        estimate = maximum_deviation_random(a, sampler, estimate_samples, z_0, metric=metric, nominal=nominal, nominal_trajectory=nominal_trajectory)
    end
    # Test the hypothesis, formulating new ones if we fail to accept it
    while true
        max_seen = maximum_deviation_random(a, sampler, K, z_0, metric=metric,
                                            estimate=estimate, nominal=nominal,
                                            nominal_trajectory=nominal_trajectory)
        if max_seen <= estimate
            break
        end
        estimate = max_seen + ϵ
    end
    estimate
end

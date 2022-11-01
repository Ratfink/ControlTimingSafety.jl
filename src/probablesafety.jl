"""
    maximum_deviation_random(a, sampler, z_0, samples; dims, estimate, nominal,
                             nominal_trajectory)

Calculate `samples` random behaviors using `sampler`, and the corresponding trajectories of
the [`Automaton`](@ref) `a` from the initial state `z_0`.  Return the maximum deviation
from a nominal trajectory.  The parameters `dims`, `nominal`, and `nominal_trajectory` are
as in [`deviation`](@ref).  If `estimate` is specified, stop early if this deviation
estimate is exceeded, returning the exceeding deviation.
"""
Base.@propagate_inbounds function maximum_deviation_random(a::Automaton,
        sampler::RealTimeScheduling.SamplerUniformMissRow, samples::Int64, z_0;
        dims=axes(z_0,1),
        estimate::Union{Float64, Nothing}=nothing,
        nominal::AbstractVector{Int64}=ones(Int64,sampler.H),
        nominal_trajectory::Union{AbstractArray{Float64}, Nothing}=nothing)
    @boundscheck length(nominal) == sampler.H || throw(DimensionMismatch("nominal must have length sampler.H"))
    @boundscheck dims ⊆ axes(z_0, 1) || throw(ArgumentError("All entries of dims must be valid indices to the first dimension of z_0"))
    @boundscheck nominal_trajectory === nothing || size(nominal_trajectory) == (size(z_0,1), sampler.H+1) || throw(DimensionMismatch("nominal_trajectory must have size (size(z_0,1), sampler.H+1)"))
    schedules = rand(sampler, samples)
    if nominal_trajectory === nothing
        nominal_trajectory = evol(a, z_0, nominal)
    end
    maxdev = 0.0
    for s in schedules
        ev = evol(a, z_0, 2 .- s)
        dist = deviation(a, z_0, ev, nominal_trajectory=nominal_trajectory, dims=dims)
        maxdev = maximum([dist; maxdev])
        if estimate !== nothing && maxdev > estimate
            return maxdev
        end
    end
    maxdev
end

"""
    estimate_deviation(a, sampler, z_0, c, B; dims, nominal, nominal_trajectory)

Estimate the deviation possible in the [`Automaton`](@ref) `a`, using `sampler` to generate
random behaviors.  The initial state is given as `z_0`.  `c` specifies the confidence with
which the returned value is asserted to be an upper bound, and `B` is the Bayes factor used
in hypothesis testing.  The parameters `dims`, `nominal`, and `nominal_trajectory` are as
in [`deviation`](@ref).

For more information, see Bineet Ghosh et al., "Statistical Hypothesis Testing of Controller
Implementations Under Timing Uncertainties," RTCSA 2022. 
DOI: [10.1109/RTCSA55878.2022.00008](https://doi.org/10.1109/RTCSA55878.2022.00008)
"""
Base.@propagate_inbounds function estimate_deviation(a::Automaton,
        sampler::RealTimeScheduling.SamplerUniformMissRow, z_0, c::Float64, B::Float64;
        dims=axes(z_0,1),
        nominal::AbstractVector{Int64}=ones(Int64,sampler.H),
        nominal_trajectory::Union{AbstractArray{Float64}, Nothing}=nothing)
    @boundscheck length(nominal) == sampler.H || throw(DimensionMismatch("nominal must have length sampler.H"))
    @boundscheck dims ⊆ axes(z_0, 1) || throw(ArgumentError("All entries of dims must be valid indices to the first dimension of z_0"))
    @boundscheck nominal_trajectory === nothing || size(nominal_trajectory) == (size(z_0,1), sampler.H+1) || throw(DimensionMismatch("nominal_trajectory must have size (size(z_0,1), sampler.H+1)"))
    K = ceil(Int64, -log(c, B+1))
    estimate = maximum_deviation_random(a, sampler, 20, z_0, dims=dims)
    while true
        max_seen = maximum_deviation_random(a, sampler, K, z_0, dims=dims,
                                            estimate=estimate, nominal=nominal,
                                            nominal_trajectory=nominal_trajectory)
        if max_seen <= estimate
            break
        end
        estimate = max_seen
    end
    estimate
end

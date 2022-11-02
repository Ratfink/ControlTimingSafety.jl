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
        sampler::RealTimeScheduling.SamplerUniformMissRow, samples::Int64, z_0;
        estimate::Union{Float64, Nothing}=nothing, metric::PreMetric=Euclidean(),
        nominal::AbstractVector{Int64}=ones(Int64,sampler.H),
        nominal_trajectory::Union{AbstractArray{Float64}, Nothing}=nothing)
    @boundscheck length(nominal) == sampler.H || throw(DimensionMismatch("nominal must have length sampler.H"))
    @boundscheck nominal_trajectory === nothing || size(nominal_trajectory) == (size(z_0,1), sampler.H+1) || throw(DimensionMismatch("nominal_trajectory must have size (size(z_0,1), sampler.H+1)"))
    # Compute the nominal trajectory if not provided
    if nominal_trajectory === nothing
        nominal_trajectory = evol(a, z_0, nominal)
    end
    Cnom = a.C * nominal_trajectory'

    # Pre-allocate memory for things
    s = falses(sampler.H)
    beh = zeros(Int64, sampler.H)
    z = zeros(sampler.H + 1, a.nz)
    z[1,:] = z_0
    Cz = zeros(size(a.C, 1), sampler.H + 1)
    dist = zeros(Float64, length(nominal)+1)

    maxdev = 0.0
    for _ in 1:samples
        # Generate a random behavior
        rand!(s, sampler)
        beh .= 2 .- s
        # Compute its evolution
        evol!(a, z, beh)
        # Compute the output
        mul!(Cz, a.C, z')
        # Compute the distance between the output and nominal trajectory
        colwise!(dist, metric, Cz, Cnom)
        # If it's the biggest we've seen, store it
        maxdev = max(maximum(dist), maxdev)
        # Bail out early if we exceeded the estimate
        if estimate !== nothing && maxdev > estimate
            return maxdev
        end
    end
    maxdev
end

"""
    estimate_deviation(a, sampler, z_0, c, B; nominal, nominal_trajectory)

Estimate the deviation possible in the [`Automaton`](@ref) `a`, using `sampler` to generate
random behaviors.  The initial state is given as `z_0`.  `c` specifies the confidence with
which the returned value is asserted to be an upper bound, and `B` is the Bayes factor used
in hypothesis testing.  The parameters `metric`, `nominal`, and `nominal_trajectory` are as
in [`deviation`](@ref).

For more information, see Bineet Ghosh et al., "Statistical Hypothesis Testing of Controller
Implementations Under Timing Uncertainties," RTCSA 2022. 
DOI: [10.1109/RTCSA55878.2022.00008](https://doi.org/10.1109/RTCSA55878.2022.00008)
"""
Base.@propagate_inbounds function estimate_deviation(a::Automaton,
        sampler::RealTimeScheduling.SamplerUniformMissRow, z_0, c::Float64, B::Float64;
        metric::PreMetric=Euclidean(), nominal::AbstractVector{Int64}=ones(Int64,sampler.H),
        nominal_trajectory::Union{AbstractArray{Float64}, Nothing}=nothing)
    @boundscheck length(nominal) == sampler.H || throw(DimensionMismatch("nominal must have length sampler.H"))
    @boundscheck nominal_trajectory === nothing || size(nominal_trajectory) == (size(z_0,1), sampler.H+1) || throw(DimensionMismatch("nominal_trajectory must have size (size(z_0,1), sampler.H+1)"))

    # Compute the number of samples to draw for the required confidence and Bayes factor
    K = ceil(Int64, -log(c, B+1))
    # Take our initial estimate from a small random sample
    estimate = maximum_deviation_random(a, sampler, 20, z_0)
    # Test the hypothesis, formulating new ones if we fail to accept it
    while true
        max_seen = maximum_deviation_random(a, sampler, K, z_0, metric=metric,
                                            estimate=estimate, nominal=nominal,
                                            nominal_trajectory=nominal_trajectory)
        if max_seen <= estimate
            break
        end
        estimate = max_seen
    end
    estimate
end

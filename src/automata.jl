import RealTimeScheduling

"""
    Automaton(Φ, T, μ, l_int)

A transducer automaton.

See also the constructors [`hold_kill`](@ref), [`hold_skip_next`](@ref),
[`zero_kill`](@ref), and [`zero_skip_next`](@ref), which create `Automaton` objects for
common deadline miss handling strategies.  For use in notebooks, etc., see
[`strat_map`](@ref) and [`strat_names`](@ref).
"""
struct Automaton
    """
    `Φ`: dynamics matrices (output alphabet).  Array of square matrices of size `nz`×`nz`.
    """
    Φ::Vector{AbstractMatrix{Float64}}
    """
    `T`: transition function.  `T[l,a]` is a location in `1:L`, or `missing` if no
    transition.
    """
    T::Matrix{Union{Missing, Int64}}
    """
    `μ`: output function.  `μ[l,a]` is an index into `Φ`, or `missing` if no transition.
    """
    μ::Matrix{Union{Missing, Int64}}
    """
    `l_int`: initial location in `1:L`.
    """
    l_int::Int64
    """
    `C`: Output matrix, i.e. `y[t] = C * z[t]`.
    """
    C::Matrix{Float64}

    function Automaton(Φ, T, μ, l_int, C)
        nz = size(Φ[1], 1)

        @boundscheck all(ϕ -> size(ϕ) == (nz, nz), Φ) || throw(DimensionMismatch("All matrices in Φ must be square and equal size"))
        @boundscheck size(T) == size(μ) || throw(DimensionMismatch("T and μ must be of equal size"))
        @boundscheck l_int ∈ axes(T, 1) || throw(ArgumentError("l_int must be a valid index to the first dimension of T"))
        @boundscheck all(t -> ismissing(t) || t ∈ axes(T, 1), T) || throw(ArgumentError("All entries of T must be valid indices to the first dimension of T or missing"))
        @boundscheck all(t -> ismissing(t) || t ∈ axes(Φ, 1), μ) || throw(ArgumentError("All entries of μ must be valid indices to Φ or missing"))
        @boundscheck size(C, 2) == nz || throw(DimensionMismatch("Rows of C must equal columns of all matrices in Φ"))

        new(Φ, T, μ, l_int, C)
    end
end

nlocations(a::Automaton) = size(a.T, 1)
nactions(a::Automaton) = size(a.T, 2)

function Base.getproperty(a::Automaton, s::Symbol)
    if s === :L
        return Base.oneto(nlocations(a))
    elseif s === :A
        return Base.oneto(nactions(a))
    elseif s === :nz
        return size(a.Φ[1], 1)
    else
        return getfield(a, s)
    end
end

function Base.propertynames(a::Automaton, private::Bool=false)
    (fieldnames(typeof(a))..., :L, :A, :nz)
end

"""
    Automaton(a, l_int)

Construct a copy of the [`Automaton`](@ref) `a` with a new `l_int`.
"""
function Automaton(a::Automaton, l_int::Int64)
    @boundscheck l_int ∈ a.L || throw(ArgumentError("l_int must be a valid location in a.L"))

    @inbounds Automaton(a.Φ, a.T, a.μ, l_int, a.C)
end

"""
    hold_skip_next(sysd, K, [c])

Construct an [`Automaton`](@ref) for the given discrete-time state space model `sysd` with
feedback gain matrix `K`, following the Hold&Skip-Next strategy.  If `c` is specified, the
resulting `Automaton` will permit only that weakly hard constraint.
"""
function hold_skip_next(sysd::AbstractStateSpace{<:Discrete},
        K::AbstractMatrix{Float64})
    _hold_skip_next(sysd, K, _skip_next()...)
end

function hold_skip_next(sysd::AbstractStateSpace{<:Discrete},
        K::AbstractMatrix{Float64}, c::RealTimeScheduling.WeaklyHardConstraint)
    _hold_skip_next(sysd, K, _skip_next(c)...)
end

function _hold_skip_next(sysd::AbstractStateSpace{<:Discrete},
        K::AbstractMatrix{Float64}, T::AbstractMatrix{Union{Int64,Missing}},
        μ::AbstractMatrix{Union{Int64,Missing}})
    A, B, C, _ = ssdata(sysd)
    p, r = size(B)

    # Split apart the pieces of K, if necessary
    K_x = -K[:,1:p]
    if size(K, 2) == p + r
        K_u = -K[:,p+1:p+r]
    else
        K_u = zeros((r, r))
    end

    # Put it all together, with the Φ matrices
    Automaton([# Φ_HH
               [A  zeros(p, p)  B;
                zeros(p, 2p + r);
                K_x  zeros(r, p)  K_u],
               # Φ_MH
               [A  zeros(p, p)  B;
                zeros(p, 2p + r);
                zeros(r, p)  K_x  K_u],
               # Φ_HM
               [A  zeros(p, p)  B;
                I(p)  zeros(p, p + r);
                zeros(r, 2p)  I(r)],
               # Φ_MM
               [A  zeros(p, p)  B;
                zeros(p, p)  I(p)  zeros(p, r);
                zeros(r, 2p)  I(r)]],
              T, μ, 1, [C zeros(size(C, 1), p + r)])
end

function _skip_next()
    T = Union{Missing,Int64}[1 2
                             1 2]
    μ = Union{Missing,Int64}[1 3
                             2 4]
    return T, μ
end

function _skip_next(c::RealTimeScheduling.MissRow)
    if c.miss == 0
        T = [1 missing]
        μ = [1 missing]
    else
        T = [ones(Int64, (c.miss+1, 1)) [collect(2:c.miss+1); missing]]
        μ = [1 3;
             repeat([2 4], c.miss-1, 1);
             2 missing]
    end
    return T, μ
end

function _skip_next(c::RealTimeScheduling.MeetAny)
    if c.meet == c.window
        T, μ = _skip_next(RealTimeScheduling.MissRow(0))
    elseif c.meet == 0
        T, μ = _skip_next()
    else
        L = 2^c.window
		T = zeros(Union{Missing, Int64}, (L, 2))
		μ = zeros(Union{Missing, Int64}, (L, 2))

        for i = 0:L-1
            # hit
			T[i+1, 1] = ((i << 1) & (L - 1) | 1) + 1
            # miss
			T[i+1, 2] = ((i << 1) & (L - 1)) + 1

            # hit
			μ[i+1, 1] = 2 - (i & 1)
            # miss
			μ[i+1, 2] = 4 - (i & 1)

            # Check for constraint
			for j = 1:2
				if count_ones(T[i+1, j] - 1) < c.meet
					T[i+1, j] = missing
					μ[i+1, j] = missing
				end
			end
		end

        # TODO: figure out how to remove this
		T = L + 1 .- T
		reverse!(T, dims=1)
		reverse!(μ, dims=1)
    end
    return T, μ
end

"""
    zero_skip_next(sysd, K, [c])

Construct an [`Automaton`](@ref) for the given discrete-time state space model `sysd` with
feedback gain matrix `K`, following the Zero&Skip-Next strategy.  If `c` is specified, the
resulting `Automaton` will permit only that weakly hard constraint.
"""
function zero_skip_next(sysd::AbstractStateSpace{<:Discrete},
        K::AbstractMatrix{Float64})
    _zero_skip_next(sysd, K, _skip_next()...)
end

function zero_skip_next(sysd::AbstractStateSpace{<:Discrete},
        K::AbstractMatrix{Float64}, c::RealTimeScheduling.WeaklyHardConstraint)
    _zero_skip_next(sysd, K, _skip_next(c)...)
end

function _zero_skip_next(sysd::AbstractStateSpace{<:Discrete},
        K::AbstractMatrix{Float64}, T::AbstractMatrix{Union{Int64,Missing}},
        μ::AbstractMatrix{Union{Int64,Missing}})
    A, B, C, _ = ssdata(sysd)
    p, r = size(B)

    # Split apart the pieces of K, if necessary
    K_x = -K[:,1:p]
    if size(K, 2) == p + r
        K_u = -K[:,p+1:p+r]
    else
        K_u = zeros((r, r))
    end

    # Put it all together, with the Φ matrices
    Automaton([# Φ_HH
               [A  zeros(p, p)  B;
                zeros(p, 2p + r);
                K_x  zeros(r, p)  K_u],
               # Φ_MH
               [A  zeros(p, p)  B;
                zeros(p, 2p + r);
                zeros(r, p)  K_x  K_u],
               # Φ_HM
               [A  zeros(p, p)  B;
                I(p)  zeros(p, p + r);
                zeros(r, 2p + r)],
               # Φ_MM
               [A  zeros(p, p)  B;
                zeros(p, p)  I(p)  zeros(p, r);
                zeros(r, 2p + r)]],
              T, μ, 1, [C zeros(size(C, 1), p + r)])
end

"""
    hold_kill(sysd, K, [c])

Construct an [`Automaton`](@ref) for the given discrete-time state space model `sysd` with
feedback gain matrix `K`, following the Hold&Kill strategy.  If `c` is specified, the
resulting `Automaton` will permit only that weakly hard constraint.
"""
function hold_kill(sysd::AbstractStateSpace{<:Discrete}, K::AbstractMatrix{Float64})
    _hold_kill(sysd, K, _kill()...)
end

function hold_kill(sysd::AbstractStateSpace{<:Discrete}, K::AbstractMatrix{Float64},
        c::RealTimeScheduling.WeaklyHardConstraint)
    _hold_kill(sysd, K, _kill(c)...)
end

function _hold_kill(sysd::AbstractStateSpace{<:Discrete}, K::AbstractMatrix{Float64},
        T::AbstractMatrix{Union{Int64,Missing}}, μ::AbstractMatrix{Union{Int64,Missing}})
    A, B, C, _ = ssdata(sysd)
    p, r = size(B)

    # Split apart the pieces of K, if necessary
    K_x = -K[:,1:p]
    if size(K, 2) == p + r
        K_u = -K[:,p+1:p+r]
    else
        K_u = zeros((r, r))
    end

    # Put it all together, with the Φ matrices
    Automaton([# Φ_H
               [A  B;
                K_x  K_u],
               # Φ_M
               [A  B;
                zeros(r, p)  I(r)]],
              T, μ, 1, [C zeros(size(C, 1), r)])
end

function _kill()
    T = Union{Missing,Int64}[1 1]
    μ = Union{Missing,Int64}[1 2]
    return T, μ
end

function _kill(c::RealTimeScheduling.MissRow)
    if c.miss == 0
        T = [1 missing]
        μ = [1 missing]
    else
        T = [ones(Int64, (c.miss+1, 1)) [collect(2:c.miss+1); missing]]
        μ = [1 2;
             repeat([1 2], c.miss-1, 1);
             1 missing]
    end
    return T, μ
end

function _kill(c::RealTimeScheduling.MeetAny)
    if c.meet == c.window
        T, μ = _kill(RealTimeScheduling.MissRow(0))
    elseif c.meet == 0
        T, μ = _kill()
    else
		L = 2^c.window
		T = zeros(Union{Missing, Int64}, (L, 2))
		μ = zeros(Union{Missing, Int64}, (L, 2))
		for i = 0:L-1
            # hit
			T[i+1, 1] = ((i << 1) & (L - 1) | 1) + 1
            # miss
			T[i+1, 2] = ((i << 1) & (L - 1)) + 1

			# hit
            μ[i+1, 1] = 1
            # miss
			μ[i+1, 2] = 2

            # Check for constraint
			for j = 1:2
				if count_ones(T[i+1, j] - 1) < c.meet
					T[i+1, j] = missing
					μ[i+1, j] = missing
				end
			end
		end
        
        # TODO: figure out how to remove this
		T = L + 1 .- T
		reverse!(T, dims=1)
		reverse!(μ, dims=1)
    end
    return T, μ
end

"""
    zero_kill(sysd, K, [c])

Construct an [`Automaton`](@ref) for the given discrete-time state space model `sysd` with
feedback gain matrix `K`, following the Zero&Kill strategy.  If `c` is specified, the
resulting `Automaton` will permit only that weakly hard constraint.
"""
function zero_kill(sysd::AbstractStateSpace{<:Discrete}, K::AbstractMatrix{Float64})
    _zero_kill(sysd, K, _kill()...)
end

function zero_kill(sysd::AbstractStateSpace{<:Discrete}, K::AbstractMatrix{Float64},
        c::RealTimeScheduling.WeaklyHardConstraint)
    _zero_kill(sysd, K, _kill(c)...)
end

function _zero_kill(sysd::AbstractStateSpace{<:Discrete}, K::AbstractMatrix{Float64},
        T::AbstractMatrix{Union{Int64,Missing}}, μ::AbstractMatrix{Union{Int64,Missing}})
    A, B, C, _ = ssdata(sysd)
    p, r = size(B)

    # Split apart the pieces of K, if necessary
    K_x = -K[:,1:p]
    if size(K, 2) == p + r
        K_u = -K[:,p+1:p+r]
    else
        K_u = zeros((r, r))
    end

    # Put it all together, with the Φ matrices
    Automaton([# Φ_H
               [A  B;
                K_x  K_u],
               # Φ_M
               [A  B;
                zeros(r, p + r)]],
              T, μ, 1, [C zeros(size(C, 1), r)])
end

"""
A mapping from human-readable names to [`Automaton`](@ref) constructor functions.

See also [`strat_names`](@ref).
"""
strat_map = Dict(
    "Hold&Skip-Next" => hold_skip_next,
    "Zero&Skip-Next" => zero_skip_next,
    "Hold&Kill" => hold_kill,
    "Zero&Kill" => zero_kill
)
"""
A stable, sorted list of the human-readable names in [`strat_map`](@ref).
"""
strat_names = sort([keys(strat_map)...])


"""
    evol_final(a, z_0, input)

Same as [`evol`](@ref), but returns `(z, l)`, where `z` is a matrix of states over time,
and `l` is the final location in the automaton.
"""
function evol_final(a::Automaton, z_0::AbstractVector{Float64}, input::AbstractVector{Int64})
    @boundscheck length(z_0) == a.nz || throw(DimensionMismatch("z_0 must have length a.nz"))
    z = zeros(length(input) + 1, a.nz)
    z[1,:] = z_0
    evol_final!(a, z, input)
end

"""
    evol_final!(a, z, input)

Same as [`evol_final`](@ref), but writes to the input matrix `z`, whose first row is `z_0`.
"""
function evol_final!(a::Automaton, z::AbstractMatrix{Float64}, input::AbstractVector{Int64})
    @boundscheck size(z,1) == length(input)+1 || throw(DimensionMismatch("z must have size (length(input)+1, a.nz)"))
    @boundscheck size(z,2) == a.nz || throw(DimensionMismatch("z must have size (length(input)+1, a.nz)"))
    l = a.l_int
    # For each time step and input action
    for (t, in) in enumerate(input)
        # Get the dynamics matrix
        μ = a.μ[l, in]
        # If we hit a missing transition, return the states that we reached,
        # and a missing final location to signal the problem to the caller.
        if ismissing(μ)
            return z[1:t,:], missing
        end
        # Apply the dynamics
        mul!(view(z,t+1,:), a.Φ[μ], view(z,t,:))
        # Transition to the new location
        l = a.T[l, in]
    end
    z, l
end

"""
    evol(a, z_0, input)

Simulate a run of [`Automaton`](@ref) `a` from the initial state `z_0`, using the sequence
of scheduler events `input`.

Returns `z`, a matrix of states over time.

See also [`evol_final`](@ref), which additionally returns the final location in the
automaton.
"""
evol(a::Automaton, z_0::AbstractVector{Float64}, input::AbstractVector{Int64}) = evol_final(a, z_0, input)[1]

"""
    evol!(a, z, input)

Same as [`evol`](@ref), but writes to the input matrix `z`, whose first row is `z_0`.
"""
evol!(a::Automaton, z::AbstractMatrix{Float64}, input::AbstractVector{Int64}) = evol_final!(a, z, input)[1]

"""
    augment(a::Automaton, x)

Pad the state vector (or matrix of column vectors) `x` to as many dimensions as in the
augmented state space of `a`.
"""
augment(a::Automaton, x::AbstractMatrix) = [x; zeros(a.nz - size(x, 1), size(x, 2))]
augment(a::Automaton, x::AbstractVector) = [x; zeros(a.nz - size(x, 1))]

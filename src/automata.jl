
"""
    Automaton(Φ, T, μ, l_int)

A transducer automaton.

# Fields

* `Φ`: dynamics matrices (output alphabet).  Array of square matrices of size `nz`×`nz`.
* `T`: transition function.  `T[l,a]` is a location in `1:L`, or `missing` if no transition.
* `μ`: output function.  `μ[l,a]` is an index into `Φ`, or `missing` if no transition.
* `l_int`: initial location in `1:L`.
"""
struct Automaton
    Φ::Vector{AbstractMatrix{Float64}}
    T::Matrix{Union{Missing, Int64}}
    μ::Matrix{Union{Missing, Int64}}
    l_int::Int64

    function Automaton(Φ, T, μ, l_int)
        nz = size(Φ[1], 1)

        @boundscheck all(t -> t == (nz, nz), size.(Φ)) || throw(DimensionMismatch("All matrices in Φ must be square and equal size"))
        @boundscheck size(T) == size(μ) || throw(DimensionMismatch("T and μ must be of equal size"))
        @boundscheck l_int ∈ axes(T, 1) || throw(ArgumentError("l_int must be a valid index to the first dimension of T"))
        @boundscheck all(t -> t == 0 || t ∈ axes(T, 1), T) || throw(ArgumentError("All entries of T must be valid indices to the first dimension of T"))
        @boundscheck all(t -> t == 0 || t ∈ axes(Φ, 1), μ) || throw(ArgumentError("All entries of μ must be valid indices to Φ"))

        new(Φ, T, μ, l_int)
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
    Automaton_lint(a, l_int)

Construct a copy of the `Automaton` `a` with a new `l_int`.
"""
function Automaton_lint(a::Automaton, l_int::Int64)
    @boundscheck l_int ∈ a.L || throw(ArgumentError("l_int must be a valid location in a.L"))

    @inbounds Automaton(a.Φ, a.T, a.μ, l_int)
end

"""
    hold_skip_next(sysd, K, miss=nothing)

Construct an `Automaton` for the given discrete-time state space model `sysd` with feedback gain matrix `K`, following the Hold&Skip-Next strategy, with at most `miss` deadline misses.
"""
function hold_skip_next(sysd, K, miss=nothing)
    p, r = size(sysd.B)

    # Split apart the pieces of K, if necessary
    K_x = -K[:,1:p]
    if size(K, 2) == p + r
        K_u = -K[:,p+1:p+r]
    else
        K_u = zeros((r, r))
    end

    # Define the automaton's structure
    if miss === nothing || miss == -1
        T = [1 2;
             1 2]
        μ = [1 3;
             2 4]
    elseif miss == 0
        T = [1 missing]
        μ = [1 missing]
    else
        T = [ones(Int64, (miss+1, 1)) [collect(2:miss+1); missing]]
        μ = [1 3;
             repeat([2 4], miss-1, 1);
             2 missing]
    end

    # Put it all together, with the Φ matrices
    Automaton([# Φ_HH
               [sysd.A  zeros(p, p)  sysd.B;
                zeros(p, 2p + r);
                K_x  zeros(r, p)  K_u],
               # Φ_MH
               [sysd.A  zeros(p, p)  sysd.B;
                zeros(p, 2p + r);
                zeros(r, p)  K_x  K_u],
               # Φ_HM
               [sysd.A  zeros(p, p)  sysd.B;
                I(p)  zeros(p, p + r);
                zeros(r, 2p)  I(r)],
               # Φ_MM
               [sysd.A  zeros(p, p)  sysd.B;
                zeros(p, p)  I(p)  zeros(p, r);
                zeros(r, 2p)  I(r)]],
              T, μ, 1)
end

"""
    zero_skip_next(sysd, K, miss=nothing)

Construct an `Automaton` for the given discrete-time state space model `sysd` with feedback gain matrix `K`, following the Zero&Skip-Next strategy, with at most `miss` deadline misses.
"""
function zero_skip_next(sysd, K, miss=nothing)
    p, r = size(sysd.B)

    # Split apart the pieces of K, if necessary
    K_x = -K[:,1:p]
    if size(K, 2) == p + r
        K_u = -K[:,p+1:p+r]
    else
        K_u = zeros((r, r))
    end

    # Define the automaton's structure
    if miss === nothing || miss == -1
        T = [1 2;
             1 2]
        μ = [1 3;
             2 4]
    elseif miss == 0
        T = [1 missing]
        μ = [1 missing]
    else
        T = [ones(Int64, (miss+1, 1)) [collect(2:miss+1); missing]]
        μ = [1 3;
             repeat([2 4], miss-1, 1);
             2 missing]
    end

    # Put it all together, with the Φ matrices
    Automaton([# Φ_HH
               [sysd.A  zeros(p, p)  sysd.B;
                zeros(p, 2p + r);
                K_x  zeros(r, p)  K_u],
               # Φ_MH
               [sysd.A  zeros(p, p)  sysd.B;
                zeros(p, 2p + r);
                zeros(r, p)  K_x  K_u],
               # Φ_HM
               [sysd.A  zeros(p, p)  sysd.B;
                I(p)  zeros(p, p + r);
                zeros(r, 2p + r)],
               # Φ_MM
               [sysd.A  zeros(p, p)  sysd.B;
                zeros(p, p)  I(p)  zeros(p, r);
                zeros(r, 2p + r)]],
              T, μ, 1)
end

"""
    hold_kill(sysd, K, miss=nothing)

Construct an `Automaton` for the given discrete-time state space model `sysd` with feedback gain matrix `K`, following the Hold&Kill strategy, with at most `miss` deadline misses.
"""
function hold_kill(sysd, K, miss=nothing)
    p, r = size(sysd.B)

    # Split apart the pieces of K, if necessary
    K_x = -K[:,1:p]
    if size(K, 2) == p + r
        K_u = -K[:,p+1:p+r]
    else
        K_u = zeros((r, r))
    end

    # Define the automaton's structure
    if miss === nothing || miss == -1
        T = [1 1]
        μ = [1 2]
    elseif miss == 0
        T = [1 missing]
        μ = [1 missing]
    else
        T = [ones(Int64, (miss+1, 1)) [collect(2:miss+1); missing]]
        μ = [1 2;
             repeat([1 2], miss-1, 1);
             1 missing]
    end

    # Put it all together, with the Φ matrices
    Automaton([# Φ_H
               [sysd.A  sysd.B;
                K_x  K_u],
               # Φ_M
               [sysd.A  sysd.B;
                zeros(r, p)  I(r)]],
              T, μ, 1)
end

"""
    zero_kill(sysd, K, miss=nothing)

Construct an `Automaton` for the given discrete-time state space model `sysd` with feedback gain matrix `K`, following the Zero&Kill strategy, with at most `miss` deadline misses.
"""
function zero_kill(sysd, K, miss=nothing)
    p, r = size(sysd.B)

    # Split apart the pieces of K, if necessary
    K_x = -K[:,1:p]
    if size(K, 2) == p + r
        K_u = -K[:,p+1:p+r]
    else
        K_u = zeros((r, r))
    end

    # Define the automaton's structure
    if miss === nothing || miss == -1
        T = [1 1]
        μ = [1 2]
    elseif miss == 0
        T = [1 missing]
        μ = [1 missing]
    else
        T = [ones(Int64, (miss+1, 1)) [collect(2:miss+1); missing]]
        μ = [1 2;
             repeat([1 2], miss-1, 1);
             1 missing]
    end

    # Put it all together, with the Φ matrices
    Automaton([# Φ_H
               [sysd.A  sysd.B;
                K_x  K_u],
               # Φ_M
               [sysd.A  sysd.B;
                zeros(r, p + r)]],
              T, μ, 1)
end

strat_map = Dict(
    "Hold&Skip-Next" => hold_skip_next,
    "Zero&Skip-Next" => zero_skip_next,
    "Hold&Kill" => hold_kill,
    "Zero&Kill" => zero_kill
)
strat_names = sort([keys(strat_map)...])


"""
    evol_final(a, z_0, input)

Same as `evol`, but returns `(z, l)`, where `z` is a matrix of states over time, and `l` is the final location in the automaton.
"""
function evol_final(a::Automaton, z_0, input)
    z = zeros(a.nz, length(input) + 1)
    z[:,1] = z_0
    l = a.l_int
    # For each time step and input action
    for (t, in) in enumerate(input)
        # Get the dynamics matrix
        μ = a.μ[l, in]
        # If we hit a missing transition, return the states that we reached,
        # and a missing final location to signal the problem to the caller.
        if ismissing(μ)
            return z[:,1:t]', missing
        end
        # Apply the dynamics
        z[:,t+1] = a.Φ[μ] * z[:,t]
        # Transition to the new location
        l = a.T[l, in]
    end
    z', l
end

"""
    evol(a, z_0, input)

Simulate a run of automaton `a` from the initial state `z_0`, using the sequence of scheduler events `input`.

Returns `z`, a matrix of states over time.

See also `evol_final`, which additionally returns the final location in the automaton.
"""
evol(a::Automaton, z_0, input) = evol_final(a, z_0, input)[1]

"""
    augment(a, x)

Pad the state vector (or matrix of column vectors) `x` to as many dimensions as in the augmented state space of `a`.
"""
augment(a::Automaton, x) = [x; zeros(a.nz - size(x, 1), size(x, 2))]

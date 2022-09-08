
"""
    Automaton(Φ, T, μ, l_int)

A transducer automaton.

# Fields

* `Φ`: dynamics matrices (output alphabet).  Array of square matrices of size `nz`×`nz`.
* `T`: transition function.  `T[l,a]` is a location in `1:L`, or `missing` if no transition.
* `μ`: output function.  `μ[l,a]` is an index into `Φ`, or `missing` if no transition.
* `l_int`: initial location in `1:L`.
* `L`: locations.  Legal locations are in the range `1:L`.
* `A`: actions (input alphabet).  Legal actions are in the range `1:A`.
* `nz`: number of dimensions in the augmented state space.
"""
struct Automaton
    Φ::Vector{AbstractMatrix{Float64}}
    T::Matrix{Union{Missing, Int64}}
    μ::Matrix{Union{Missing, Int64}}
    l_int::Int64
    L::Int64
    A::Int64
    nz::Int64

    function Automaton(Φ, T, μ, l_int)
        nz = size(Φ[1], 1)

        @boundscheck all(t -> t == (nz, nz), size.(Φ)) || throw(DimensionMismatch("All matrices in Φ must be square and equal size"))
        @boundscheck size(T) == size(μ) || throw(DimensionMismatch("T and μ must be of equal size"))
        @boundscheck l_int ∈ axes(T, 1) || throw(ArgumentError("l_int must be a valid index to the first dimension of T"))
        @boundscheck all(t -> t == 0 || t ∈ axes(T, 1), T) || throw(ArgumentError("All entries of T must be valid indices to the first dimension of T"))
        @boundscheck all(t -> t == 0 || t ∈ axes(Φ, 1), μ) || throw(ArgumentError("All entries of μ must be valid indices to Φ"))

        new(Φ, T, μ, l_int, size(T)..., nz)
    end
end

"""
    Automaton_lint(a, l_int)

Construct a copy of the `Automaton` `a` with a new `l_int`.
"""
function Automaton_lint(a::Automaton, l_int::Int64)
    @boundscheck l_int ∈ axes(a.T, 1) || throw(ArgumentError("l_int must be a valid index to the first dimension of T"))

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
        T = [2 1;
             2 1]
        μ = [3 1;
             4 2]
    elseif miss == 0
        T = [missing 1]
        μ = [missing 1]
    else
        T = [[collect(2:miss+1); missing] ones(Int64, (miss+1, 1))]
        μ = [3 1;
             repeat([4 2], miss-1, 1);
             missing 2]
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
        T = [2 1;
             2 1]
        μ = [3 1;
             4 2]
    elseif miss == 0
        T = [missing 1]
        μ = [missing 1]
    else
        T = [[collect(2:miss+1); missing] ones(Int64, (miss+1, 1))]
        μ = [3 1;
             repeat([4 2], miss-1, 1);
             missing 2]
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
        μ = [2 1]
    elseif miss == 0
        T = [missing 1]
        μ = [missing 1]
    else
        T = [[collect(2:miss+1); missing] ones(Int64, (miss+1, 1))]
        μ = [2 1;
             repeat([2 1], miss-1, 1);
             missing 1]
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
        μ = [2 1]
    elseif miss == 0
        T = [missing 1]
        μ = [missing 1]
    else
        T = [[collect(2:miss+1); missing] ones(Int64, (miss+1, 1))]
        μ = [2 1;
             repeat([2 1], miss-1, 1);
             missing 1]
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

"""
    evol_final(a, z_0, input)

Same as `evol`, but returns `(z, l)`, where `z` is a matrix of states over time, and `l` is the final location in the automaton.
"""
function evol_final(a::Automaton, z_0, input)
    t_max = length(input)
    z = zeros(size(z_0, 1), t_max + 1)
    z[:,1] = z_0
    l = a.l_int
    # For each time step
    for t = 1:t_max
        # Get the dynamics matrix
        μ = a.μ[l, input[t]]
        # If we hit a missing transition, return the states that we reached,
        # and a missing final location to signal the problem to the caller.
        if ismissing(μ)
            return z[:,1:t]', missing
        end
        # Apply the dynamics
        z[:,t+1] = a.Φ[μ] * z[:,t]
        # Transition to the new location
        l = a.T[l, input[t]]
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

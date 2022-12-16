"""
    schedule_xghtc(constraints; slotsize=1, H=100)

Generate a schedule for a set of weakly hard constraints. The schedule assumes that all tasks are synchronous and have equal periods. At most `slotsize` tasks may be scheduled in a single period. The schedule has total length `H` periods.
    
Shengjie Xu, Bineet Ghosh, Clara Hobbs, P.S. Thiagarajan, and Samarjit Chakraborty, "Safety-Aware Flexible Schedule Synthesis for Cyber-Physical Systems using Weakly-Hard Constraints." ASP-DAC 2023.
"""
function schedule_xghtc(constraints::Vector{RealTimeScheduling.WeaklyHardConstraint}; slotsize::Int64=1, H::Int64=100)
    # Check if "utilization" is greater than available slot size
    utilization = sum(c -> c.meet/c.window, constraints)
    if utilization > slotsize
        return Vector{Int64}()
    end

    # Create the scheduler automaton from individual constraints
    controllers = _controller_automaton.(constraints)
    a = _scheduler_automaton(controllers, slotsize=slotsize)

    # Initialize the list of current states from the initial state of
    # scheduler automaton
    current_states = Dict{Int64, Cons{Int64}}(a.l_0 => list(a.l_0))

    # Traverse the tree until 
    #   (1) the required number of steps is reached,
    #   (2) a cycle is found, or 
    #   (3) no more states are valid for exploration
    for step in 1:steps
        next_states = Dict{Int64, Cons{Int64}}()
        for (l, path) in pairs(current_states), σ in a.Σ
            l_new = a.T(l, σ)
            if !a.Q_f(l_new) || haskey(next_states, l_new)
                continue
            elseif l_new in path
                # Found cycle -> Case (2)
                return reverse(cons(l_new, path))
            end
            next_states[l_new] = cons(l_new, path)
        end

        if isempty(next_states)
            # No more valid states -> Case (3)
            return Vector{Int64}()
        end

        current_states = next_states
    end

    for (l, path) in pairs(current_states)
        if a.Q_f(l)
            # If the outer loop ends and accepting state is found -> Case (1)
            return reverse(path)
        end
    end

    # If the outer loop ends and accepting state is not found -> Case (1)
    Vector{Int64}()
end

"""
    function synthesize_constraints(sysd, K, z_0, d_max, maxwindow, n, t)

Find all `MeetAny` weakly hard constraints with window size at most `maxwindow` that guarantees the deviation upper bound is at most `d_max`. The system is specified by [`Automaton`](@ref) `a` and initial state is `z_0`. `n` and `t` are as in [`bounded_runs_iter`](@ref).
"""
function synthesize_constraints(sysd::AbstractStateSpace{<:Discrete},
    K::AbstractMatrix{Float64}, z_0::AbstractVecOrMat, d_max::Float64,
    maxwindow::Int64, n::Int64, t::Int64)

    # TODO: 
    #   1) Test to make sure it works the same
    #   2) Discuss: should additional safe constraints be included in result? 
    #      Or only the ones that are actually checked & useful for scheduling
    safe_constraints = RealTimeScheduling.MeetAny[]

    # Do not need to go through all O(maxwindow^2) constraints,
    # see paper for optimization argument
    meet = 1
    for window in 2:maxwindow
        while meet < window
            constraint = RealTimeScheduling.MeetAny(meet, window)
            a = hold_kill(sysd, K, constraint)
            # Check if the deviation bound is within the safety margin
            if maximum(bounded_runs_iter(a, z_0, n, t, safety_margin=d_max)) <= d_max
                # All constraints with (m, window) where m >= meet are valid
                push!(safe_constraints, constraint)
                break
            end
            meet += 1
        end
    end

    safe_constraints
end

# Helper functions

"""
    _undigit(d[, base=2])

Convert a list of digits to a number. Default base=2.
```julia
julia> _undigit([1, 0, 0])
4
```
"""
function _undigit(d; base=2)
    s = zero(eltype(d))
    mult = one(eltype(d))
    for val in reverse(d)
        s += val * mult
        mult *= base
    end
    return s
end

"""
    _digits_b2r(x[, pad])

Digits **b**ase **2**, **r**everse
A shortcut for `digits(x, base=2, pad=pad) |> reverse`
"""
function _digits_b2r(x::Int, pad::Int=0)
    digits(x, base=2, pad=pad) |> reverse
end

"""
    _state_separation(l, B[, indigits=false])

_state_separation takes a number `l` representing the overall state of
a `_SynthesizedAutomaton` and an array `B` representing the number of bits
for each comprising `Automaton` of the `_SynthesizedAutomaton` to compute a 
list of individual states. For example

```julia
julia> _state_separation(6, [1, 2])
[1, 2]
```
The explanation of the example is as follows
```
6 -> [1, 1, 0]
[1, 1, 0] / [1, 2] -> [[1], [1, 0]]
[[1], [1, 0]] -> [1, 2]
```
"""
function _state_separation(l, B; indigits=false)
    @assert l < 2^sum(B) "l has more bits than the sum of B"
    
    bits = digits(l, base=2, pad=sum(B)) |> reverse
    index_points = cumsum(B)
    bits_separated = [bits[i+1:j] for (i, j) in zip([0; index_points[1:end-1]], index_points)]
    if indigits
        state_separated = bits_separated
    else
        state_separated = map(_undigit, bits_separated)
    end
end


"""
Struct for an abstract controller automaton which enforces only the weakly 
hard constraint for that controller. i.e., it is not concerned with the
dynamical matrix of the system.
"""
struct _AbstractAutomaton
    L::Int64			# # of locations. Legal locations are in the range 0:L-1.
    Σ::Int64			# # of actions. Legal actions are in the range 0:Σ-1.
    T::Function			# Transition function. T(l,σ) is a location in 0:L-1.
    l_0::Int64			# Initial location in L.
    Q_f::Function		# Function that returns whether a location is final
end

"""
Struct for a synthesized automaton from multiple abstract automata.
"""
struct _SynthesizedAutomaton
    N::Int64			# # of comprising automata
    B::Vector{Int64}	# List of bits for all comprising controllers
    L::Int64			# Locations. Legal locations are in the range 0:L-1.
    Σ::Vector{Int64}	# List of actions.
    T::Function			# Transition function.  T(l,σ) is a location in 0:L-1.
    l_0::Int64			# Initial location in L.
    Q_f::Function		# Function that returns whether a location is final.
end

"""
Build an abstract controller automaton for the given weakly hard constraint
"""
function _controller_automaton(c::RealTimeScheduling.WeaklyHardConstraint)
    # Define the automaton's structure. State l=0 is the dead state (trapping)
    if c.meet == 0
        # No requirement. Always valid.
        L = 2
        T = (l, σ) -> 1
    elseif c.meet == c.window
        # No misses allowed.
        L = 2
        T = (l, σ) -> l & σ
    else
        # Build full automaton. Dead state is l = 1
        L = 2^c.window
        T = function (l, σ)
            l_new = (l << 1) & (L - 1) | σ
            l == 0 || count_ones(l_new) < c.meet ? 0 : l_new
        end
    end
    
    # Put it all together. Starting position is L-1 since there is no miss from the beginning
    _AbstractAutomaton(L, 2, T, L-1, !iszero)
end

"""
Build a scheduler automaton from a given list of controller automata
"""
function _scheduler_automaton(controllers::Vector{_AbstractAutomaton}; slotsize::Int64=1)
    # Converting tuple to array with collect()
    N = length(controllers)
    B = map(a -> ceil(Int64, log2(a.L)), controllers) |> collect
    L = 2^sum(B)

    all_actions = 0:2^length(controllers)-1
    Σ = filter(σ -> count_ones(σ) <= slotsize, all_actions)

    function T(l::Int, σ::Int)
        @assert l < L "Illegal location"
        @assert σ in Σ "Illegal action"
        states = _state_separation(l, B)
        actions = _digits_b2r(σ, N)
        # @info σ, actions
        new_locations = map((controller, l, σ) -> controller.T(l, σ), controllers, states, actions)
        # @info states, actions, new_locations
        new_location_bits = map((x, y) -> _digits_b2r(x, y), new_locations, B)
        # @info new_location_bits
        result = Iterators.flatten(new_location_bits) |> collect |> _undigit
        # @info _digits_b2r(result)
        result
    end

    function Q_f(l::Int)
        states = _state_separation(l, B);
        all(map((c, l) -> c.Q_f(l), controllers, states))
    end

    _SynthesizedAutomaton(N, B, L, Σ, T, L-1, Q_f)
end

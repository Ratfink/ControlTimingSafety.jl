"""
    schedule_xghtc(constraints; slotsize=1, H=100)

Generate a schedule for a set of weakly hard constraints. The schedule assumes that all tasks are synchronous and have equal periods. At most `slotsize` tasks may be scheduled in a single period. The schedule has total length `H` periods.
    
Shengjie Xu, Bineet Ghosh, Clara Hobbs, P.S. Thiagarajan, and Samarjit Chakraborty, "Safety-Aware Flexible Schedule Synthesis for Cyber-Physical Systems using Weakly-Hard Constraints." ASP-DAC 2023.
"""
function schedule_xghtc(constraints::Vector{RealTimeScheduling.WeaklyHardConstraint}; slotsize::Integer=1, H::Integer=100)

end

"""
    function synthesize_constraints(sysd, K, z_0, d_max, maxwindow, n, t)

Find all `MeetAny` weakly hard constraints with window size at most `maxwindow` that guarantees the deviation upper bound is at most `d_max`. The system is specified by [`Automaton`](@ref) `a` and initial state is `z_0`. `n` and `t` are as in [`bounded_runs_iter`](@ref).
"""
function synthesize_constraints(sysd::AbstractStateSpace{<:Discrete},
    K::AbstractMatrix{Float64}, z_0::AbstractVecOrMat, d_max::Float64,
    maxwindow::Integer, n::Integer, t::Integer)

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
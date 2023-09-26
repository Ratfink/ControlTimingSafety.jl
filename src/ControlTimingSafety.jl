module ControlTimingSafety

using LinearAlgebra
using Random
using ControlSystemsBase
using Distances

export Automaton, nlocations, nactions
export hold_skip_next, zero_skip_next, hold_kill, zero_kill
export strat_map, strat_names
export evol_final, evol, evol_final!, evol!, augment
include("automata.jl")

export corners_from_bounds, merge_bounds, merge_bounds!
export bounded_runs, bounded_runs_iter, deviation
include("safety.jl")

export maximum_deviation_random, estimate_deviation
include("probablesafety.jl")

export schedule_xghtc, synthesize_constraints
export estimate_constraints, verify_schedule, schedule_to_sequence
export devub, devest
include("schedule_synthesis.jl")

end

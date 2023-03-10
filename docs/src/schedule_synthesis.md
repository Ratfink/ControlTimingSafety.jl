# Schedule Synthesis

The main functionality of this section is divided into two functions
1. [`synthesize_constraints`](@ref): given the dynamics and safety margin of **one**
   control task, find all weakly-hard constraints (up to a maximum window size `maxwindow`)
   with which the system behaves safely,
2. [`estimate_constraints`](@ref): similar to [`synthesize_constraints`](@ref), but uses 
   the statistical method [`estimate_deviation`](@ref) to find the deviation upper bound of
   each weakly-hard constraints.
3. [`schedule_xghtc`](@ref): given the lists of safe weakly-hard constraints for **all** 
   control tasks, synthesizes a schedule so that they all behave safely (if such schedules 
   exist).

```@docs
synthesize_constraints
schedule_xghtc
```

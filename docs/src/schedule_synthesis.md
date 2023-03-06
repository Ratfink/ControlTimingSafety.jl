# Schedule Synthesis

The main functionality of this section is divided into two functions
1. [`synthesize_constraints`](@ref): given the dynamics and safety margin of **one**
   control task, find all weakly-hard constraints (up to a maximum window size `maxwindow`)
   with which the system behaves safely, and
2. [`schedule_xghtc`](@ref): given the lists of safe weakly-hard constraints for **all** 
   control tasks, synthesizes a schedule so that they all behave safely (if such schedules 
   exist).

```@docs
synthesize_constraints
synthesize_constraints_deviation
schedule_xghtc
schedule_optimization
```

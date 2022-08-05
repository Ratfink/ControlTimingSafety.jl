### A Pluto.jl notebook ###
# v0.19.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ a0d254fa-beba-424a-b77a-13c9bfce2db9
begin
	using LinearAlgebra
	using ControlSystems
	using Random
	using Plots
	gr()
	using PlutoUI
	using LaTeXStrings
	using QuadGK
	using Distances
	using SparseArrays
	using OffsetArrays
	using Distributions

	using CSV
	using Tables
	using Printf
	
	md"""Import needed packages."""
end

# ╔═╡ afeab312-2d5c-11ec-36d7-8fa4f9a66fa8
md"""
# Control with Deadline Misses

In this notebook, we examine the effects of deadline misses on control systems.  We consider a linear time-invariant plant, with a feedback controller implemented on a processor.  This controller runs periodically, and may experience deadline misses.  Different strategies may be employed when deadlines are missed.  We describe these strategies as finite-state transducer automata, mapping strings of deadline hits and misses to strings of dynamics matrices.
$(TableOfContents())

## Automata Basics

To begin our implementation, we first define an automaton struct.
"""

# ╔═╡ e0bea5c9-c711-4703-ac3a-4df5e0042e6d
struct Automaton
	# L: locations.  Legal locations are in the range 1:L.
	L::Int64
	# A: actions (input alphabet).  Legal actions are in the range 1:A.
	A::Int64
	# Φ: dynamics matrices (output alphabet).  Array of square matrices of equal size.
	Φ::Vector{AbstractMatrix{Float64}}
	# T: transition function.  T[l,a] is a location in 1:L, or missing.
	T::Matrix{Union{Missing, Int64}}
	# μ: output function.  μ[l,a] is a character from Φ.
	μ::Matrix{Union{Missing, Int64}}
	# l_int: initial location in L.
	l_int::Int64
end

# ╔═╡ 89dbc50b-2ce7-46c7-bdec-73166abef4ae
md"""
For convenience later on, we'll also want this little function to make a new Automaton that's an exact copy of another, but with a different initial location.
"""

# ╔═╡ ae811d13-c813-40f3-ba7a-878568291252
Automaton_lint(other::Automaton, l_int::Int64) = Automaton(other.L, other.A, other.Φ, other.T, other.μ, l_int)

# ╔═╡ 7431a2e3-5b31-43b0-a40b-ab0a977e919d
md"""
Now we define a function that runs an automaton, given an initial augmented state and string of actions.  It returns a matrix whose columns are the state vector at each time step, and the final location in the automaton.  This information allows us to continue the evolution in a second call to Evol, which we will use later.
"""

# ╔═╡ dac9640f-eac8-401d-aada-62c864a758ef
function Evol(z_0, automaton, input)
    t_max = size(input, 1)
    z = zeros(size(z_0, 1), t_max + 1)
    z[:,1] = z_0
	l = automaton.l_int
	# For each time step
    for t = 1:t_max
		# Get the dynamics matrix
		μ = automaton.μ[l, input[t]]
		# If we hit a missing transition, return the states that we reached,
		# and a missing final location to signal the problem to the caller.
		if ismissing(μ)
			return z[:,1:t]', missing
		end
		# Apply the dynamics
        #z[:,t+1] = automaton.Φ[:,:,μ] * z[:,t]
		z[:,t+1] = automaton.Φ[μ] * z[:,t]
		# Transition to the new location
		l = automaton.T[l, input[t]]
    end
    z', l
end

# ╔═╡ 264e5c66-5a92-4a4c-8c0c-056e077856ab
md"""
One more utility function next, to pad a state $x$ (or array of states) into an augmented state $z$ used by an automaton.  This will be convenient when we use the automata later on.
"""

# ╔═╡ e86df2ab-534b-471a-90b7-43c7f980abff
Augment(x, automaton) = [x; zeros(size(automaton.Φ[1],1) - size(x,1), size(x,2))]

# ╔═╡ f5a869b9-b2a8-41a7-a503-8523ae5842c6
md"""
Now we are ready to run some simulations.  But we can't just yet, because we don't have any automata to simulate!  We can fix that by implementing a function to construct a Hold&Skip-Next automaton, given a discrete state-space model and a feedback gain matrix K.

Under this strategy, on a deadline miss, the old control input is *held*, while the job that missed is allowed to keep running and the *next* jobs are *skipped* until it completes.

!!! note

    The semantics of Hold&Skip-Next are following [^maggiostability], but the implementation is different, allowing for any number of sequential deadline misses.  This makes the augmented matrices much smaller, resulting in much faster computation.
"""

# ╔═╡ 8c489d40-d44d-48d0-845a-8afbd5161c0b
function HoldAndSkipNext(sysd, K, miss=nothing)
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
		L = 2
		T = [2 1;
			 2 1]
		μ = [3 1;
			 4 2]
	elseif miss == 0
		L = 1
		T = [missing 1]
		μ = [missing 1]
	else
		L = miss + 1
		T = [[collect(2:L); missing] ones(Int64, (L, 1))]
		μ = [3 1;
			 repeat([4 2], L-2, 1);
			 missing 2]
	end
	
	# Put it all together, with the Φ matrices
	Automaton(L, 2, [
			# Φ_HH
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

# ╔═╡ 1d610dac-8989-4bba-99fd-f248fcf10e8c
md"""
We can create automata for numerous other strategies as well.  Here we implement Zero&Skip-Next, which is similar to Hold&Skip-Next, except that when a deadline is missed, the control input applied in the next period is *zero*, rather than holding the current control input.
"""

# ╔═╡ 3d095e19-4f82-4c68-a104-2b7d94467849
function ZeroAndSkipNext(sysd, K, miss=nothing)
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
		L = 2
		T = [2 1;
			 2 1]
		μ = [3 1;
			 4 2]
	elseif miss == 0
		L = 1
		T = [missing 1]
		μ = [missing 1]
	else
		L = miss + 1
		T = [[collect(2:L); missing] ones(Int64, (L, 1))]
		μ = [3 1;
			 repeat([4 2], L-2, 1);
			 missing 2]
	end
	
	# Put it all together, with the Φ matrices
	Automaton(L, 2, [
		    # Φ_HH
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

# ╔═╡ a92710a5-9606-4ba6-88a4-962bca2a2c3e
md"""
There are other approaches that can be taken as well.  Instead of Skip-Next, we can follow the *Kill* strategy, where a job that misses its deadline is killed and the next job is released as usual.
"""

# ╔═╡ fdd6db9c-2d42-4073-b7ca-22d6396b86ab
function HoldAndKill(sysd, K, miss=nothing)
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
		L = 1
		T = [1 1]
		μ = [2 1]
	elseif miss == 0
		L = 1
		T = [missing 1]
		μ = [missing 1]
	else
		L = miss + 1
		T = [[collect(2:L); missing] ones(Int64, (L, 1))]
		μ = [2 1;
			 repeat([2 1], L-2, 1);
			 missing 1]
	end
	
	# Put it all together, with the Φ matrices
	Automaton(L, 2, [
		    # Φ_H
			[sysd.A  sysd.B;
			 K_x  K_u],
			# Φ_M
			[sysd.A  sysd.B;
	         zeros(r, p)  I(r)]],
		T, μ, 1)
end

# ╔═╡ fcb00b19-7eb1-40ca-b7a6-8fbf61397172
function ZeroAndKill(sysd, K, miss=nothing)
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
		L = 1
		T = [1 1]
		μ = [2 1]
	elseif miss == 0
		L = 1
		T = [missing 1]
		μ = [missing 1]
	else
		L = miss + 1
		T = [[collect(2:L); missing] ones(Int64, (L, 1))]
		μ = [2 1;
			 repeat([2 1], L-2, 1);
			 missing 1]
	end
	
	# Put it all together, with the Φ matrices
	Automaton(L, 2, [
		    # Φ_H
			[sysd.A  sysd.B;
			 K_x  K_u],
			# Φ_M
			[sysd.A  sysd.B;
	         zeros(r, p + r)]],
		T, μ, 1)
end

# ╔═╡ 2481a5d6-71dd-4530-80c0-a6b4eae4a13b
function HoldAndContinue(sysd, K, miss=nothing)
	# FIXME: miss must not be used :v:v
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
		L = 1
		T = [1 1 1]
		μ = [1 2 3]
	elseif miss == 0
		L = 1
		T = [missing 1]
		μ = [missing 1]
	else
		L = miss + 1
		T = [[collect(2:L); missing] ones(Int64, (L, 1))]
		μ = [2 1;
			 repeat([2 1], L-2, 1);
			 missing 1]
	end
	
	# Put it all together, with the Φ matrices
	Automaton(L, 3, [
		    # Φ_1
			[sysd.A  zeros(p, p)  sysd.B;
			 I(p)  zeros(p, p+r);
			 K_x  zeros(r, p)  K_u],
			# Φ_2
			[sysd.A  zeros(p, p)  sysd.B;
			 I(p)  zeros(p, p+r);
			 zeros(r, p)  K_x  K_u],
			# Φ_H
			[sysd.A  zeros(p, p)  sysd.B;
			 I(p)  zeros(p, p+r);
			 zeros(r, 2p)  I(r)]],
		T, μ, 1)
end

# ╔═╡ 9cddd940-d4e3-44e1-ac01-de314f9e8517
function ZeroAndContinue(sysd, K, miss=nothing)
	# FIXME: miss must not be used :v:v
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
		L = 1
		T = [1 1 1]
		μ = [1 2 3]
	elseif miss == 0
		L = 1
		T = [missing 1]
		μ = [missing 1]
	else
		L = miss + 1
		T = [[collect(2:L); missing] ones(Int64, (L, 1))]
		μ = [2 1;
			 repeat([2 1], L-2, 1);
			 missing 1]
	end
	
	# Put it all together, with the Φ matrices
	Automaton(L, 3, [
		    # Φ_1
			[sysd.A  zeros(p, p)  sysd.B;
			 I(p)  zeros(p, p+r);
			 K_x  zeros(r, p)  K_u],
			# Φ_2
			[sysd.A  zeros(p, p)  sysd.B;
			 I(p)  zeros(p, p+r);
			 zeros(r, p)  K_x  K_u],
			# Φ_H
			[sysd.A  zeros(p, p)  sysd.B;
			 I(p)  zeros(p, p+r);
			 zeros(r, 2p+r)]],
		T, μ, 1)
end

# ╔═╡ 1317e2b7-2b8a-4bb8-aee8-bcd6f8d38a03
md"""
Lastly, we define a mapping from human-readable strategy names to their constructor functions.
"""

# ╔═╡ bddaa61b-851f-4759-889a-b71adff1ca57
begin
	strat_map = Dict(
		"Hold&Skip-Next" => HoldAndSkipNext,
		"Zero&Skip-Next" => ZeroAndSkipNext,
		"Hold&Kill" => HoldAndKill,
		"Zero&Kill" => ZeroAndKill,
		"Hold&Continue" => HoldAndContinue,
		"Zero&Continue" => ZeroAndContinue
	)
	strat_names = sort([keys(strat_map)...])
end

# ╔═╡ 9772aa30-71d2-4083-a3ad-d842b0fd2b13
md"""
## Demonstration of Automata

Now, to create an instance of a Hold&Skip-Next automaton, we create a two-dimensional LTI system: an RC network, from [^gabelroberts] (pp. 160–161).
"""

# ╔═╡ 03aa61fd-67ca-45f5-9eef-7464901598cc
sys = let
	r_1 = 100000
	r_2 = 500000
	r_3 = 200000
	c_1 = 0.000002
	c_2 = 0.000010
	A = [-1/c_1 * (1/r_1 + 1/r_2)  1/(r_2*c_1)
	     1/(r_2*c_2)               -1/c_2 * (1/r_2 + 1/r_3)]
	#B = [1/(r_1*c_1)  0
	#     0            1/(r_3*c_2)]
	B = [1/(r_1*c_1)
	     1/(r_3*c_2)]
	C = [1  -1]
	D = 0
	#A = [10 0; -2 -1]
	#B = [5 1; 4 10]
	#C = [0 0]
	#D = [0 0]=#
#=	J = 0.1;
b = 0.1;
K = 0.1;
R = 0.1;
L = -0.5;
A = [-b/J   K/J
    -K/L   -R/L];
B = [0
    1/L];
C = [1   0];
D = 0;=#
	#=R = 0.025
	ωel = 2000*pi
	Ld = 0.0001
	Lq = 0.00012
	A = [-R/Ld  Lq * ωel / Ld;  -Ld * ωel / Lq   -R / Lq]
	B = [1 / Ld   0;  0   1/Lq]
	C = [0 0]
	D = [0 0]=#
	#=A = [-0.313 56.7 0; -0.0139 -0.426 0; 0 56.7 0]
	B = [0.232; 0.0203; 0]
	C = [0 0 1]
	D = [0]=#
	A = [-0.313 56.7 0; -0.0139 -0.426 0; 0 56.7 0];
	B = [0.232; 0.0203; 0];
	C = [0 0 1];
	D = [0];
	ss(A, B, C, D)
end

# ╔═╡ 0cd0ceff-ebc5-4069-9bb9-14fcc667a1ef
md"""
We then discretize the model with a constant period $h$ (in seconds), and no sensor-to-actuator delay.

|     |                                                        |
|----:|:-------------------------------------------------------|
| $h$ | $(@bind ctrl_delay NumberField(0:0.05:1, default=0.1)) |
"""

# ╔═╡ acc7b261-8713-4565-9589-f59af764175c
#sysd = c2d(sys, ctrl_delay)
#sysd = ss([0.996 0.075; -0.052 0.996], [0.100 0.003; -0.003 0.083], [0 0], [0 0], 0.00001)
sysd = ss([1 0.13; 0 1], [0.02559055118110236; 0.39370078740157477], [1 0], [0], 0.020)

# ╔═╡ e07b2e53-e0eb-4e00-a25f-869d98d04cbe
md"""
Next, we discretize the model again with the same period, but now with a sensor-to-actuator delay of one period, which is assumed throughout the remainder of this work.  This allows a software task that implements the controller to complete in any amount of time less than $h$, and for the actuator to apply the new control input at the next sampling instant.  This is referred to in the literature as the *logical execution time* (LET) paradigm.

Note that this model uses an augmented, three-dimensional state space.  The code for constructing this augmented model is listed below.
"""

# ╔═╡ 1faf171b-9a38-42dc-867f-4fd576228bde
function c2da(sysc::AbstractStateSpace{<:ControlSystems.Continuous}, Ts::Real, d::Real)
    f_Φ = s -> ℯ^(sysc.A * s)
    Φ = f_Φ(Ts)
    Γ_0 = quadgk(f_Φ, 0, Ts-d)[1] * sysc.B
    Γ_1 = quadgk(f_Φ, Ts-d, Ts)[1] * sysc.B
    Φ_a = [Φ  Γ_1;
		   zeros(size(sysc.B, 2), size(sysc.A, 2)+size(sysc.B, 2))]
    Γ_a = [Γ_0;
		   I]
    C_a = [sysc.C  zeros(size(sysc.C,1), size(sysc.B,2))]
    D_a = zeros(size(C_a,1), size(Γ_a,2))
    ss(Φ_a, Γ_a, C_a, D_a, Ts)
end

# ╔═╡ e1c05296-ddb4-4da0-b972-d952ab3a5ef2
sysd_delay = c2da(sys, ctrl_delay, ctrl_delay)

# ╔═╡ fc69b358-f22c-450f-9ab0-5ed9ed564207
md"""
We now create an optimal controller for the delayed system using LQR.

|              |                                                        |
|-------------:|:-------------------------------------------------------|
|   State cost | $(@bind state_cost NumberField(0:0.1:10, default=2))   |
| Control cost | $(@bind control_cost NumberField(0:0.1:10, default=1)) |
"""

# ╔═╡ 72650fe6-cd13-46f1-9ce7-10831a27c396
K = let
	#Q = [state_cost 0 0;
	#	 0 state_cost 0;
	#	 0 0 state_cost]
	#R = [control_cost 0; 0 0][1:1,1:1]  # Hack to make a 1x1 matrix.  [1] makes a vector, which has no isposdef method
	#lqr(sysd_delay, Q, R)
	#lqr(sys, [1 0; 0 1], [1 0; 0 1])
	#=-[0 0 5 0 -4 0;
	 0 0 1 7 -3 7;
	 0 0 0 0 0 0;
	 0 0 0 0 0 0]=#
	#dlqr(sysd, Array(I(2)), Array(I(2)))
	#Q = [0 0 0 0; 0 0 0 0; 0 0 state_cost 0; 0 0 0 0]
	#R = control_cost;
	#dlqr(sysd_delay, Q, R)
	[0.293511 0.440267]
end

# ╔═╡ 06c7c6d1-719e-4c9a-bf28-939288c335e9
md"""
And finally, we simulate the delayed system under a user-selectable strategy, both with deadline misses and without.
The nominal behavior, where no deadlines are missed, is shown in blue; other colors show trajectories with probabilistic deadline misses.  Try changing the miss strategy and simulation parameters below!

|                          |      |
|-------------------------:|:-----|
|                 Strategy | $(@bind sim_strategy Select(strat_names)) |
| Orange $P(\mathrm{hit})$ | $(@bind hit_prob_1 Slider(0:0.05:1, default=0.90, show_value=true)) |
|    Red $P(\mathrm{hit})$ | $(@bind hit_prob_2 Slider(0:0.05:1, default=0.75, show_value=true)) |
| Purple $P(\mathrm{hit})$ | $(@bind hit_prob_3 Slider(0:0.05:1, default=0.35, show_value=true)) |
|          Simulation time | $(@bind sim_time Slider(5:5:500, default=500, show_value=true)) |
|                 $x_1[0]$ | $(@bind x_0_1 NumberField(-10:0.1:10, default=1)) |
|                 $x_2[0]$ | $(@bind x_0_2 NumberField(-10:0.1:10, default=1)) |

|                                                |
|:----------------------------------------------:|
| $(@bind recompute_plot Button("🎲 Randomize")) |
"""

# ╔═╡ 9ed43a4f-f1c3-4299-81ff-2345f4b9338a
md"""
Through this example, the effects of deadline misses can be observed.  When the probability of a hit is very low (say, 0.05), there are commonly very long strings of deadline misses, followed by a hit.  Once that hit occurs, the control input that is applied was computed with stale data, and may be far from the best control input to apply at present.  This leads to behavior that is hard to predict.  Our aim is to determine how far any trajectory can get from the nominal behavior.
"""

# ╔═╡ bb865586-2365-410a-be80-04b65c935381
md"""
### Example 2: Random sequences of actions

We next show a second example in which the same system is run for a large number of hit/miss sequences, selected at random according to a uniform distribution with the constraint of no more than $n$ consecutive misses.
"""

# ╔═╡ d9d3bb42-43e0-4b68-8135-6447294de688
md"""
|                 |                                                              |
|----------------:|:-------------------------------------------------------------|
|        Strategy | $(@bind strategy_2 Select(strat_names))                      |
|      Max misses | $(@bind maxmiss_2 Slider(-1:16, default=5, show_value=true)) |
| Simulation time | $(@bind H_2 Slider(5:5:500, default=500, show_value=true))   |
|        $x_1[0]$ | $(@bind x₀_1_2 NumberField(-10:0.1:10, default=1))          |
|        $x_2[0]$ | $(@bind x₀_2_2 NumberField(-10:0.1:10, default=1))          |

|                                                 |
|:-----------------------------------------------:|
| $(@bind recompute_plot_2 Button("🎲 Randomize")) |
"""

# ╔═╡ 9179f1a5-26df-4e20-8e88-cbf6df6c24e6
md"""
In order for this random sampling to work, we need a way to perform this uniform sampling.  Fortunately, Julia provides a nice API for sampling random values from different types of data and distributions.  We next create a `Sampler` type for sampling strings from an `Automaton` object.
"""

# ╔═╡ f7161723-6e58-4ddd-ac74-f2b1488e6d2e
struct SamplerAutomatonInput <: Random.Sampler{Int64}
	l::Matrix{BigInt}
	function SamplerAutomatonInput(k_miss, H)
		l = zeros(BigInt, k_miss+2, H+1)
		l[1:k_miss+1, 1] .= 1
		for i = 2:H+1
			for q = 1:k_miss+2
				if q == k_miss+2
					l[q, i] = 2 * l[q, i-1]
				else
					l[q, i] = l[q+1, i-1] + l[1, i-1]
				end
			end
		end
		new(l)
	end
end

# ╔═╡ ffdc8ef4-2eea-4ed3-a57d-ea0f14219f7b
md"""
We next define a method of the function `Random.rand` that allows us to sample from this distribution.
"""

# ╔═╡ 14847de2-a3a8-4206-8bf7-64a9896fa750
function Random.rand(rng::AbstractRNG, sp::SamplerAutomatonInput)
	k_miss = size(sp.l, 1) - 2
	H = size(sp.l, 2) - 1
	q = 0
	ret = zeros(Int64, H)
	for i = 1:H
		d = BigFloat(sp.l[q + 1, H - i + 2])
		if q == k_miss + 1
			prob_one = sp.l[q + 1, H - i + 1] / d
		else
			prob_one = sp.l[1, H - i + 1] / d
		end
		rand_bit = Random.rand(Binomial(1, Float64(prob_one)))
		ret[i] = rand_bit
		if rand_bit == 1
			if q != k_miss + 1
				q = 0
			end
		elseif q != k_miss + 1
			q += 1
		end
	end
	ret .+ 1
end

# ╔═╡ ff4a5ea5-96f5-4054-ab8b-681aef3ad017
let
	# Recompute when the button is clicked
	recompute_plot
	
	hsn = strat_map[sim_strategy](sysd, K)
	
	hit_prob = [1, hit_prob_1, hit_prob_2, hit_prob_3]
	line_colors = [:blue, :orange, :red, :purple]
	z_0 = Augment([x_0_1, x_0_2], hsn)
	
	z = [0]
	plot(title="Trajectories for different hit probabilities", xlabel=L"x_1", ylabel=L"x_2", legend=:bottomright)
	for i = size(hit_prob, 1):-1:1
	    seq = (rand(sim_time) .<= hit_prob[i]) .+ 1
		z, _ = Evol(z_0, hsn, seq)
	    plot!(z[:,1], z[:,2], label=L"P(\mathrm{hit}) = %$(hit_prob[i])", linecolor=line_colors[i])
		#plot!(0:sim_time, z[:,1], label=L"P(\mathrm{hit}) = %$(hit_prob[i])", linecolor=line_colors[i])
	end
	plot!()
end

# ╔═╡ 5406484b-7c43-41df-8f63-3e1d16fde3aa
# ╠═╡ disabled = true
#=╠═╡
let
	# Recompute when the button is clicked
	recompute_plot_2

	T = strat_map[strategy_2](sysd, K)
	z₀ = Augment([x₀_1_2, x₀_2_2], T)

	plot(xlabel=L"x_1", ylabel=L"x_2", zlabel=L"t", legend=:topleft, aspectratio=1, format=:png, thickness_scaling=2, size=(1360,906), camera=(70,15))
	
	# Calculate nominal trajectory
	seq = ones(Int64, H_2) .+ 1
	nominal, _ = Evol(z₀, T, seq)

	junk = 0

	# Compute random trajectories
	sp = SamplerAutomatonInput(maxmiss_2, H_2)
	n_samples = 100
	trj = Array{Float64}(undef, n_samples, H_2+1, size(T.Φ[1], 1))
	dist = Array{Float64}(undef, n_samples, H_2+1)
	max_possible = (H_2 / (maxmiss_2 + 1)) * maxmiss_2 + H_2 % (maxmiss_2 + 1)
	for i = axes(trj, 1)
		accepted = false
		while !accepted
			seq = rand(sp)
			if rand(Float64) > sum(seq .== 1) / (max_possible * .5)
				accepted = true
			else
				junk += 1
			end
		end
		trj[i,:,:], _ = Evol(z₀, T, seq)
		dist[i,:] = colwise(Euclidean(), trj[i,:,1:2]', nominal[:,1:2]')
	end
	dev = maximum(dist, dims=2)
	@info junk

	# Plot safety pipe
	θ = LinRange(0, 2π, 360)
	pipe_radius = sort(dev, dims=1)[end-1] - 0.01
	@info pipe_radius
	for i in axes(nominal, 1)
		plot!(cos.(θ) * pipe_radius .+ nominal[i,1], sin.(θ) * pipe_radius .+ nominal[i,2], repeat([i,], size(θ,1)), label=(i == 1) ? "Safety pipe" : "", seriestype=[:shape,], c=:lightblue, linecolor=:lightblue)
	end

	# Plot random trajectories
	label_good = false
	# Plot the good trajectories first
	for i = axes(trj, 1)
		if dev[i] < pipe_radius
			plot!(trj[i,:,1], trj[i,:,2], 1:H_2+1, label=(!label_good) ? "Random" : "", color=:green, opacity=10/n_samples)
			label_good = true
		end
	end
	# Plot the bad trajectories on top of the good ones
	label_bad = false
	nom_crosses = zeros(Float64, H_2+1)
	for i = axes(trj, 1)
		if dev[i] >= pipe_radius
			t_viol = argmax(dist[i,:], dims=1)
			#label_bad || plot!(cos.(θ) * pipe_radius .+ nominal[t_viol,1], sin.(θ) * pipe_radius .+ nominal[t_viol,2], repeat([t_viol,], size(θ,1)), label="Radius violated", style=:dot, linecolor=:gray)
			crosses = [(d >= pipe_radius) ? 1 : 0 for d in dist[i,:]]
			plot!(trj[i,:,1], trj[i,:,2], 1:H_2+1, label=(!label_bad) ? "Violation" : "", color=:red, markershape=:x, markeralpha=crosses)
			nom_crosses .+= crosses
			label_bad = true
		end
	end

	# Finally, plot the nominal trajectory on top of everything else
	plot!(nominal[:,1], nominal[:,2], 1:H_2+1, label="Nominal", color=:black, linewidth=3, marker=:x, markeralpha=nom_crosses)
end
  ╠═╡ =#

# ╔═╡ 8f202fe5-008b-4105-8322-8f595fb6b545
md"""
## Bounded Tree Reachable Set

In this section, we use a *bounded tree method* to compute a reachable set for the system under any sequence of deadline hits and misses.  In this method, we start with an axis-aligned box initial set, and simulate its evolution for $n$ steps, for each of the $|\mathcal{A}|^n$ possible sequences of scheduler actions.  We then take an axis-aligned bounding box of the resulting sets, and repeat this process.

This algorithm more generally supports running any automaton, including those that implement different deadline miss strategies, and weakly-hard constraints.  This gives us more power to model scheduler behavior.

We begin by defining a function to run one iteration of this process.  For the iteration to be sound, each iteration must return one bounding box for each final location in the automaton.  Then, a separate tree will be run for each final location output by the previous round, so this function must also take as a parameter the initial location to use.
"""

# ╔═╡ f61f27cb-48d3-4f95-8547-5b0a3a2b6394
function BoundedTree(automaton, bounds, n)
	q = size(bounds, 1)
	
	# Bounding boxes for each final location, time step
	ret = Array{Float64}(undef, automaton.L, n+1, q, 2) * NaN
	
	# For each corner
	for corner in Base.product(eachrow(bounds)...)
		# For each hit/miss sequence
		for seq in Base.product(eachcol(repeat(1:automaton.A, outer=[1,n]))...)
			# Compute the evolution for this corner and sequence
			z, l_fin = Evol([c for c in corner], automaton, seq)
			if ismissing(l_fin)
				continue
			end
			# Calculate min and max for this final location at each time step
			ret[l_fin,:,:,1] = minimum(x->(isnan(x) || isinf(x)) ? Inf : x, cat(z, ret[l_fin,:,:,1], dims=3), dims=3)
			ret[l_fin,:,:,2] = maximum(x->(isnan(x) || isinf(x)) ? -Inf : x, cat(z, ret[l_fin,:,:,2], dims=3), dims=3)
		end
	end
	
	ret
end

# ╔═╡ 5ae87150-d2cb-4e84-b7c9-bb255b4f36e0
md"""
While this function gets the job done, it's far from the most efficient way to implement the bounded tree search.  In particular, it computes the entire trajectory for each sequence of actions, when most of this work could be shared between sequences.  This means the complexity of this first implementation is actually $O(|\mathcal{A}|^n \cdot n)$, which is asymptotically faster-growing than the $O(|\mathcal{A}|^n)$ growth of the number of sequences.

Instead, we implement a second function performing the same job, but instead doing a preorder traversal of the automaton.  This removes the redundant computation, giving an asymptotic speedup by bringing us down to the expected time complexity.
"""

# ╔═╡ 97c1f398-a62f-401b-802d-283c8dc64307
function BoundedTreeFast(automaton, bounds, n)
	# Stack
	z = Matrix{Float64}(undef, n+1, size(bounds, 1))
	loc = Vector{Int64}(undef, n+1)
	act = Vector{Int64}(undef, n+1)
	
	# Bounding boxes for each final location, time step
	ret = Array{Float64}(undef, automaton.L, n+1, size(bounds, 1), 2) * NaN

	# For each corner
	for corner in Base.product(eachrow(bounds)...)
		# Create the stack frame for time 0
		z[1,:] = [c for c in corner]
		loc[1] = automaton.l_int
		act[1] = 1
		# Initialize the stack pointer
		sp = 1
		# While we haven't popped all the way out
		while sp >= 1
			# If we've reached a leaf
			if sp == n+1
				# Calculate min and max for this final location at each time step
				ret[loc[sp],:,:,1] = minimum(x->(isnan(x) || isinf(x)) ? Inf : x, cat(z, ret[loc[sp],:,:,1], dims=3), dims=3)
				ret[loc[sp],:,:,2] = maximum(x->(isnan(x) || isinf(x)) ? -Inf : x, cat(z, ret[loc[sp],:,:,2], dims=3), dims=3)
				sp -= 1
			# If we're out of actions from this step
			elseif act[sp] > automaton.A
				sp -= 1
			# If the transition is missing
			elseif ismissing(automaton.T[loc[sp], act[sp]])
				# Try the next transition
				act[sp] = act[sp] + 1
			# If the transition is present
			else
				z[sp+1,:] = automaton.Φ[automaton.μ[loc[sp], act[sp]]] * z[sp,:]
				loc[sp+1] = automaton.T[loc[sp], act[sp]]
				act[sp+1] = 1
				act[sp] = act[sp] + 1
				sp = sp + 1
			end
		end
	end
	ret
end

# ╔═╡ 0eb09c86-d48e-4a1e-b7ea-5f1c96df2667
begin
	#bounds = [10 10; 10 10]
	#bounds = [-7 -7; -4 -4]
	#bounds = [-6 -6; -8 -8]
	#bounds = [9 9; 9 9]
	#bounds = [-5 -5; 2 2]
	#bounds = let θ = 7π/16; [cos(θ) cos(θ); sin(θ) sin(θ)] end
	bounds = [2 2; 2 2]
	#bounds = [0.1 0.1; 0 0]
	
	md"""
Using this function, we can compute a bounding box of the reachable set over $n$ time steps.  For example, using the system and controller from the previous section:

|                     |                                                                      |
|--------------------:|:---------------------------------------------------------------------|
|            Strategy | $(@bind sim_strategy_one_bounded_tree Select(strat_names)) |
|                 $n$ | $(@bind n_one_bounded_tree Slider(1:16, default=8,  show_value=true)) |
|          Max misses | $(@bind maxmiss_one_bounded_tree Slider(-1:16, default=5, show_value=true)) |
| Sample trajectories | $(@bind sample_traj CheckBox(default=false))                         |
|         Unit circle | $(@bind unit_circ CheckBox(default=true))                         |
"""
end

# ╔═╡ 8e1164d6-9559-4103-b4e6-577166d0ab2b
begin
	points_one_bounded_tree = let
		automaton = strat_map[sim_strategy_one_bounded_tree](sysd, K, maxmiss_one_bounded_tree)
		augbounds = Augment(bounds, automaton)
		BoundedTreeFast(automaton, augbounds, n_one_bounded_tree)
	end
	
	md"""
	The reachable sets are computed in this cell and shown below.  Try enabling the sample trajectories, verifying that they remain inside the boxes for each time step.
	
	!!! warning
	
	    This is a brute-force reachable set for all hit/miss sequences of length $n$, with time complexity $O(|\mathcal{A}|^n)$.  It may take a couple minutes to compute with large values set on the slider above.
	"""
end

# ╔═╡ da2987de-5bd6-46a8-ba21-58f5d2a2f397
md"""
Notice that the time taken to compute this reachable set grows exponentially in $n$, so smaller values will result in shorter execution times in the bounded tree iteration.  We implement this iteration next.
"""

# ╔═╡ a7761898-3e69-4371-9315-8282fe07407d
md"""
The reachable sets computed by the bounded tree method are shown below.  Try adjusting the parameters used for the model: the number of time steps per tree $n$, and the total number of time steps $t$.

|                     |                                                              |
|--------------------:|:-------------------------------------------------------------|
|            Strategy | $(@bind sim_strategy_trees Select(strat_names))              |
|                 $n$ | $(@bind n_trees Slider(1:24, default=5, show_value=true))    |
|          Max misses | $(@bind max_miss_trees Slider(-1:16, default=4, show_value=true)) |
|                 $t$ | $(@bind time_trees Slider(5:5:1000, default=100, show_value=true)) |
| Sample trajectories | $(@bind sample_traj_trees CheckBox(default=false))           |
|         Unit circle | $(@bind unit_circ_trees CheckBox(default=true))              |
"""

# ╔═╡ 883adf52-b18c-4f24-9189-601b5a50d7c2
md"""
## Computing an Upper-Bound on Deviation

Given reachable sets at each time step, it is straightforward to compute a bound to the deviation from a nominal trajectory.  We can simply compute the system evolution from the corners of the initial set using the nominal behavior, then calculate the distance between these corners at each time step and the corners of the corresponding reachable set.  Taking the maximum yields our deviation bound $\mathbf{d}_\mathit{ub}$.
"""

# ╔═╡ 421ecae9-b1b3-4a9e-af01-6cc9befa5884
md"""
Using this function, we can plot the deviation at every time step, as shown below for the bounded tree reachable sets computed in the previous section.  The maximum deviation $\mathbf{d}_\mathit{ub}$ is highlighted on the plot.
"""

# ╔═╡ 7c7ea3d4-fe32-4560-b501-7845265dff52
md"""
## Appendix
"""

# ╔═╡ 0363259a-2096-4825-8cec-bcb906f20a23
md"""
### Bibliography

[^gabelroberts]: Robert Gabel and Richard Roberts.
	"Signals and Linear Systems."
	Second edition.
	Wiley, 1980.
[^maggiostability]: Martina Maggio, Arne Hamann, Eckart Mayer-John, and Dirk Ziegenbein.
	"Control-system stability under consecutive deadline misses constraints."
	ECRTS 2020.
	Online: [https://drops.dagstuhl.de/opus/volltexte/2020/12384/pdf/LIPIcs-ECRTS-2020-21.pdf](https://drops.dagstuhl.de/opus/volltexte/2020/12384/pdf/LIPIcs-ECRTS-2020-21.pdf)
"""

# ╔═╡ 6baa69a2-74af-4c85-8368-49cfa35628b8
md"""
### Miscellaneous Code
"""

# ╔═╡ 6a82ce98-0a17-4f13-8979-94c1eccdd5b4
function merge_bounds(b)
	mins = minimum(x->(isnan(x) || isinf(x)) ? Inf : x, b[:,:,:,1], dims=1)
	maxs = maximum(x->(isnan(x) || isinf(x)) ? -Inf : x, b[:,:,:,2], dims=1)
	cat(dims=4, mins, maxs)[1,:,:,:]
end

# ╔═╡ d7c5ea66-4b3f-4920-9eda-858466a076b3
function BoundedTreeIter(automaton, bounds, n, t)
	q = size(bounds, 1)
	
	# Dimensions: time, augmented state, min/max
	all_bounds = Array{Float64}(undef, n*(t+1)+1, q, 2)
	all_bounds[1,:,:] = bounds
	
	bounds = BoundedTreeFast(automaton, bounds, n)
	all_bounds[1:n+1,:,:] = merge_bounds(bounds)[:,:,1:end]
	
	# Dimensions: initial location, final location, time, augmented state, min/max
	new_bounds = Array{Any}(undef, automaton.L, automaton.L, n+1, q, 2)
	for i in 1:t
		# Simulate each box from previous iteration
		for i in 1:automaton.L
			a = Automaton_lint(automaton, i)
			new_bounds[i,:,:,:,:] = BoundedTreeFast(a, bounds[i,end,:,:], n)
		end
		# Merge resulting boxes from these simulations
		for i in 1:automaton.L
			bounds[i,:,:,:] = merge_bounds(new_bounds[:,i,:,:,:])
		end
		# Save the bounds
		all_bounds[n*i+2:n*(i+1)+1,:,:] = merge_bounds(bounds)[2:end,:,:]
	end
	all_bounds
end

# ╔═╡ b7e416f5-897b-4e45-bd84-a5ef25e7c43d
begin
	t_trees = div(time_trees, n_trees)
	all_bounds = let
		automaton = strat_map[sim_strategy_trees](sysd, K, max_miss_trees)
		augbounds = Augment(bounds, automaton)
		BoundedTreeIter(automaton, augbounds, n_trees, t_trees)
	end
	
	md"""
	The reachable sets are computed in this cell and shown below.
	
	!!! warning
	
		This algorithm has time complexity $O(t |\mathcal{A}|^n)$, so for large values of both parameters, this could take an hour or more.
	"""
end

# ╔═╡ eb79be36-73ba-423d-9d77-a42ab4007246
function corners_from_bounds(bounds; cycle=false, dims=nothing)
	if dims == nothing
		dims = axes(bounds, 1)
	end
	ldims = length(dims)

	corners = cat(reshape([[c...] for c in Base.product(eachrow(bounds[dims,:])...)], 2^ldims)..., dims=2)
	if cycle
		gray(x) = x ⊻ (x >> 1)
		[corners[:,gray.(0:2^ldims-1) .+ 1] corners[:,1]]
	else
		corners
	end
end

# ╔═╡ 28932e44-1762-4ca7-9698-80eeaffa99a1
let
	plot(title="Reachable set for $(n_one_bounded_tree) time steps", xlabel="x_1", ylabel="x_2", legend=:bottomright)
	merged = merge_bounds(points_one_bounded_tree)
	for k = 0:n_one_bounded_tree
		corners = corners_from_bounds(merged[begin+k,:,:], cycle=true)
		plot!(corners[1,:], corners[2,:], label=L"x[%$(k)]")
	end
	
	if sample_traj
		hsn = HoldAndSkipNext(sysd, K)
	    x, _ = Evol(Augment([bounds[1,1], bounds[2,1]], hsn), hsn, ones(Int64, n_one_bounded_tree)*2)
		plot!(x[:,1], x[:,2], label="Sample", linecolor=:blue, marker=:circle)
	    x, _ = Evol(Augment([bounds[1,1], bounds[2,2]], hsn), hsn, ones(Int64, n_one_bounded_tree)*2)
		plot!(x[:,1], x[:,2], label="Sample", linecolor=:blue, marker=:circle)
	    x, _ = Evol(Augment([bounds[1,2], bounds[2,1]], hsn), hsn, ones(Int64, n_one_bounded_tree)*2)
		plot!(x[:,1], x[:,2], label="Sample", linecolor=:blue, marker=:circle)
	    x, _ = Evol(Augment([bounds[1,2], bounds[2,2]], hsn), hsn, ones(Int64, n_one_bounded_tree)*2)
		plot!(x[:,1], x[:,2], label="Sample", linecolor=:blue, marker=:circle)
	end
	
	θ = LinRange(0, 2π, 500)
	if unit_circ
		plot!(cos.(θ), sin.(θ), linecolor=:gray, label="", xlims=[-1, 1], ylims=[-1, 1], aspect_ratio=1)
	end
	plot!()
end

# ╔═╡ 2280b3d6-d3ab-4049-b796-a2a036b55471
let
	plot(title="Reachable set for $(time_trees) time steps, $(n_trees) per tree", xlabel="x_1", ylabel="x_2", legend=:bottomright)
	for t in 1:time_trees+1
		corners = corners_from_bounds(all_bounds[t,:,:], cycle=true)
		plot!(corners[1,:], corners[2,:], label=L"x[%$(t-1)]")
	end
	
	
	if sample_traj_trees
		hsn = HoldAndSkipNext(sysd, K)
	    x, _ = Evol(Augment([bounds[1,1], bounds[2,1]], hsn), hsn, ones(Int64, n_trees*t_trees)*2)
		plot!(x[:,1], x[:,2], label="Sample", linecolor=:blue, marker=:circle)
	    x, _ = Evol(Augment([bounds[1,1], bounds[2,2]], hsn), hsn, ones(Int64, n_trees*t_trees)*2)
		plot!(x[:,1], x[:,2], label="Sample", linecolor=:blue, marker=:circle)
	    x, _ = Evol(Augment([bounds[1,2], bounds[2,1]], hsn), hsn, ones(Int64, n_trees*t_trees)*2)
		plot!(x[:,1], x[:,2], label="Sample", linecolor=:blue, marker=:circle)
	    x, _ = Evol(Augment([bounds[1,2], bounds[2,2]], hsn), hsn, ones(Int64, n_trees*t_trees)*2)
		plot!(x[:,1], x[:,2], label="Sample", linecolor=:blue, marker=:circle)
	end
	
	θ = LinRange(0, 2π, 500)
	if unit_circ_trees
		plot!(cos.(θ), sin.(θ), linecolor=:gray, label="", xlims=[-1, 1], ylims=[-1, 1], aspect_ratio=1)
	end
	plot!(legend=false)
end

# ╔═╡ 728953ad-b271-4ac0-bdc7-e8c71597b092
function deviation(automaton, bounds, reachable; dims=axes(bounds,1), metric=Euclidean(), nominal=repeat([2],size(reachable,1)-1))
	if nominal == nothing
		# Assume a nominal behavior of all hits
		nominal = ones(Int64, size(reachable, 1) - 1) .+ 1
	end
	
	if dims == nothing
		# Assume dimensions 1 and 2 are the plant state
		# TODO: this is mostly just hardcoded for now, ignoring the dims variable.
		#   Partially the fault of corners_from_bounds.
		dims = (1, 2)
	end
	
	# Dimensions: state variables, points, time
	reachable_corners = cat([corners_from_bounds(reachable[t,:,:], dims=dims) for t in axes(reachable, 1)]..., dims=3)
	
	# Dimensions: state variables, points, time
	ev = Array{Float64}(undef, length(dims), 2^size(bounds,1), size(reachable, 1))
	corners = corners_from_bounds(bounds, dims=axes(bounds,1))
	for (i, c) in enumerate(eachcol(corners))
		e, _ = Evol(c, automaton, nominal)
		ev[:,i,:] = e'[dims,:]
	end
	
	# Compute Hausdorff distance at each time step
	H = Array{Float64}(undef, size(ev, 3))
	for t in axes(ev, 3)
		dist = pairwise(metric, reachable_corners[:,:,t], ev[:,:,t])
		H_row = maximum(minimum(dist, dims=1))
		H_col = maximum(minimum(dist, dims=2))
		H[t] = maximum((H_row, H_col))
	end
	H
end

# ╔═╡ 6d98e67e-e933-4448-8ec1-90cc0ac6f96a
let
	automaton = strat_map[sim_strategy_trees](sysd, K, max_miss_trees)
	d = deviation(automaton, Augment(bounds, automaton), all_bounds, dims=[2])[1:time_trees+1]
	i = argmax(d)
	v = maximum(d)
	
	plot(xlabel="Time", ylabel="Deviation", guidefontsize=18, tickfontsize=14, legendfontsize=14)
	if v < 1e3
		#=uls = CSV.File("/home/cghobbs/grad/Jittery-Scheduler/uls_zk.csv", header=false)
		uls_mat = Tables.matrix(uls)
		plot!(0:size(uls_mat,2)-1, uls_mat[1,:], label="Uncertain Linear Systems", linewidth=2)
		rr = CSV.File("/home/cghobbs/grad/Jittery-Scheduler/fsm_zk_3.csv", header=false)
		rr_mat = Tables.matrix(rr)
		plot!(0:size(rr_mat,2)-1, rr_mat[1,:], label="Recurrence Relations", linestyle=:dash, linewidth=2)=#
		plot!(0:size(d,1)-1, d, label="Bounded Runs", #=linestyle=:dashdot,=# linewidth=2)
		plot!([0, 75, 150], [0.1, 0.1, 0.1], label="", color=:black, series_annotations=["", text("Safety bound", :bottom), ""])
	else
		# If the max is big, use a log scale.  Don't plot the first two points because they're guaranteed to be zero.
		plot!(2:size(d,1)-1, d[3:end], label="Deviation", yaxis=:log)
	end
	#scatter!([i-1, uls_mat[2,1]], [v, uls_mat[1,Int64(uls_mat[2,1])+1]], label="", color=:black, series_annotations=[text(" ($(i-1), $(@sprintf("%.2f", v))) ", (i < n_trees*t_trees/2) ? :left : :right, 14), text(" ($(Int64(uls_mat[2,1])), $(@sprintf("%.2f", uls_mat[1,Int64(uls_mat[2,1])+1]))) ", :left, :bottom, 14)])
	scatter!([i-1], [v], label="", color=:black, series_annotations=[text(" ($(i-1), $(@sprintf("%.4f", v))) ", (i < n_trees*t_trees/2) ? :left : :right, 14)])
	scatter!([11-1], [d[11]], label="", color=:black, series_annotations=[text(" ($(11-1), $(@sprintf("%.4f", d[11]))) ", (i < n_trees*t_trees/2) ? :left : :right, 14)])
	scatter!([21-1], [d[21]], label="", color=:black, series_annotations=[text(" ($(21-1), $(@sprintf("%.4f", d[21]))) ", (i < n_trees*t_trees/2) ? :left : :right, 14)])
	plot!(legend=:topright)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
ControlSystems = "a6e380b2-a6ca-5380-bf3e-84a91bcd477e"
Distances = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
OffsetArrays = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[compat]
CSV = "~0.10.2"
ControlSystems = "~0.11.12"
Distances = "~0.10.7"
Distributions = "~0.25.41"
LaTeXStrings = "~1.3.0"
OffsetArrays = "~1.10.8"
Plots = "~1.22.2"
PlutoUI = "~0.7.11"
QuadGK = "~2.4.2"
Tables = "~1.6.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "1d6835607e9f214cb4210310868f8cf07eb0facc"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.1.34"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "5e98d6a6aa92e5758c4d58501b7bf23732699fa3"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.2"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[CPUSummary]]
deps = ["Hwloc", "IfElse", "Static"]
git-tree-sha1 = "ba19d1c8ff6b9c680015033c66802dd817a9cf39"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.1.7"

[[CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "9519274b50500b8029973d241d32cfbf0b127d97"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.2"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "54fc4400de6e5c3e27be6047da2ef6ba355511f8"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.6"

[[ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "6b6f04f93710c71550ec7e16b650c1b9a612d0b6"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.16.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[ControlSystems]]
deps = ["Colors", "DSP", "DelayDiffEq", "DiffEqCallbacks", "IterTools", "LaTeXStrings", "LinearAlgebra", "MacroTools", "MatrixEquations", "MatrixPencils", "OrdinaryDiffEq", "Plots", "Polynomials", "Printf", "Random", "SparseArrays"]
git-tree-sha1 = "d8c59dbc710104fac6a1e8ee92fca14ffdadf759"
uuid = "a6e380b2-a6ca-5380-bf3e-84a91bcd477e"
version = "0.11.12"

[[DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "fe2287966e085df821c0694df32d32b6311c6f4c"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.4"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "NonlinearSolve", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "UnPack"]
git-tree-sha1 = "fd0ef97b21b6eea22a917ada02ebe38a85e08197"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.33.0"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[DiffEqBase]]
deps = ["ArrayInterface", "ChainRulesCore", "DataStructures", "Distributions", "DocStringExtensions", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "IterativeSolvers", "LabelledArrays", "LinearAlgebra", "Logging", "MuladdMacro", "NonlinearSolve", "Parameters", "PreallocationTools", "Printf", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "StaticArrays", "Statistics", "SuiteSparse", "ZygoteRules"]
git-tree-sha1 = "ad27076f769f812e64310fcaf03658b68259cb85"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.71.0"

[[DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "LinearAlgebra", "NLsolve", "OrdinaryDiffEq", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "e57ecaf9f7875714c164ccca3c802711589127cf"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.20.1"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "9bc5dac3c8b6706b58ad5ce24cffd9861f07c94f"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.9.0"

[[Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "5863b0b10512ed4add2b5ec07e335dc6121065a5"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.41"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[ExponentialUtilities]]
deps = ["ArrayInterface", "LinearAlgebra", "Printf", "Requires", "SparseArrays"]
git-tree-sha1 = "3e1289d9a6a54791c1ee60da0850f4fd71188da6"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.11.0"

[[ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "463cb335fa22c4ebacfd1faba5fde14edb80d96c"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.5"

[[FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[FastBroadcast]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "26be48918640ce002f5833e8fc537b2ba7ed0234"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.1.8"

[[FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "04d13bfa8ef11720c24e4d840c0033d145537df7"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.17"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8756f9935b7ccc9064c6eef0bff0ad643df733a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.7"

[[FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "6eae72e9943d8992d14359c32aed5f892bda1569"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.10.0"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "1bd6fc0c344fc0cbee1f42f8d2e7ec8253dda2d2"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.25"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[FunctionWrappers]]
git-tree-sha1 = "241552bc2209f0fa068b6415b1942cc0aa486bcc"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.2"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "0c603255764a1fa0b61752d2bec14cfbd18f7fe8"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+1"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "d189c6d2004f63fd3c91748c458b09f26de0efaa"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.61.0"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "aa22e1ee9e722f1da183eb33370df4c1aeb6c2cd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.63.1+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "d727758173afef0af878b29ac364a0eca299fc6b"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.5.1"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "3965a3216446a6b020f0d48f1ba94ef9ec01720d"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.6"

[[Hwloc]]
deps = ["Hwloc_jll"]
git-tree-sha1 = "92d99146066c5c6888d5a3abc871e6a214388b91"
uuid = "0e44f5e4-bd66-52a0-8798-143a42290a1d"
version = "2.0.0"

[[Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d8bccde6fc8300703673ef9e1383b11403ac1313"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.7.0+0"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "8d70835a3759cdd75881426fced1508bb7b7e1b6"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.1"

[[IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "323a38ed1952d30586d0fe03412cde9399d3618b"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.5.0"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "22df5b96feef82434b07327e2d3c770a9b21e023"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "41158dee1d434944570b02547d404e075da15690"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.7.0"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LinearMaps]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "dbb14c604fc47aa4f2e19d0ebb7b6416f3cfa5f5"
uuid = "7a12625a-238d-50fd-b39a-03d52299707e"
version = "3.5.1"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "5455aef09b40e5020e1520f551fa3135040d4ed0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+2"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[ManualMemory]]
git-tree-sha1 = "9cb207b18148b2199db259adfa923b45593fe08e"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.6"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MatrixEquations]]
deps = ["LinearAlgebra", "LinearMaps"]
git-tree-sha1 = "49b86537cefbfbf4b86afa3cc17226f1781d4389"
uuid = "99c1a7ee-ab34-5fd5-8076-27c950a045f4"
version = "2.1.0"

[[MatrixPencils]]
deps = ["LinearAlgebra", "Polynomials", "Random"]
git-tree-sha1 = "caf22dd55c613095638a113b39f38e20016e2c36"
uuid = "48965c70-4690-11ea-1f13-43a2532b2fa8"
version = "1.6.7"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "29714d0a7a8083bba8427a4fbfb00a540c681ce7"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.3"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MuladdMacro]]
git-tree-sha1 = "c6190f9a7fc5d9d5915ab29f2134421b12d24a68"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.2"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "73deac2cbae0820f43971fad6c08f6c4f2784ff2"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.3.2"

[[NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[NaNMath]]
git-tree-sha1 = "f755f36b19a5116bb580de457cda0c140153f283"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.6"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[NonlinearSolve]]
deps = ["ArrayInterface", "FiniteDiff", "ForwardDiff", "IterativeSolvers", "LinearAlgebra", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "UnPack"]
git-tree-sha1 = "b61c51cd5b9d8b197dfcbbf2077a0a4e1505278d"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "0.3.14"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[OrdinaryDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "Logging", "MacroTools", "MuladdMacro", "NLsolve", "Polyester", "RecursiveArrayTools", "Reexport", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "73d43522c1c756c739ef70a51941b0b028682994"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "5.60.2"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "ee26b350276c51697c9c2d88a072b339f9f03d73"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.5"

[[Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "92f91ba9e5941fc781fecf5494ac1da87bdac775"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.0"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "6f1b25e8ea06279b5689263cc538f51331d7ca17"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.3"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "e7523dd03eb3aaac09f743c23c1a553a8c834416"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.22.7"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "ae6145ca68947569058866e443df69587acc1806"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.32"

[[Polyester]]
deps = ["ArrayInterface", "IfElse", "ManualMemory", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities", "VectorizationBase"]
git-tree-sha1 = "3ced65f2f182e5b5335a573eaa98f883eba3678b"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.3.9"

[[Polynomials]]
deps = ["Intervals", "LinearAlgebra", "MutableArithmetics", "RecipesBase"]
git-tree-sha1 = "f184bc53e9add8c737e50fa82885bc3f7d70f628"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "2.0.24"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "db3a23166af8aebf4db5ef87ac5b00d36eb771e2"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.0"

[[PreallocationTools]]
deps = ["ArrayInterface", "ForwardDiff", "LabelledArrays"]
git-tree-sha1 = "361c1f60ffdeeddf02f36b463ab8b138194e5f25"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.1.1"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "7ad0dfa8d03b7bcf8c597f59f5292801730c55b8"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.1"

[[RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "ChainRulesCore", "DocStringExtensions", "FillArrays", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "5144e1eafb2ecc75765888a4bdcd3a30a6a08b14"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.24.1"

[[RecursiveFactorization]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6761a5d1f9646affb2a369ff932841fff77934a3"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.1.0"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "40c1c606543c0130cd3673f0dd9e11f2b5d76cd0"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.26.0"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "15dfe6b103c2a993be24404124b8791a09460983"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.11"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "fca29e68c5062722b5b4435594c3d1ba557072a3"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.7.1"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "75c89362201983c500dd34923b015dbecdae7a90"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.20.0"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e6bf188613555c78062842777b116905a9f9dd49"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.0"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "a8f30abc7c64a39d389680b74e749cf33f872a70"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.3.3"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "2884859916598f974858ff01df7dfc6c708dd895"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.3"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "d88665adc9bcf45903013af0982e2fd05ae3d0a6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "51383f2d367eb3b444c961d485c565e4c0cf4ba0"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.14"

[[StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f35e1879a71cca95f4826a14cdbf0b9e253ed918"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.15"

[[StrideArraysCore]]
deps = ["ArrayInterface", "ManualMemory", "Requires", "ThreadingUtilities", "VectorizationBase"]
git-tree-sha1 = "9ab16bda5fe1212e0af0bea80f1d11096aeb3248"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.1.18"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "d21f2c564b21a202f4677c0fba5b5ee431058544"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.4"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "884539ba8c4584a3a8173cb4ee7b61049955b79c"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.4.7"

[[TimeZones]]
deps = ["Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "0f1017f68dc25f1a0cb99f4988f78fe4f2e7955f"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.7.1"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "Hwloc", "IfElse", "Libdl", "LinearAlgebra", "Static"]
git-tree-sha1 = "5e6e23728d6c8d26d2826f6cb2cd21892a958a43"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.20.38"

[[VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "66d72dc6fcc86352f01676e8f0f698562e60510f"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.23.0+0"

[[WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "c69f9da3ff2f4f02e811c3323c22e5dfcb584cfa"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.1"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─afeab312-2d5c-11ec-36d7-8fa4f9a66fa8
# ╠═e0bea5c9-c711-4703-ac3a-4df5e0042e6d
# ╟─89dbc50b-2ce7-46c7-bdec-73166abef4ae
# ╠═ae811d13-c813-40f3-ba7a-878568291252
# ╟─7431a2e3-5b31-43b0-a40b-ab0a977e919d
# ╠═dac9640f-eac8-401d-aada-62c864a758ef
# ╟─264e5c66-5a92-4a4c-8c0c-056e077856ab
# ╠═e86df2ab-534b-471a-90b7-43c7f980abff
# ╟─f5a869b9-b2a8-41a7-a503-8523ae5842c6
# ╟─8c489d40-d44d-48d0-845a-8afbd5161c0b
# ╟─1d610dac-8989-4bba-99fd-f248fcf10e8c
# ╟─3d095e19-4f82-4c68-a104-2b7d94467849
# ╟─a92710a5-9606-4ba6-88a4-962bca2a2c3e
# ╟─fdd6db9c-2d42-4073-b7ca-22d6396b86ab
# ╟─fcb00b19-7eb1-40ca-b7a6-8fbf61397172
# ╠═2481a5d6-71dd-4530-80c0-a6b4eae4a13b
# ╠═9cddd940-d4e3-44e1-ac01-de314f9e8517
# ╟─1317e2b7-2b8a-4bb8-aee8-bcd6f8d38a03
# ╠═bddaa61b-851f-4759-889a-b71adff1ca57
# ╟─9772aa30-71d2-4083-a3ad-d842b0fd2b13
# ╠═03aa61fd-67ca-45f5-9eef-7464901598cc
# ╟─0cd0ceff-ebc5-4069-9bb9-14fcc667a1ef
# ╠═acc7b261-8713-4565-9589-f59af764175c
# ╟─e07b2e53-e0eb-4e00-a25f-869d98d04cbe
# ╠═1faf171b-9a38-42dc-867f-4fd576228bde
# ╠═e1c05296-ddb4-4da0-b972-d952ab3a5ef2
# ╟─fc69b358-f22c-450f-9ab0-5ed9ed564207
# ╠═72650fe6-cd13-46f1-9ce7-10831a27c396
# ╟─06c7c6d1-719e-4c9a-bf28-939288c335e9
# ╠═ff4a5ea5-96f5-4054-ab8b-681aef3ad017
# ╟─9ed43a4f-f1c3-4299-81ff-2345f4b9338a
# ╟─bb865586-2365-410a-be80-04b65c935381
# ╟─d9d3bb42-43e0-4b68-8135-6447294de688
# ╠═5406484b-7c43-41df-8f63-3e1d16fde3aa
# ╟─9179f1a5-26df-4e20-8e88-cbf6df6c24e6
# ╠═f7161723-6e58-4ddd-ac74-f2b1488e6d2e
# ╟─ffdc8ef4-2eea-4ed3-a57d-ea0f14219f7b
# ╠═14847de2-a3a8-4206-8bf7-64a9896fa750
# ╟─8f202fe5-008b-4105-8322-8f595fb6b545
# ╟─f61f27cb-48d3-4f95-8547-5b0a3a2b6394
# ╟─5ae87150-d2cb-4e84-b7c9-bb255b4f36e0
# ╠═97c1f398-a62f-401b-802d-283c8dc64307
# ╠═0eb09c86-d48e-4a1e-b7ea-5f1c96df2667
# ╟─8e1164d6-9559-4103-b4e6-577166d0ab2b
# ╟─28932e44-1762-4ca7-9698-80eeaffa99a1
# ╟─da2987de-5bd6-46a8-ba21-58f5d2a2f397
# ╠═d7c5ea66-4b3f-4920-9eda-858466a076b3
# ╠═a7761898-3e69-4371-9315-8282fe07407d
# ╠═b7e416f5-897b-4e45-bd84-a5ef25e7c43d
# ╟─2280b3d6-d3ab-4049-b796-a2a036b55471
# ╟─883adf52-b18c-4f24-9189-601b5a50d7c2
# ╠═728953ad-b271-4ac0-bdc7-e8c71597b092
# ╟─421ecae9-b1b3-4a9e-af01-6cc9befa5884
# ╠═6d98e67e-e933-4448-8ec1-90cc0ac6f96a
# ╟─7c7ea3d4-fe32-4560-b501-7845265dff52
# ╟─0363259a-2096-4825-8cec-bcb906f20a23
# ╟─6baa69a2-74af-4c85-8368-49cfa35628b8
# ╠═a0d254fa-beba-424a-b77a-13c9bfce2db9
# ╠═6a82ce98-0a17-4f13-8979-94c1eccdd5b4
# ╠═eb79be36-73ba-423d-9d77-a42ab4007246
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

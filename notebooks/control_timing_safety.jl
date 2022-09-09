### A Pluto.jl notebook ###
# v0.19.11

#> [frontmatter]
#> title = "ControlTimingSafety.jl Demo"
#> date = "2022-09-09"
#> description = "An illustration of the methods from Hobbs et al., \"Safety Analysis of Embedded Controllers under Implementation Platform Timing Uncertainties,\" TCAD 2022."

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

# ╔═╡ 7e822d4c-3043-11ed-040a-f50d3cf62833
begin
    import Pkg
	# Add the ControlTimingSafety.jl dev package 😃
	# This makes Pluto stop managing packages on its own, so you'll have to have all the packages listed below installed in your default environment.
    Pkg.develop(path="../")
    # Instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

	using Revise
	using ControlTimingSafety
    using Plots, PlutoUI, LaTeXStrings
	using Random, Distributions
	using Distances
	using ControlSystems, LinearAlgebra
end

# ╔═╡ b80a4d32-2a7b-4a59-8d04-d515471e86b4
md"""
# ControlTimingSafety.jl Demo

This notebook provides a demonstration of the ControlTimingSafety package, which implements the transducer automata and bounded runs iteration algorithm from [^hobbssafety].
$(TableOfContents())

!!! note

	This notebook is designed to be run directly from the ControlTimingSafety.jl repository, and disables Pluto's package management to achieve this.  Please take a look at the Pkg cell in [the appendix](#35b0c50d-4452-4003-b1f8-a6552b04df0b) and make sure all the necessary packages are installed in your default environment.
"""

# ╔═╡ c46333f1-fcba-4d8e-9d8c-015fef5a6eb9
md"""
## Demonstration of Transducer Automata

We begin by creating a two-dimensional LTI system, namely, the RC Network model from [^hobbssafety], in turn based on a circuit in [^gabelroberts].
"""

# ╔═╡ 8c2d4565-5326-4fd3-8c3b-1491de1040ff
sys = let
	r₁ = 100000
	r₂ = 500000
	r₃ = 200000
	c₁ = 0.000002
	c₂ = 0.000010
	A = [-1/c₁ * (1/r₁ + 1/r₂)  1/(r₂*c₁)
	     1/(r₂*c₂)              -1/c₂ * (1/r₂ + 1/r₃)]
	B = [1/(r₁*c₁)
	     1/(r₃*c₂)]
	C = [1  -1]
	D = 0
	ss(A, B, C, D)
end

# ╔═╡ c6481e62-b578-4270-b899-c322dae1d009
md"""
We then discretize the model with a constant period $h$ (in seconds), and no sensor-to-actuator delay.

|     |                                                        |
|----:|:-------------------------------------------------------|
| $h$ | $(@bind ctrl_delay NumberField(0:0.05:1, default=0.1)) |
"""

# ╔═╡ b7f9fb0e-00ae-402b-ac7f-fc0ac7c70933
sysd = c2d(sys, ctrl_delay)

# ╔═╡ 9a8d3bdf-12ba-4b96-a230-429b88d9ef84
md"""
We next create a version of this model with a sensor-to-actuator delay of one sampling period.  This is achieved using an augmented state space $z[t] = \big[x[t]^T\ \ u[t-1]\big]^T$.  The $u[t-1]$ entry is stored using the new $B$ matrix.  The original $B$ matrix from `sysd` appears in the third column of the new $A$ matrix, so that the control input can be correctly applied to the system state.
"""

# ╔═╡ 12823d77-c26a-4d34-98cd-81f01e4b42e9
sysd_delay = ss([sysd.A sysd.B; 0 0 0], [0; 0; 1], [sysd.C 0], sysd.D, sysd.Ts)

# ╔═╡ e76a5efa-2d14-47eb-8dac-b9327ac208c6
md"""
We now create an optimal controller for the delayed system using LQR.  The default values (state cost 2.0, control cost 1.0) gives the controller used in [^hobbssafety].  Try changing the values to see the effects below!

|              |                                                        |
|-------------:|:-------------------------------------------------------|
|   State cost | $(@bind state_cost NumberField(0:0.1:10, default=2))   |
| Control cost | $(@bind control_cost NumberField(0:0.1:10, default=1)) |
"""

# ╔═╡ 4f9ffd6b-a409-4e9b-ba24-57b93a645973
K = let
	Q = [state_cost 0 0;
		 0 state_cost 0;
		 0 0 state_cost]
	R = [control_cost;;]
	lqr(sysd_delay, Q, R)
end

# ╔═╡ 36cafaaa-e244-4404-8837-620a4285e76e
md"""
### Example 1: Effects of deadline misses

We are now ready to simulate the delayed system under a user-selectable strategy, both with deadline misses and without.
The nominal behavior, where no deadlines are missed, is shown in blue; other colors show trajectories with probabilistic deadline misses.  Try changing the miss strategy and simulation parameters below!

|                          |      |
|-------------------------:|:-----|
|                 Strategy | $(@bind sim_strategy Select(strat_names)) |
| Orange $P(\mathrm{hit})$ | $(@bind hit_prob_1 Slider(0:0.05:1, default=0.90, show_value=true)) |
|    Red $P(\mathrm{hit})$ | $(@bind hit_prob_2 Slider(0:0.05:1, default=0.65, show_value=true)) |
| Purple $P(\mathrm{hit})$ | $(@bind hit_prob_3 Slider(0:0.05:1, default=0.05, show_value=true)) |
|          Simulation time | $(@bind sim_time Slider(5:5:500, default=500, show_value=true)) |
|                 $x_1[0]$ | $(@bind x_0_1 NumberField(-10:0.1:10, default=1)) |
|                 $x_2[0]$ | $(@bind x_0_2 NumberField(-10:0.1:10, default=1)) |

|                                                |
|:----------------------------------------------:|
| $(@bind recompute_plot Button("🎲 Randomize")) |
"""

# ╔═╡ 8dde22c7-0af1-4ffa-9cc6-fc8aa4f0d096
md"""
### Example 2: Random sampling of action sequences

We next show a second example in which the same system is run for a large number of hit/miss sequences, selected at random with the constraint of no more than $n$ consecutive misses.  This example creates a plot similar to Figure 1 in [^hobbssafety].
"""

# ╔═╡ 06d49fa1-942c-4871-9516-291a9de2fd1b
md"""
|                   |                                                              |
|------------------:|:-------------------------------------------------------------|
|          Strategy | $(@bind strategy_2 Select(strat_names))                      |
|    Max misses $n$ | $(@bind maxmiss_2 Slider(-1:16, default=5, show_value=true)) |
| Number of samples | $(@bind n_samples_2 NumberField(5:5:100, default=100))       |
|   Simulation time | $(@bind H_2 Slider(5:5:500, default=100, show_value=true))   |
|          $x_1[0]$ | $(@bind x₀_1_2 NumberField(-10:0.1:10, default=-1))           |
|          $x_2[0]$ | $(@bind x₀_2_2 NumberField(-10:0.1:10, default=-1))           |
|    Camera azimuth | $(@bind azimuth_2 Slider(-180:1:180, default=35, show_value=true)) |
|  Camera elevation | $(@bind elevation_2 Slider(-90:1:90, default=25, show_value=true)) |

|                                                 |
|:-----------------------------------------------:|
| $(@bind recompute_plot_2 Button("🎲 Randomize")) |
"""

# ╔═╡ ba54ee0e-4919-4c35-a48b-4b53511eb4be
md"""
The code to compute this plot follows in the remainder of this section.  We first simply compute the nominal trajectory consisting of all deadline hits.
"""

# ╔═╡ c822e252-552b-48bd-8269-266d3d30f7ea
nominal_2 = let	
	T = strat_map[strategy_2](sysd, K)
	z₀ = augment(T, [x₀_1_2, x₀_2_2])
	evol(T, z₀, ones(Int64, H_2))
end

# ╔═╡ a057b28a-de4f-4b8f-87c1-4d32e3f38242
md"""
Next, we compute all the random trajectories.  An extra check is done to bias towards more interesting-looking trajectories, making the sampling in the plot not truly uniform.
"""

# ╔═╡ 5926e7e5-1e4b-408e-a02a-6941bfb95845
md"""
We next compute the distance between each point of each random trajectory and the corresponding point of the nominal trajectory.
"""

# ╔═╡ 79f22b98-8bf7-465b-b900-81b22a4b0ca8
md"""
Using this distance, we compute the deviation between each random trajectory and the nominal one.
"""

# ╔═╡ d9ed868d-ee5b-4470-8edc-3b6adb429500
md"""
In order for this random sampling to work, we need a way to perform this uniform sampling.  Fortunately, Julia provides a nice API for sampling random values from different types of data and distributions.  We next create a `Sampler` type for uniformly sampling strings with at most $k$ misses in a row, of length $H$.
"""

# ╔═╡ ea5c98dc-0496-43ec-a8f1-2566e0590c93
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

# ╔═╡ 0bb32ddb-77c6-4196-90cd-9ecd27b9bf17
md"""
We next define a method of the function `Random.rand` that allows us to sample from this distribution.
"""

# ╔═╡ 275a5d3f-2936-408d-9705-a807df9f0b79
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
	2 .- ret
end

# ╔═╡ 2ccd0233-8503-4fcb-be27-f8bf78df3219
let
	# Recompute when the button is clicked
	recompute_plot
	
	automaton = strat_map[sim_strategy](sysd, K)
	
	hit_prob = [1, hit_prob_1, hit_prob_2, hit_prob_3]
	line_colors = [:blue, :orange, :red, :purple]
	z_0 = augment(automaton, [x_0_1, x_0_2])
	
	z = [0]
	plot(title="Trajectories for different hit probabilities", xlabel=L"x_1", ylabel=L"x_2", legend=:bottomright)
	for i = size(hit_prob, 1):-1:1
	    seq = (rand(sim_time) .>= hit_prob[i]) .+ 1
		z = evol(automaton, z_0, seq)
	    plot!(z[:,1], z[:,2], label=L"P(\mathrm{hit}) = %$(hit_prob[i])", linecolor=line_colors[i])
	end
	plot!()
end

# ╔═╡ a6b57e9d-0d61-41f1-96ff-4502daa0ba4a
trj_2 = let
	# Recompute when the button is clicked
	recompute_plot_2

	T = strat_map[strategy_2](sysd, K)
	z₀ = augment(T, [x₀_1_2, x₀_2_2])

	# Compute random trajectories
	sp = SamplerAutomatonInput(maxmiss_2, H_2)
	seq = ones(Int64, H_2)
	trj = Array{Float64}(undef, n_samples_2, H_2+1, size(T.Φ[1], 1))
	max_possible = (H_2 / (maxmiss_2 + 1)) * maxmiss_2 + H_2 % (maxmiss_2 + 1)
	for i = axes(trj, 1)
		accepted = false
		while !accepted
			seq = rand(sp)
			if rand(Float64) > sum(seq .== 2) / (max_possible * .5)
				accepted = true
			end
		end
		trj[i,:,:] = evol(T, z₀, seq)
	end
	trj
end

# ╔═╡ 95f31b7c-051f-4b15-8edd-0628faf94b5d
dist_2 = let	
	dist = Array{Float64}(undef, n_samples_2, H_2+1)
	for i in axes(trj_2, 1)
		dist[i,:] = colwise(Euclidean(), trj_2[i,:,1:2]', nominal_2[:,1:2]')
	end
	dist
end

# ╔═╡ f2ddd147-d1c7-4453-8997-76752c5bc515
dev_2 = maximum(dist_2, dims=2)

# ╔═╡ 05221746-efbb-45f5-be39-dbeb09996088
# ╠═╡ skip_as_script = true
#=╠═╡
let
	plot(xlabel=L"x_1", ylabel=L"x_2", zlabel=L"t", legend=:topleft, format=:png, thickness_scaling=2, size=(1360,906), camera=(azimuth_2,elevation_2))
	
	# Plot safety pipe
	θ = LinRange(0, 2π, 40)
	pipe_radius = sort(dev_2, dims=1)[end-1] - 0.01
	circ_x = cos.(θ) * pipe_radius
	circ_y = sin.(θ) * pipe_radius
	for i in axes(nominal_2, 1)
		plot!(circ_x .+ nominal_2[i,1], circ_y .+ nominal_2[i,2], repeat([i,], size(θ,1)), label=(i == 1) ? "Safety pipe" : "", seriestype=[:shape,], c=:lightblue, linecolor=:lightblue)
	end

	# Plot random trajectories
	label_good = false
	# Plot the good trajectories first
	for i = axes(trj_2, 1)
		if dev_2[i] < pipe_radius
			plot!(trj_2[i,:,1], trj_2[i,:,2], 1:H_2+1, label=(!label_good) ? "Random" : "", color=:green, opacity=10/n_samples_2)
			label_good = true
		end
	end
	# Plot the bad trajectories on top of the good ones
	label_bad = false
	nom_crosses = zeros(Float64, H_2+1)
	for i = axes(trj_2, 1)
		if dev_2[i] >= pipe_radius
			t_viol = argmax(dist_2[i,:], dims=1)
			#label_bad || plot!(cos.(θ) * pipe_radius .+ nominal_2[t_viol,1], sin.(θ) * pipe_radius .+ nominal_2[t_viol,2], repeat([t_viol,], size(θ,1)), label="Radius violated", style=:dot, linecolor=:gray)
			crosses = [(d >= pipe_radius) ? 1 : 0 for d in dist_2[i,:]]
			plot!(trj_2[i,:,1], trj_2[i,:,2], 1:H_2+1, label=(!label_bad) ? "Violation" : "", color=:red, markershape=:x, markeralpha=crosses)
			nom_crosses .+= crosses
			label_bad = true
		end
	end

	# Finally, plot the nominal trajectory on top of everything else
	plot!(nominal_2[:,1], nominal_2[:,2], 1:H_2+1, label="Nominal", color=:black, linewidth=3, marker=:x, markeralpha=nom_crosses)
end
  ╠═╡ =#

# ╔═╡ c9856b8f-4341-4e33-9a73-3b1570fb24fa
md"""
## Bounded Tree Reachable Set

In this section, we use a *bounded tree method* to compute a reachable set for the system under any sequence of deadline hits and misses.  In this method, we start with an axis-aligned box initial set, and simulate its evolution for $n$ steps, for each of the $O(|\mathcal{A}|^n)$ possible sequences of scheduler actions.  We then take an axis-aligned bounding box of the resulting sets, and repeat this process.

This algorithm supports running any automaton, including those that implement different deadline miss strategies, and weakly-hard constraints.  This gives us power to model a variety of scheduler behaviors that have been considered in the literature.

The ControlTimingSafety package provides the implementation of this algorithm.  Following [^hobbssafety], it is broken into two functions: one that simulates all possible runs for a bounded length of time, and another that runs the former iteratively.  This section demonstrates both.
"""

# ╔═╡ fbf988b0-7574-4dbb-929f-3d7df117d9dc
md"""
### Example 3: Bounded Runs algorithm

This section illustrates the Bounded Runs algorithm, used to explore all possible sequences of actions for a small time horizon.  The example uses the same plant model and controller from the previous section, with a small initial set of states.

|                         |                                                                             |
|------------------------:|:----------------------------------------------------------------------------|
|                Strategy | $(@bind sim_strategy_bounded_runs Select(strat_names))                      |
|                     $n$ | $(@bind n_bounded_runs Slider(1:16, default=8,  show_value=true)) |
|              Max misses | $(@bind maxmiss_bounded_runs Slider(-1:16, default=5, show_value=true)) |
| Show nominal trajectory | $(@bind show_nom_traj_bounded_runs CheckBox(default=false))                                |
"""

# ╔═╡ 0843dd0e-0091-402f-9f28-fb112232f140
bounds = [2; 2;; 2.5; 2.5]

# ╔═╡ ef358d22-a530-4881-a282-b850d3655427
begin
	points_bounded_runs = let
		automaton = strat_map[sim_strategy_bounded_runs](sysd, K, maxmiss_bounded_runs)
		augbounds = augment(automaton, bounds)
		bounded_runs(automaton, augbounds, n_bounded_runs)
	end
	
	md"""
	The reachable sets are computed in this cell and shown below.  Try enabling the sample trajectories, verifying that they remain inside the boxes for each time step.
	
	!!! warning
	
	    This is effectively a brute-force reachable set calculation for all hit/miss sequences of length $n$, with time complexity $O(|\mathcal{A}|^n)$.  It may take several minutes to compute with large values set on the sliders above.
	"""
end

# ╔═╡ c4a6821b-deb9-4a2b-89db-09b382c72b1c
let
	plot(title="Reachable set for $(n_bounded_runs) time steps", xlabel="x_1", ylabel="x_2", legend=:bottomright)
	merged = merge_bounds(points_bounded_runs)
	for k = 0:n_bounded_runs
		corners = corners_from_bounds(merged[begin+k,:,:], cycle=true, dims=1:2)
		plot!(corners[1,:], corners[2,:], label=L"x[%$(k)]")
	end
	
	if show_nom_traj_bounded_runs
		hsn = hold_skip_next(sysd, K)
	    x = evol(hsn, augment(hsn, [bounds[1,1], bounds[2,1]]), ones(Int64, n_bounded_runs))
		plot!(x[:,1], x[:,2], label="Nominal", linecolor=:blue, marker=:circle)
	    x = evol(hsn, augment(hsn, [bounds[1,1], bounds[2,2]]), ones(Int64, n_bounded_runs))
		plot!(x[:,1], x[:,2], label="Nominal", linecolor=:blue, marker=:circle)
	    x = evol(hsn, augment(hsn, [bounds[1,2], bounds[2,1]]), ones(Int64, n_bounded_runs))
		plot!(x[:,1], x[:,2], label="Nominal", linecolor=:blue, marker=:circle)
	    x = evol(hsn, augment(hsn, [bounds[1,2], bounds[2,2]]), ones(Int64, n_bounded_runs))
		plot!(x[:,1], x[:,2], label="Nominal", linecolor=:blue, marker=:circle)
	end
	plot!()
end

# ╔═╡ 35b0c50d-4452-4003-b1f8-a6552b04df0b
md"""
## Appendix
"""

# ╔═╡ 8b1dbd59-f31a-4833-9628-d62aba43cd28
md"""
### Bibliography

[^gabelroberts]: Robert Gabel and Richard Roberts.
	"Signals and Linear Systems."
	Second edition.
	Wiley, 1980.
[^hobbssafety]: Clara Hobbs, Bineet Ghosh, Shengjie Xu, Parasara Sridhar Duggirala, and Samarjit Chakraborty.
	"Safety Analysis of Embedded Controllers under Implementation Platform Timing Uncertainties."
	TCAD 2022.
[^maggiostability]: Martina Maggio, Arne Hamann, Eckart Mayer-John, and Dirk Ziegenbein.
	"Control-system stability under consecutive deadline misses constraints."
	ECRTS 2020.
	Online: [https://drops.dagstuhl.de/opus/volltexte/2020/12384/pdf/LIPIcs-ECRTS-2020-21.pdf](https://drops.dagstuhl.de/opus/volltexte/2020/12384/pdf/LIPIcs-ECRTS-2020-21.pdf)
"""

# ╔═╡ Cell order:
# ╟─b80a4d32-2a7b-4a59-8d04-d515471e86b4
# ╟─c46333f1-fcba-4d8e-9d8c-015fef5a6eb9
# ╟─8c2d4565-5326-4fd3-8c3b-1491de1040ff
# ╟─c6481e62-b578-4270-b899-c322dae1d009
# ╟─b7f9fb0e-00ae-402b-ac7f-fc0ac7c70933
# ╟─9a8d3bdf-12ba-4b96-a230-429b88d9ef84
# ╟─12823d77-c26a-4d34-98cd-81f01e4b42e9
# ╟─e76a5efa-2d14-47eb-8dac-b9327ac208c6
# ╟─4f9ffd6b-a409-4e9b-ba24-57b93a645973
# ╟─36cafaaa-e244-4404-8837-620a4285e76e
# ╟─2ccd0233-8503-4fcb-be27-f8bf78df3219
# ╟─8dde22c7-0af1-4ffa-9cc6-fc8aa4f0d096
# ╟─06d49fa1-942c-4871-9516-291a9de2fd1b
# ╟─05221746-efbb-45f5-be39-dbeb09996088
# ╟─ba54ee0e-4919-4c35-a48b-4b53511eb4be
# ╟─c822e252-552b-48bd-8269-266d3d30f7ea
# ╟─a057b28a-de4f-4b8f-87c1-4d32e3f38242
# ╟─a6b57e9d-0d61-41f1-96ff-4502daa0ba4a
# ╟─5926e7e5-1e4b-408e-a02a-6941bfb95845
# ╟─95f31b7c-051f-4b15-8edd-0628faf94b5d
# ╟─79f22b98-8bf7-465b-b900-81b22a4b0ca8
# ╟─f2ddd147-d1c7-4453-8997-76752c5bc515
# ╟─d9ed868d-ee5b-4470-8edc-3b6adb429500
# ╠═ea5c98dc-0496-43ec-a8f1-2566e0590c93
# ╟─0bb32ddb-77c6-4196-90cd-9ecd27b9bf17
# ╠═275a5d3f-2936-408d-9705-a807df9f0b79
# ╟─c9856b8f-4341-4e33-9a73-3b1570fb24fa
# ╟─fbf988b0-7574-4dbb-929f-3d7df117d9dc
# ╟─0843dd0e-0091-402f-9f28-fb112232f140
# ╟─ef358d22-a530-4881-a282-b850d3655427
# ╟─c4a6821b-deb9-4a2b-89db-09b382c72b1c
# ╟─35b0c50d-4452-4003-b1f8-a6552b04df0b
# ╠═7e822d4c-3043-11ed-040a-f50d3cf62833
# ╟─8b1dbd59-f31a-4833-9628-d62aba43cd28

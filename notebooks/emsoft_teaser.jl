### A Pluto.jl notebook ###
# v0.19.9

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

# ╔═╡ 19ac7356-34fc-11ed-06ab-67dbac77f1dc
begin
    import Pkg
	# Add the ControlTimingSafety.jl dev package 😃
	# This makes Pluto stop managing packages on its own, so you'll have to have all the packages listed below installed in your default environment.
    Pkg.develop(path="../")
    # Instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

	using Revise
	using ControlTimingSafety
    using Plots, Plots.PlotMeasures, PlutoUI, LaTeXStrings
	using Random, Distributions
	using Distances
	using ControlSystems, LinearAlgebra
	using Interpolations
end

# ╔═╡ a25ca2c0-b14a-4e26-9de7-4d17b4d40fe8
md"""
# EMSOFT 2022 Teaser Animations

This notebook contains the code to generate the animations in the EMSOFT 2022 teaser video for the paper "Safety Analysis of Embedded Controllers under
Implementation Platform Timing Uncertainties" **TODO: Link to the finished video**.  They're all created as separate video clips, and edited together with a voiceover in **TODO an open-source video editing program**.
"""

# ╔═╡ bde0e250-f700-40be-b9dc-eae26763b1d1
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

# ╔═╡ fe7ab3cf-5bf4-4af8-896d-af1838a089fd
md"*h* = $(@bind ctrl_delay NumberField(0:0.05:1, default=0.1))"

# ╔═╡ 6dac77d6-090a-45b6-acdf-0bcda8445f03
sysd = c2d(sys, ctrl_delay)

# ╔═╡ df7b4e48-122f-4676-be5f-d53e3ed3c354
sysd_delay = ss([sysd.A sysd.B; 0 0 0], [0; 0; 1], [sysd.C 0], sysd.D, sysd.Ts)

# ╔═╡ 63d725ee-734e-4953-abdc-4ba398534bf6
md"""
|              |                                                        |
|-------------:|:-------------------------------------------------------|
|   State cost | $(@bind state_cost NumberField(0:0.1:10, default=2))   |
| Control cost | $(@bind control_cost NumberField(0:0.1:10, default=1)) |
"""

# ╔═╡ 348989fe-4d30-4fbd-a2f3-c0db26b4eed6
K = let
	Q = [state_cost 0 0;
		 0 state_cost 0;
		 0 0 state_cost]
	R = [control_cost;;]
	lqr(sysd_delay, Q, R)
end

# ╔═╡ c1f773c8-aaed-4da8-bb17-132ed5b82ebd
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

# ╔═╡ af628411-9be0-4ed4-b31f-976d5d110344
animsteps(f) = (H_2+1) .^ ((0:f) ./ f)

# ╔═╡ 427002a6-5f53-42b1-be5d-ddef2932e000
md"## Computation of trajectories"

# ╔═╡ 706d1c80-4b29-4f7a-b85a-7559bc579507
nominal_2 = let	
	T = strat_map[strategy_2](sysd, K)
	z₀ = augment(T, [x₀_1_2, x₀_2_2])
	evol(T, z₀, ones(Int64, H_2))
end

# ╔═╡ 1c0d8251-b8a8-4e68-a825-fd8174dc6a04
md"## Automaton input sampling code"

# ╔═╡ 288202fc-c40f-4ddc-ba88-3ded8b3f7dd3
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

# ╔═╡ 11c07e8d-44b0-400f-a1aa-dcb7694cec9b
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
		rand_bit = Random.rand(rng, Binomial(1, Float64(prob_one)))
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

# ╔═╡ c69777bc-a7d0-4ed4-9837-0cf166fc9dd3
trj_2 = let
	# Recompute when the button is clicked
	recompute_plot_2

	T = strat_map[strategy_2](sysd, K)
	z₀ = augment(T, [x₀_1_2, x₀_2_2])

	# Compute random trajectories
	xoshiro = Xoshiro(31337)
	sp = SamplerAutomatonInput(maxmiss_2, H_2)
	seq = ones(Int64, H_2)
	trj = Array{Float64}(undef, n_samples_2, H_2+1, size(T.Φ[1], 1))
	max_possible = (H_2 / (maxmiss_2 + 1)) * maxmiss_2 + H_2 % (maxmiss_2 + 1)
	for i = axes(trj, 1)
		accepted = false
		while !accepted
			seq = rand(xoshiro, sp)
			if rand(xoshiro, Float64) > sum(seq .== 2) / (max_possible * .5)
				accepted = true
			end
		end
		trj[i,:,:] = evol(T, z₀, seq)
	end
	trj
end

# ╔═╡ b730daed-6bb1-47a1-a4b8-ee1e606f3dda
dist_2 = let	
	dist = Array{Float64}(undef, n_samples_2, H_2+1)
	for i in axes(trj_2, 1)
		dist[i,:] = colwise(Euclidean(), trj_2[i,:,1:2]', nominal_2[:,1:2]')
	end
	dist
end

# ╔═╡ 53434f2b-1733-4e5b-b174-821c4717fd4d
dev_2 = maximum(dist_2, dims=2)

# ╔═╡ ffa3854f-3633-4b54-a8b7-c449137c54b9
pipe_radius = sort(dev_2, dims=1)[end-1] - 0.01

# ╔═╡ 46a30bcf-85a0-44a9-9428-fe7206a4f10d
strat_var_2 = let
	# Recompute when the button is clicked
	recompute_plot_2

	automata = [hold_kill(sysd, K), hold_skip_next(sysd, K), zero_kill(sysd, K), zero_skip_next(sysd, K)]
	z₀ = [augment(a, [x₀_1_2, x₀_2_2]) for a in automata]

	# Compute random trajectories
	xoshiro = Xoshiro(31337)
	sp = SamplerAutomatonInput(maxmiss_2, H_2)
	seq = rand(xoshiro, sp)
	trj = Array{Float64}(undef, length(automata), H_2+1, 2)
	max_possible = (H_2 / (maxmiss_2 + 1)) * maxmiss_2 + H_2 % (maxmiss_2 + 1)
	for i in axes(trj, 1)
		trj[i,:,:] = evol(automata[i], z₀[i], seq)[:,1:2]
	end
	trj
end

# ╔═╡ 7ac8bae8-feb6-4789-b487-365e970656c7
md"## Appendix"

# ╔═╡ c435311f-72a7-4dc9-915d-ff5112686bde
md"### Colors"

# ╔═╡ 6b58db92-abc1-4123-a805-48d40c0914f0
carolina_blue = colorant"#4B9CD3"

# ╔═╡ e9b228fa-4674-4e18-a04c-196642246ab3
navy = colorant"#13294B"

# ╔═╡ 53f962e1-8988-47b5-9e24-ba3b4b273768
black = colorant"#151515"

# ╔═╡ 40f14ec9-cce9-4439-aa20-ead6946acb07
white = colorant"#FFFFFF"

# ╔═╡ d623d388-b218-45ef-b1dd-00c6d0ba2f7d
let
	itp = interpolate(nominal_2, BSpline(Linear()))

	anim = @animate for i in animsteps(30)
		fi = Int64(floor(i))
		p = plot(
			xlabel=L"x_1", ylabel=L"x_2", zlabel=L"t", guidefontsize=14,
			xlims=(-1.25,0.25), ylims=(-1.25,0.25), zlims=(0,100),
			thickness_scaling=2, size=(1920,1080),
			camera=(azimuth_2,elevation_2),
			legend=:topleft, legendfontsize=14, legendfontfamily="helvetica bold", legend_font_color=white, background_color_legend=navy,
			fontfamily="computer modern",
			background_color=carolina_blue, foreground_color=navy, foreground_color_guide=black, margin=-60px)
		
		# Finally, plot the nominal trajectory on top of everything else
		plot!(p, [itp[begin:fi,1]..., itp[i,1]] , [itp[begin:fi,2]..., itp[i,2]], [0:fi-1..., i-1], label="Nominal", color=white, linewidth=3)
		@info i
	end
	webm(anim, fps=30)
end

# ╔═╡ 801c8332-292d-4c8d-8098-93574cc754b6
gray = colorant"#F8F8F8"

# ╔═╡ e9c25a03-2619-4556-a187-d0c02de54f1a
basin_slate = colorant"#4F758B"

# ╔═╡ fc1073f6-d2c1-4e50-8d9b-1968ac88f1d2
campus_sandstone = colorant"#F4E8DD"

# ╔═╡ 319485a7-b5ac-4de0-8b6c-b3c1d4a44cce
longleaf_pine = colorant"#00594C"

# ╔═╡ 33a6152a-2fcd-4f88-9d04-1b7074149a34
azalea_pink = colorant"#EF426F"

# ╔═╡ 998ce331-09c5-4a6f-abe1-bf61c81bf6e7
tile_teal = colorant"#00A5AD"

# ╔═╡ 71032cbe-91e2-4056-808e-4917244c4747
sunburst_yellow = colorant"#FFD100"

# ╔═╡ 73084196-b1f1-4bc7-9190-989b2dcae677
davie_green = colorant"#C4D600"

# ╔═╡ 24a7cc48-5c21-4e7b-9120-2a9939817bb7
let
	plot(format=:png,
		xlabel=L"x_1", ylabel=L"x_2", zlabel=L"t", guidefontsize=14,
		xlims=(-1.25,0.25), ylims=(-1.25,0.25), zlims=(0,100),
		thickness_scaling=2, size=(1920,1080),
		camera=(azimuth_2,elevation_2),
		legend=:topleft, legendfontsize=14, legendfontfamily="helvetica bold", legend_font_color=white, background_color_legend=navy,
		fontfamily="computer modern",
		background_color=carolina_blue, foreground_color=navy, foreground_color_guide=black, margin=-60px)
	
	# Plot safety pipe
	θ = LinRange(0, 2π, 40)
	circ_x = cos.(θ) * pipe_radius
	circ_y = sin.(θ) * pipe_radius
	for i in axes(nominal_2, 1)
		plot!(circ_x .+ nominal_2[i,1], circ_y .+ nominal_2[i,2], repeat([i-1,], size(θ,1)), label=(i == 1) ? "Safety pipe" : "", seriestype=[:shape,], c=:lightblue, linecolor=basin_slate)
	end

	# Plot random trajectories
	label_good = false
	# Plot the good trajectories first
	for i = axes(trj_2, 1)
		if dev_2[i] < pipe_radius
			plot!(trj_2[i,:,1], trj_2[i,:,2], 0:H_2, label=(!label_good) ? "Misses" : "", color=davie_green#=, opacity=10/n_samples_2=#)
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
			plot!(trj_2[i,:,1], trj_2[i,:,2], 0:H_2, label=(!label_bad) ? "Violation" : "", color=azalea_pink, markershape=:xcross, markeralpha=crosses)
			nom_crosses .+= crosses
			label_bad = true
		end
	end

	# Finally, plot the nominal trajectory on top of everything else
	plot!(nominal_2[:,1], nominal_2[:,2], 0:H_2, label="Nominal", color=white, linewidth=3, marker=:xcross, markeralpha=nom_crosses)
end

# ╔═╡ 4aa065f7-2c8a-4f6e-bd0e-102318d13a99
let
	anim = @animate for i in range(1, stop=size(trj_2, 1), length=30)
		fi = Int64(floor(i))
		p = plot(
			xlabel=L"x_1", ylabel=L"x_2", zlabel=L"t", guidefontsize=14,
			xlims=(-1.25,0.25), ylims=(-1.25,0.25), zlims=(0,100),
			thickness_scaling=2, size=(1920,1080),
			camera=(azimuth_2,elevation_2),
			legend=:topleft, legendfontsize=14, legendfontfamily="helvetica bold", legend_font_color=white, background_color_legend=navy,
			fontfamily="computer modern",
			background_color=carolina_blue, foreground_color=navy, foreground_color_guide=black, margin=-60px)

		# Plot random trajectories
		# Fudge the legend
		plot!(p, nominal_2[:,1] , nominal_2[:,2], 0:H_2, label="Misses", color=davie_green)
		# Plot the good trajectories first
		for j = 1:Int64(round(i))
			if dev_2[j] < pipe_radius
				plot!(trj_2[j,:,1], trj_2[j,:,2], 0:H_2, label="", color=davie_green, opacity=40/n_samples_2)
			end
		end
		
		# Finally, plot the nominal trajectory on top of everything else
		plot!(p, nominal_2[:,1] , nominal_2[:,2], 0:H_2, label="Nominal", color=white, linewidth=3)
		@info i
	end
	webm(anim, fps=30)
end

# ╔═╡ b0a87bbe-b0de-4602-8bf5-b43632c0c068
let
	anim = @animate for i in animsteps(30)
		fi = Int64(floor(i))
		p = plot(
			xlabel=L"x_1", ylabel=L"x_2", zlabel=L"t", guidefontsize=14,
			xlims=(-1.25,0.25), ylims=(-1.25,0.25), zlims=(0,100),
			thickness_scaling=2, size=(1920,1080),
			camera=(azimuth_2,elevation_2),
			legend=:topleft, legendfontsize=14, legendfontfamily="helvetica bold", legend_font_color=white, background_color_legend=navy,
			fontfamily="computer modern",
			background_color=carolina_blue, foreground_color=navy, foreground_color_guide=black, margin=-60px)
		
		# Plot safety pipe
		θ = LinRange(0, 2π, 40)
		circ_x = cos.(θ) * pipe_radius
		circ_y = sin.(θ) * pipe_radius
		for j in 1:Int64(floor(i))
			plot!(circ_x .+ nominal_2[j,1], circ_y .+ nominal_2[j,2], repeat([j-1,], size(θ,1)), label=(j == 1) ? "Safety pipe" : "", seriestype=[:shape,], c=:lightblue, linecolor=basin_slate)
		end

		# Plot random trajectories
		# Fudge the legend
		plot!(p, nominal_2[:,1] , nominal_2[:,2], 0:H_2, label="Misses", color=davie_green)
		# Plot the good trajectories first
		for i = axes(trj_2, 1)
			if dev_2[i] < pipe_radius
				plot!(trj_2[i,:,1], trj_2[i,:,2], 0:H_2, label="", color=davie_green, opacity=40/n_samples_2)
				label_good = true
			end
		end
		
		# Finally, plot the nominal trajectory on top of everything else
		plot!(p, nominal_2[:,1] , nominal_2[:,2], 0:H_2, label="Nominal", color=white, linewidth=3)
		@info i
	end
	webm(anim, fps=30)
end

# ╔═╡ 4ec5556f-ba0e-4e72-aaa8-0b7078ec584d
let
	nom_crosses = zeros(Float64, H_2+1)
	anim = @animate for i in animsteps(30)
		fi = Int64(floor(i))
		p = plot(
			xlabel=L"x_1", ylabel=L"x_2", zlabel=L"t", guidefontsize=14,
			xlims=(-1.25,0.25), ylims=(-1.25,0.25), zlims=(0,100),
			thickness_scaling=2, size=(1920,1080),
			camera=(azimuth_2,elevation_2),
			legend=:topleft, legendfontsize=14, legendfontfamily="helvetica bold", legend_font_color=white, background_color_legend=navy,
			fontfamily="computer modern",
			background_color=carolina_blue, foreground_color=navy, foreground_color_guide=black, margin=-60px)
		
		# Plot safety pipe
		θ = LinRange(0, 2π, 40)
		circ_x = cos.(θ) * pipe_radius
		circ_y = sin.(θ) * pipe_radius
		for j in axes(nominal_2, 1)
			plot!(circ_x .+ nominal_2[j,1], circ_y .+ nominal_2[j,2], repeat([j-1,], size(θ,1)), label=(j == 1) ? "Safety pipe" : "", seriestype=[:shape,], c=:lightblue, linecolor=(nom_crosses[j] > 0) ? navy : basin_slate)
		end

		# Plot random trajectories
		# Fudge the legend
		plot!(p, nominal_2[:,1] , nominal_2[:,2], 0:H_2, label="Misses", color=davie_green)
		# Plot the good trajectories first
		for i = axes(trj_2, 1)
			if dev_2[i] < pipe_radius
				plot!(trj_2[i,:,1], trj_2[i,:,2], 0:H_2, label="", color=davie_green, opacity=40/n_samples_2)
				label_good = true
			end
		end
		
		# Finally, plot the nominal trajectory on top of everything else
		plot!(p, nominal_2[:,1] , nominal_2[:,2], 0:H_2, label="Nominal", color=white, linewidth=3, marker=:xcross, markeralpha=nom_crosses)

		# Plot the bad trajectories on top of the good ones
		label_bad = false
		for j = axes(trj_2, 1)
			if dev_2[j] >= pipe_radius
				itp = interpolate(trj_2[j,:,:], BSpline(Linear()))
				t_viol = argmax(dist_2[j,:], dims=1)
				crosses = [(d >= pipe_radius) ? 1 : 0 for d in dist_2[j,:]]
				plot!(p, [itp[begin:fi,1]..., itp[i,1]] , [itp[begin:fi,2]..., itp[i,2]], [0:fi-1..., i-1], label=label=(!label_bad) ? "Violation" : "", color=azalea_pink, markershape=:xcross, markeralpha=[crosses[begin:fi]..., 0])
				nom_crosses[begin:fi] .+= crosses[begin:fi]
				label_bad = true
			end
		end
		@info i
	end
	webm(anim, fps=30)
end

# ╔═╡ cc9bd794-b5ea-4b1c-8b31-362628922af0
let
	anim = @animate for i in animsteps(30)
		fi = Int64(floor(i))
		p = plot(
			xlabel=L"x_1", ylabel=L"x_2", zlabel=L"t", guidefontsize=14,
			xlims=(-1.25,0.25), ylims=(-1.25,0.25), zlims=(0,100),
			thickness_scaling=2, size=(1920,1080),
			camera=(azimuth_2,elevation_2),
			legend=:topleft, legendfontsize=14, legendfontfamily="helvetica bold", legend_font_color=white, background_color_legend=navy,
			fontfamily="computer modern",
			background_color=carolina_blue, foreground_color=navy, foreground_color_guide=black, margin=-60px)
		
		# Plot the nominal trajectory under everything else
		plot!(p, nominal_2[:,1] , nominal_2[:,2], 0:H_2, label="Nominal", color=white, linewidth=3)

		# Plot the missing trajectories on top of the nominal one
		colors = [longleaf_pine, davie_green, azalea_pink, sunburst_yellow]
		for j = axes(strat_var_2, 1)
			itp = interpolate(strat_var_2[j,:,:], BSpline(Linear()))
			plot!(p, [itp[begin:fi,1]..., itp[i,1]] , [itp[begin:fi,2]..., itp[i,2]], [0:fi-1..., i-1], label=strat_names[j], color=colors[j])
		end
		@info i
	end
	webm(anim, fps=30)
end

# ╔═╡ Cell order:
# ╟─a25ca2c0-b14a-4e26-9de7-4d17b4d40fe8
# ╟─bde0e250-f700-40be-b9dc-eae26763b1d1
# ╟─fe7ab3cf-5bf4-4af8-896d-af1838a089fd
# ╠═6dac77d6-090a-45b6-acdf-0bcda8445f03
# ╟─df7b4e48-122f-4676-be5f-d53e3ed3c354
# ╟─63d725ee-734e-4953-abdc-4ba398534bf6
# ╟─348989fe-4d30-4fbd-a2f3-c0db26b4eed6
# ╟─c1f773c8-aaed-4da8-bb17-132ed5b82ebd
# ╟─24a7cc48-5c21-4e7b-9120-2a9939817bb7
# ╠═af628411-9be0-4ed4-b31f-976d5d110344
# ╟─d623d388-b218-45ef-b1dd-00c6d0ba2f7d
# ╟─4aa065f7-2c8a-4f6e-bd0e-102318d13a99
# ╟─b0a87bbe-b0de-4602-8bf5-b43632c0c068
# ╟─4ec5556f-ba0e-4e72-aaa8-0b7078ec584d
# ╟─cc9bd794-b5ea-4b1c-8b31-362628922af0
# ╟─427002a6-5f53-42b1-be5d-ddef2932e000
# ╟─706d1c80-4b29-4f7a-b85a-7559bc579507
# ╟─c69777bc-a7d0-4ed4-9837-0cf166fc9dd3
# ╟─b730daed-6bb1-47a1-a4b8-ee1e606f3dda
# ╟─53434f2b-1733-4e5b-b174-821c4717fd4d
# ╟─ffa3854f-3633-4b54-a8b7-c449137c54b9
# ╠═46a30bcf-85a0-44a9-9428-fe7206a4f10d
# ╟─1c0d8251-b8a8-4e68-a825-fd8174dc6a04
# ╠═288202fc-c40f-4ddc-ba88-3ded8b3f7dd3
# ╠═11c07e8d-44b0-400f-a1aa-dcb7694cec9b
# ╟─7ac8bae8-feb6-4789-b487-365e970656c7
# ╠═19ac7356-34fc-11ed-06ab-67dbac77f1dc
# ╟─c435311f-72a7-4dc9-915d-ff5112686bde
# ╟─6b58db92-abc1-4123-a805-48d40c0914f0
# ╟─e9b228fa-4674-4e18-a04c-196642246ab3
# ╟─53f962e1-8988-47b5-9e24-ba3b4b273768
# ╟─40f14ec9-cce9-4439-aa20-ead6946acb07
# ╟─801c8332-292d-4c8d-8098-93574cc754b6
# ╟─e9c25a03-2619-4556-a187-d0c02de54f1a
# ╟─fc1073f6-d2c1-4e50-8d9b-1968ac88f1d2
# ╟─319485a7-b5ac-4de0-8b6c-b3c1d4a44cce
# ╟─33a6152a-2fcd-4f88-9d04-1b7074149a34
# ╟─998ce331-09c5-4a6f-abe1-bf61c81bf6e7
# ╟─71032cbe-91e2-4056-808e-4917244c4747
# ╟─73084196-b1f1-4bc7-9190-989b2dcae677

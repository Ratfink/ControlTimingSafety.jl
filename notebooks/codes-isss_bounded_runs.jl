### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ fc76ab0e-be9d-11ed-13cb-6bf8bf993862
begin
    import Pkg
	# Add the ControlTimingSafety.jl dev package 😃
	# This makes Pluto stop managing packages on its own, so you'll have to have all the packages listed below installed in your default environment.
    Pkg.develop(path="../")
    # Instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()
	Pkg.add("CairoMakie")

	using Revise
	using ControlTimingSafety
	using RealTimeScheduling
    using CairoMakie, PlutoUI, LaTeXStrings
	using Random, Distributions
	using Distances
	using ControlSystemsBase, LinearAlgebra
	using Printf
end

# ╔═╡ c0630ced-14fd-46e3-a18f-26840d2db65c
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

# ╔═╡ 0ccdd1fd-4530-4862-938f-acbe71926e2b
ctrl_delay = 0.02

# ╔═╡ 27f322d3-ac0b-47c6-ba0b-f471a5a4a269
sysd = c2d(sys, ctrl_delay)

# ╔═╡ 2d279d84-4c4f-4b2a-b870-29dc3e6b5849
sysd_delay = ss([sysd.A sysd.B; 0 0 0], [0; 0; 1], [sysd.C 0], sysd.D, sysd.Ts)

# ╔═╡ dc5ffeb3-3914-4383-ba25-3126d724bed0
K = let
	Q = [2 0 0;
		 0 2 0;
		 0 0 2]
	R = [1;;]
	lqr(sysd_delay, Q, R)
end

# ╔═╡ 444ec3d5-4d24-4820-9007-cac4cc815583
a = zero_kill(sysd, K)

# ╔═╡ cb335268-1816-405c-96fb-d3127844bd1b
z_0 = augment(a, [1, 1])

# ╔═╡ f7a33c98-9939-470e-a7d4-84480add797b
O = MeetAny(1, 4)

# ╔═╡ 1d28df71-d8af-4d0e-b58b-67330d56ea36
H=7

# ╔═╡ 1aae46ef-7acb-47fc-ab06-43993d2c1901
begin
	f = Figure(resolution=(1800, 800))
	ax = Axis3(f[1, 1], aspect=(3, 4, 1), xticks=0:H, xlabel="t", ylabel="x[1]", zlabel="x[2]", xlabeloffset=30, ylabeloffset=60, zlabeloffset=90, xlabelsize=24, ylabelsize=24, zlabelsize=24, azimuth=1.46π, elevation=π/24)
	evols = []
	for i in 0:2^H-1
		bits = digits(i, base=2, pad=H)
		if BitVector(bits) ⊢ O
			actions = 2 .- bits
			z = evol(a, z_0, actions)
			push!(evols, (bits, z))
			scatterlines!(0:H, z[:,1], z[:,2], color=(:black, 0.05), label=nothing)
		end
	end
	h = []
	hm = []
	hmm = []
	mmm = []
	for (bits, z) in evols
		if bits[end] == 1
			push!(h, z[end,:])
		elseif bits[end-1:end] == [1, 0]
			push!(hm, z[end,:])
		elseif bits[end-2:end] == [1, 0, 0]
			push!(hmm, z[end,:])
		else
			push!(mmm, z[end,:])
		end
	end
	for st in (h, hm, hmm, mmm)
		m = hcat(st...)
		min1, max1 = minimum(m[1,:]), maximum(m[1,:])
		min2, max2 = minimum(m[2,:]), maximum(m[2,:])
		lines!([H, H, H, H, H], [min1, min1, max1, max1, min1], [min2, max2, max2, min2, min2], overdraw=true)
	end
	save("bounded_runs.pdf", f)
	f
end

# ╔═╡ Cell order:
# ╟─c0630ced-14fd-46e3-a18f-26840d2db65c
# ╠═0ccdd1fd-4530-4862-938f-acbe71926e2b
# ╟─27f322d3-ac0b-47c6-ba0b-f471a5a4a269
# ╟─2d279d84-4c4f-4b2a-b870-29dc3e6b5849
# ╟─dc5ffeb3-3914-4383-ba25-3126d724bed0
# ╠═444ec3d5-4d24-4820-9007-cac4cc815583
# ╟─cb335268-1816-405c-96fb-d3127844bd1b
# ╟─f7a33c98-9939-470e-a7d4-84480add797b
# ╠═1d28df71-d8af-4d0e-b58b-67330d56ea36
# ╠═1aae46ef-7acb-47fc-ab06-43993d2c1901
# ╠═fc76ab0e-be9d-11ed-13cb-6bf8bf993862

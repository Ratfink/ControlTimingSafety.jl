using ControlTimingSafety
using ControlSystemsBase

sys = let
	r_1 = 100000
	r_2 = 500000
	r_3 = 200000
	c_1 = 0.000002
	c_2 = 0.000010
	A = [-1/c_1 * (1/r_1 + 1/r_2)  1/(r_2*c_1)
	     1/(r_2*c_2)               -1/c_2 * (1/r_2 + 1/r_3)]
	B = [1/(r_1*c_1)
	     1/(r_3*c_2)]
    C = [1 0; 0 1]
	D = 0

	ss(A, B, C, D)
end
K = [0.16007906919762957 0.21425701514482584 0.022142131526773]

p = 0.023

sysd = c2d(sys, p)

z0 = [100., 100., 0.]

d_max = 1000.
maxwindow = 6
n = 10
H = 100

# display(synthesize_constraints(sysd, K, z0, d_max, maxwindow, n, H, fullresults=true)[2])
for i in 1:10
	println(devub(1, 2, sysd, K, z0, d_max, n, H))
end

using ITensors
using Printf

include("./tpq_mps.jl")

# physical feature
J = 1.0 + 0.0im # Heisenberg model parameter
# J = 4.0 + 0.0im # transverse ising model parameter
# Γ = 2.0 + 0.0im # transverse ising model, transverse field

# mTPQ parameters
l = 1.0 + 0.0im
kmax = 500

# cTPQ parameters
tempstep = 0.05
numtemps = 80

# MPS parameters
n = 4 # system size
χ_aux = 1 # auxiliary site dimension
χ = 20 # bond dimension

# prepare random MPS
ψ, sites = randomTPQMPS(n, χ, χ_aux)

# prepare MPO
l_h = l_heisenberg(sites, n, J, l)
h = heisenberg(sites, n, J)
# l_h = l_transverseising(sites, n, J, Γ, l)
# h = transverseising(sites, n, J, Γ)

# output will be written in "output_ctpq.dat" and "output_mtpq.dat"
canonicalsummation(ψ, l_h, h, χ, l; auxdimension = χ_aux, kmax=kmax, filename="output")
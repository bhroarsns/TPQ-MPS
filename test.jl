using ITensors
using Printf

include("./tpq_mps.jl")

# physical feature
n = 4 # system size
J = 1.0 # Heisenberg model parameter
# J = 4.0
# Γ = 2.0

# MPS parameters
χ = 2 # auxiliary site dimension (for future updates)
ξ = 20 # bond dimension
# mTPQ parameters
l = 1.0
kmax = 500
# cTPQ parameters
tempstep = 0.05
numtemps = 80

# prepare random MPS
ψ, sites = randomTPQMPS(n, χ, ξ)

# prepare MPO
l_h = l_heisenberg(sites, n, J, l)
h = heisenberg(sites, n, J)
# l_h = l_transverseising(sites, n, J, Γ, l)
# h = transverseising(sites, n, J, Γ)
m = magnetization(sites, n)

canonicalsummation(ψ, h, ξ, l)

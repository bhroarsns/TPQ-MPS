using ITensors
using Printf

include("./tpq_mps.jl")

# physical feature
n = 16 # system size
J = 1.0 # Heisenberg model parameter
# J = 4.0
# Γ = 2.0

# MPS parameters
χ = 20 # auxiliary site dimension (for future updates)
ξ = χ # bond dimension
# mTPQ parameters
l = 1.0
kmax = 100
# cTPQ parameters
tempstep = 0.05
numtemps = 80

# prepare random MPS
ψ, sites = randomTPQMPS(n, χ, ξ)

# prepare MPO
h = heisenberg(sites, n, J)
m = magnetization(sites, n)

canonicalsummation(ψ, h, ξ, χ, l, kmax, tempstep, numtemps, m)

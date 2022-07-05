using Distributed
@everywhere using ITensors
@everywhere using Printf

@everywhere include("./tpq_mps.jl")

# physical feature
J = 1.0 + 0.0im # Heisenberg model parameter
# J = 4.0
# Γ = 2.0

# mTPQ parameters
l = 1.0 + 0.0im
kmax = 500
# cTPQ parameters
tempstep = 0.05
numtemps = 80

# # MPS parameters
# n = 4 # system size
# χ = 20 # bond dimension
# χ_aux = 20

# ψ, sites = randomTPQMPS(n, χ, χ_aux)

# MPS parameters
n = 16 # system size
χ = 20 # bond dimension

χ_aux = 20

future1 = @spawnat :any for i in 1:10
    # prepare random MPS
    local ψ, sites = randomTPQMPS(n, χ, χ_aux)
    
    # prepare MPO
    local l_h = l_heisenberg(sites, n, J, l)
    local h = heisenberg(sites, n, J)
    # l_h = l_transverseising(sites, n, J, Γ, l)
    # h = transverseising(sites, n, J, Γ)

    canonicalsummation(ψ, l_h, h, χ, l; auxdimension = χ_aux, kmax=kmax, filename="result/heisen-$χ_aux/run$i")
end

future2 = @spawnat :any for i in 26:35
    # prepare random MPS
    local ψ, sites = randomTPQMPS(n, χ, χ_aux)
    
    # prepare MPO
    local l_h = l_heisenberg(sites, n, J, l)
    local h = heisenberg(sites, n, J)
    # l_h = l_transverseising(sites, n, J, Γ, l)
    # h = transverseising(sites, n, J, Γ)

    canonicalsummation(ψ, l_h, h, χ, l; auxdimension = χ_aux, kmax=kmax, filename="result/heisen-$χ_aux/run$i")
end

wait(future1)
wait(future2)

future1 = @spawnat :any for i in 16:25
    # prepare random MPS
    local ψ, sites = randomTPQMPS(n, χ, χ_aux)
    
    # prepare MPO
    local l_h = l_heisenberg(sites, n, J, l)
    local h = heisenberg(sites, n, J)
    # l_h = l_transverseising(sites, n, J, Γ, l)
    # h = transverseising(sites, n, J, Γ)

    canonicalsummation(ψ, l_h, h, χ, l; auxdimension = χ_aux, kmax=kmax, filename="result/heisen-$χ_aux/run$i")
end

future2 = @spawnat :any for i in 41:50
    # prepare random MPS
    local ψ, sites = randomTPQMPS(n, χ, χ_aux)
    
    # prepare MPO
    local l_h = l_heisenberg(sites, n, J, l)
    local h = heisenberg(sites, n, J)
    # l_h = l_transverseising(sites, n, J, Γ, l)
    # h = transverseising(sites, n, J, Γ)

    canonicalsummation(ψ, l_h, h, χ, l; auxdimension = χ_aux, kmax=kmax, filename="result/heisen-$χ_aux/run$i")
end

wait(future1)
wait(future2)

χ_aux = 40

future1 = @spawnat :any for i in 1:25
    # prepare random MPS
    local ψ, sites = randomTPQMPS(n, χ, χ_aux)
    
    # prepare MPO
    local l_h = l_heisenberg(sites, n, J, l)
    local h = heisenberg(sites, n, J)
    # l_h = l_transverseising(sites, n, J, Γ, l)
    # h = transverseising(sites, n, J, Γ)

    canonicalsummation(ψ, l_h, h, χ, l; auxdimension = χ_aux, kmax=kmax, filename="result/heisen-$χ_aux/run$i")
end

future2 = @spawnat :any for i in 26:50
    # prepare random MPS
    local ψ, sites = randomTPQMPS(n, χ, χ_aux)
    
    # prepare MPO
    local l_h = l_heisenberg(sites, n, J, l)
    local h = heisenberg(sites, n, J)
    # l_h = l_transverseising(sites, n, J, Γ, l)
    # h = transverseising(sites, n, J, Γ)

    canonicalsummation(ψ, l_h, h, χ, l; auxdimension = χ_aux, kmax=kmax, filename="result/heisen-$χ_aux/run$i")
end

wait(future1)
wait(future2)
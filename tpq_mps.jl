export l_heisenberg, l_transverseising, heisenberg, transverseising, magnetization, randomTPQMPS, canonicalform!, canonicalsummation

using ITensors
using Printf
using Dates

function l_heisenberg(sites::Vector{Index{Int64}}, n::Int64, J::Float64, l::Float64)
    JJ = J / n

    mpo = OpSum()
    mpo += l, "I", 2
    for j = 2:n
        mpo += -JJ / 2.0, "S+", j, "S-", j+1
        mpo += -JJ / 2.0, "S-", j, "S+ ", j+1
        mpo += -JJ, "Sz", j, "Sz", j+1
    end

    return MPO(mpo, sites)
end

function l_transverseising(sites::Vector{Index{Int64}}, n::Int64, J::Float64, Γ::Float64, l::Float64)
    JJ = J / n
    ΓΓ = Γ / n

    mpo = OpSum()
    mpo += l, "I", 2
    for j = 2:n
        mpo += -JJ, "Sz", j, "Sz", j+1
        mpo += -ΓΓ, "Sx", j
    end
    mpo += -ΓΓ, "Sx", n + 1

    return MPO(mpo, sites)
end

function heisenberg(sites::Vector{Index{Int64}}, n::Int64, J::Float64)
    mpo = OpSum()
    for j = 2:n
        mpo += J / 2.0, "S+", j, "S-", j+1
        mpo += J / 2.0, "S-", j, "S+ ", j+1
        mpo += J, "Sz", j, "Sz", j+1
    end

    return MPO(mpo, sites)
end

function transverseising(sites::Vector{Index{Int64}}, n::Int64, J::Float64, Γ::Float64)
    mpo = OpSum()
    for j = 2:n
        mpo += J, "Sz", j, "Sz", j+1
        mpo += Γ, "Sx", j
    end
    mpo += Γ, "Sx", n + 1

    return MPO(mpo, sites)
end

function magnetization(sites::Vector{Index{Int64}}, n::Int64)
    mpo = OpSum()
    for j = 2:n+1
        mpo += "Sz", j
    end

    return MPO(mpo, sites)
end

function randomTPQMPS(system_size::Int64, bonddimension::Int64, auxdimension::Int64, sitetype::String="S=1/2")
    sites = siteinds(sitetype, system_size + 2)
    sites[1] = Index(auxdimension, tags="Link,l=1")
    sites[system_size + 2] = Index(auxdimension, tags="Link,l=$(system_size+1)")

    ψ = randomMPS(Complex{Float64}, sites, bonddimension)

    # replace #1 site by auxiliary site
    auxL = Index(auxdimension, tags="Site,aux-L")
    ψ[2] = ψ[1] * ψ[2]
    ψ[1] = delta(Complex{Float64}, auxL, sites[1])
    sites[1] = auxL

    # replace #n+1 site by auxiliary site
    auxR = Index(auxdimension, tags="Site,aux-R")
    ψ[system_size + 1] = ψ[system_size + 1] * ψ[system_size + 2]
    ψ[system_size + 2] = delta(Complex{Float64}, sites[system_size + 2], auxR)
    sites[system_size + 2] = auxR

    println("initial state successfully generated")
    return (ψ, sites)
end

function canonicalform!(A::MPS, bonddimension::Int64, auxdimension::Int64)
    len = length(A)
    sites = siteinds(A)

    if len != 1
        U1, P1 = factorize(A[1] * A[2], sites[1]; which_decomp="qr")
        prevlink, _ = inds(P1)
        A[1] = U1
        A[2] = P1

        for j = 2:len-1
            Uj, Pj = factorize(A[j] * A[j+1], prevlink, sites[j]; which_decomp="qr")
            prevlink, _ = inds(Pj)
            A[j] = Uj
            A[j+1] = Pj
        end

        Vn, Pn = factorize(A[len-1] * A[len], sites[len]; which_decomp="svd", maxdim=auxdimension)
        prevlink, _ = inds(Pn)
        A[len] = Vn
        A[len-1] = Pn

        for j = len-1:-1:3
            Vj, Pj = factorize(A[j - 1] * A[j], sites[j], prevlink; which_decomp="svd", maxdim = bonddimension)
            prevlink, _ = inds(Pj)
            A[j] = Vj
            A[j - 1] = Pj
        end
    end
end

function canonicalsummation(initial::MPS, l_h::MPO, hamiltonian::MPO, bonddimension::Int64, l::Float64; auxdimension::Int64=bonddimension, kmax::Int64=500, tempstep::Float64=0.05, numtemps::Int64=80, filename::String="output")
    println("")
    println("=====================================")
    println("$(Dates.now()): start calculating canonical summation")
    ψ = deepcopy(initial)
    n = length(ψ) - 2

    # operator preparation
    m = magnetization(siteinds(ψ), n)
    h2 = hamiltonian' * hamiltonian
    m2 = m' * m

    # data storage
    energys = zeros(Complex{Float64}, numtemps, 2kmax + 2)
    energysquares = zeros(Complex{Float64}, numtemps, 2kmax + 2)
    magnets = zeros(Complex{Float64}, numtemps, 2kmax + 2)
    magnetsquares = zeros(Complex{Float64}, numtemps, 2kmax + 2)
    betanorms = zeros(Complex{Float64}, numtemps, 2kmax + 2)
    factors = zeros(Complex{Float64}, numtemps, kmax + 2)
    mtpqdata = zeros(Complex{Float64}, kmax+1, 2)

    # temperatures to calculate cTPQ
    temps = (numtemps * tempstep):-tempstep:tempstep
    println("temperature range: from $(tempstep) to $(numtemps * tempstep)")

    for k = 0:kmax
        @printf "\rstep %03i" k
        
        khk = log(inner(ψ', hamiltonian, ψ))
        kh2k = log(inner(ψ'', h2, ψ))
        kmk = log(inner(ψ', m, ψ))
        km2k = log(inner(ψ'', m2, ψ))

        mtpqdata[k+1, 1] = exp(khk)
        mtpqdata[k+1, 2] = exp(kh2k)
        
        ϕ = apply(l_h, ψ)
        canonicalform!(ϕ, bonddimension, auxdimension)

        khk1 = log(inner(ψ', hamiltonian, ϕ))
        kh2k1 = log(inner(ψ'', h2, ϕ))
        kmk1 = log(inner(ψ', m, ϕ))
        km2k1 = log(inner(ψ'', m2, ϕ))
        kk1 = logdot(ψ, ϕ)

        ψ = deepcopy(ϕ)
        kk = lognorm(ψ)
        normalize!(ψ)

        for (i, t) in enumerate(temps)
            factor = factors[i, k+1]
            new_factor = factor + log(n) - log(t) - log(2k + 1)
            energys[i, 2k+1] = khk + factor
            energys[i, 2k+2] = khk1 + new_factor
            energysquares[i, 2k+1] = kh2k + factor
            energysquares[i, 2k+2] = kh2k1 + new_factor
            magnets[i, 2k+1] = kmk + factor
            magnets[i, 2k+2] = kmk1 + new_factor
            magnetsquares[i, 2k+1] = km2k + factor
            magnetsquares[i, 2k+2] = km2k1 + new_factor
            betanorms[i, 2k+1] = factor
            betanorms[i, 2k+2] = kk1 + new_factor
            factors[i, k+2] = new_factor + log(n) - log(t) - log(2k + 2) + 2.0 * kk
        end
    end
    println("\r$(Dates.now()): canonical summation finished")

    open(string(filename, "_mtpq.dat"), "w") do io
        write(io, "# temperature, energy density, specific heat, energy variance, energy, energy^2\n")
        for k in 1:kmax+1
            temp = real(n*(l - mtpqdata[k, 1])/2(k-1))
            energy = real(mtpqdata[k, 1])
            energysquare = real(mtpqdata[k, 2])
            energyvariance = energysquare - energy^2

            write(io, "$temp $(energy/n) $(energyvariance / temp^2 /n) $energyvariance $energy $energysquare\n")
        end
        println("output written in \"$(filename)_mtpq.dat\"")
    end

    open(string(filename, "_ctpq.dat"), "w") do io
        write(io, "# temperature, energy density, mean magnetization, specific heat, susceptibility, norm^2\n")
        for (i, t) in enumerate(temps)
            sortedenergys = sort(energys[i, :]; by = x -> real(x))
            sortedenergysquares = sort(energysquares[i, :]; by = x -> real(x))
            sortedmagnets = sort(magnets[i, :]; by = x -> real(x))
            sortedmagnetsquares = sort(magnetsquares[i, :]; by = x -> real(x))
            sortedbetanorms = sort(betanorms[i, :]; by = x -> real(x))

            divider = real(last(sortedbetanorms))
        
            ene = 0.0
            for e in sortedenergys
                ene += exp(e - divider)
            end

            ene2 = 0.0
            for e2 in sortedenergysquares
                ene2 += exp(e2 - divider)
            end

            mag = 0.0
            for m in sortedmagnets
                mag += exp(m - divider)
            end

            mag2 = 0.0
            for m2 in sortedmagnetsquares
                mag2 += exp(m2 - divider)
            end
        
            bet_nor = 0.0
            for b in sortedbetanorms
                bet_nor += exp(b - divider)
            end
            
            write(io, "$t $(real(ene) / n) $(real(mag) / n) $(real(ene2) / n) $(real(mag2) / n) $(real(bet_nor))\n")
        end
        println("output written in \"$(filename)_ctpq.dat\"")
    end
    println("=====================================")
    println("")
end

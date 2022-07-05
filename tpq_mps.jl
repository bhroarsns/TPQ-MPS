using ITensors
using Printf
using Dates

function l_heisenberg(sites::Vector{Index{Int64}}, n::Int64, J::Complex{Float64}, l::Complex{Float64})
    JJ = J / n

    mpo = OpSum()
    mpo += l, "I", 2
    for j = 2:n
        mpo += -JJ / (2.0 + 0.0im), "S+", j, "S-", j+1
        mpo += -JJ / (2.0 + 0.0im), "S-", j, "S+ ", j+1
        mpo += -JJ, "Sz", j, "Sz", j+1
    end

    return MPO(mpo, sites)
end

function l_transverseising(sites::Vector{Index{Int64}}, n::Int64, J::Complex{Float64}, Γ::Complex{Float64}, l::Complex{Float64})
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

function heisenberg(sites::Vector{Index{Int64}}, n::Int64, J::Complex{Float64})
    mpo = OpSum()
    for j = 2:n
        mpo += J / (2.0 + 0.0im), "S+", j, "S-", j+1
        mpo += J / (2.0 + 0.0im), "S-", j, "S+ ", j+1
        mpo += J, "Sz", j, "Sz", j+1
    end

    return MPO(mpo, sites)
end

function transverseising(sites::Vector{Index{Int64}}, n::Int64, J::Complex{Float64}, Γ::Complex{Float64})
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

    if len != 1
        for i = 1:len-2
            Qi, Ri = factorize(A[i] * A[i+1], uniqueinds(A[i], A[i+1]); which_decomp="qr", tags="Link,l=1")
            A[i] = Qi
            A[i+1] = Ri
        end

        Vn, snUn = factorize(A[len] * A[len-1], uniqueinds(A[len], A[len-1]); which_decomp="svd", maxdim=auxdimension, tags="Link,l=$(len-1)")
        A[len] = Vn
        A[len-1] = snUn

        for i = len-1:-1:3
            Vi, siUi = factorize(A[i] * A[i-1], uniqueinds(A[i], A[i-1]); which_decomp="svd", maxdim=bonddimension, tags="Link,l=$(i-1)")
            A[i] = Vi
            A[i-1] = siUi
        end
    end
end

function canonicalsummation(initial::MPS, l_h::MPO, hamiltonian::MPO, bonddimension::Int64, l::Complex{Float64}; auxdimension::Int64=bonddimension, kmax::Int64=500, tempstep::Float64=0.05, numtemps::Int64=80, filename::String="output")
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
    mtpqdata = zeros(Complex{Float64}, kmax+1, 4)

    # temperatures to calculate cTPQ
    temps = (numtemps * tempstep):-tempstep:tempstep
    println("temperature range: from $(tempstep) to $(numtemps * tempstep)")

    for k = 0:kmax
        @printf "\rstep %03i\n" k
        
        khk = log(inner(ψ', hamiltonian, ψ))
        kh2k = log(inner(ψ'', h2, ψ))
        kmk = log(inner(ψ', m, ψ))
        km2k = log(inner(ψ'', m2, ψ))

        mtpqdata[k+1, 1] = exp(khk)
        mtpqdata[k+1, 2] = exp(kh2k)
        mtpqdata[k+1, 3] = exp(kmk)
        mtpqdata[k+1, 4] = exp(km2k)
        
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
        println(io, "# k_BT, h, m, C/N, chi, delta E, delta M")
        for k in 1:kmax+1
            temp = real(n*(l - mtpqdata[k, 1])/2(k-1))
            energy = real(mtpqdata[k, 1])
            energysquare = real(mtpqdata[k, 2])
            energyvariance = energysquare - energy^2

            magnet = real(mtpqdata[k, 3])
            magnetsquare = real(mtpqdata[k, 4])
            magnetvariance = magnetsquare - magnet^2

            println(io, temp, " ", energy/n, " ", magnet/n, " ", energyvariance/temp^2/n, " ", magnetvariance/temp/n, " ", energyvariance, " ", magnetvariance)
        end
        println("output written in \"$(filename)_mtpq.dat\"")
    end

    open(string(filename, "_ctpq.dat"), "w") do io
        println(io, "# k_BT, <a|H|a>, <a|S_z|a>, <a|H^2|a>, <a|S_z^2|a>, <a|a>, N")
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
            
            println(io, t, " ", real(ene), " ", real(mag), " ", real(ene2), " ", real(mag2), " ", real(bet_nor), " ", n)
        end
        println("output written in \"$(filename)_ctpq.dat\"")
    end
    println("=====================================")
    println("")
end

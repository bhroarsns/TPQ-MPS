export heisenberg, transverseising, magnetization, randomTPQMPS, canonicalform!, canonicalsummation

using ITensors
using Printf

function l_heisenberg(sites::Vector{Index{Int64}}, n::Int64, J::Float64, l::Float64)
    JJ = J / n

    mpo = OpSum()
    mpo += l, "I", 2
    for j = 2:n
        mpo += -JJ / 2.0, "S+", j, "S-", j+1
        mpo += -JJ / 2.0, "S-", j, "S+ ", j+1
        mpo += -JJ, "Sz", j, "Sz", j+1
    end

    println("iteration operator of the Heisenberg model generated")
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

    println("iteration operator of the transverse Ising model generated")
    return MPO(mpo, sites)
end

function heisenberg(sites::Vector{Index{Int64}}, n::Int64, J::Float64)
    mpo = OpSum()
    for j = 2:n
        mpo += J / 2.0, "S+", j, "S-", j+1
        mpo += J / 2.0, "S-", j, "S+ ", j+1
        mpo += J, "Sz", j, "Sz", j+1
    end

    println("MPO of the energy density operator of the Heisenberg model generated")
    return MPO(mpo, sites)
end

function transverseising(sites::Vector{Index{Int64}}, n::Int64, J::Float64, Γ::Float64)
    mpo = OpSum()
    for j = 2:n
        mpo += J, "Sz", j, "Sz", j+1
        mpo += Γ, "Sx", j
    end
    mpo += Γ, "Sx", n + 1

    println("MPO of the energy density operator of the transverse Ising model generated")
    return MPO(mpo, sites)
end

function magnetization(sites::Vector{Index{Int64}}, n::Int64)
    mpo = OpSum()
    for j = 2:n+1
        mpo += "Sz", j
    end

    println("MPO of magnetization generated")
    return MPO(mpo, sites)
end

function randomTPQMPS(system_size::Int64, χ::Int64, ξ::Int64, sitetype::String="S=1/2")
    sites = siteinds(sitetype, system_size + 2)
    sites[1] = Index(χ, tags="Link,l=1")
    sites[system_size + 2] = Index(χ, tags="Link,l=$(system_size+1)")

    ψ = randomMPS(Complex{Float64}, sites, ξ)

    # replace #1 site by auxiliary site
    auxL = Index(χ, tags="Site,aux-L")
    ψ[2] = ψ[1] * ψ[2]
    ψ[1] = delta(Complex{Float64}, auxL, sites[1])
    sites[1] = auxL

    # replace #n+1 site by auxiliary site
    auxR = Index(χ, tags="Site,aux-R")
    ψ[system_size + 1] = ψ[system_size + 1] * ψ[system_size + 2]
    ψ[system_size + 2] = delta(Complex{Float64}, sites[system_size + 2], auxR)
    sites[system_size + 2] = auxR

    println("initial state successfully generated")
    return (ψ, sites)
end

function canonicalform!(A::MPS, ξ::Int64, χ::Int64)
    len = length(A)
    sites = siteinds(A)

    if len != 1
        U1, s1, V1 = LinearAlgebra.svd(A[1] * A[2], sites[1]; lefttags = "Link,l=1")
        prevlink, _ = inds(s1)
        A[1] = U1
        A[2] = s1 * V1

        for j = 2:len-1
            Uj, sj, Vj = LinearAlgebra.svd(A[j] * A[j + 1], sites[j], prevlink)
            prevlink, _ = inds(sj)
            A[j] = Uj
            A[j + 1] = sj * Vj
        end

        Vn, sn, Un = LinearAlgebra.svd(A[len - 1] * A[len], sites[len]; lefttags = "Link,l=$(len-1)", maxdim = χ)
        prevlink, _ = inds(sn)
        A[len] = Vn
        A[len - 1] = Un * sn

        for j = len-1:-1:3
            Vj, sj, Uj = LinearAlgebra.svd(A[j - 1] * A[j], sites[j], prevlink; lefttags = "Link,l=$(j-1)", maxdim = ξ)
            prevlink, _ = inds(sj)
            A[j] = Vj
            A[j - 1] = Uj * sj
        end
    end
end

function canonicalsummation(initial::MPS, l_h::MPO, h::MPO, ξ::Int64, l::Float64; χ::Int64=ξ, kmax::Int64=500, tempstep::Float64=0.05, numtemps::Int64=80, filename::String="output")
    println("")
    println("=====================================")
    println("start calculating canonical summation")
    ψ = deepcopy(initial)
    n = length(ψ) - 2

    # operator preparation
    m = magnetization(siteinds(ψ), n)
    h2 = h' * h
    m2 = m' * m
    # h2 = contract(h', h, truncate=false)
    # m2 = contract(m', m, truncate=false)

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
        
        khk = log(inner(ψ', h, ψ))
        kh2k = log(inner(ψ'', h2, ψ))
        kmk = log(inner(ψ', m, ψ))
        km2k = log(inner(ψ'', m2, ψ))

        mtpqdata[k+1, 1] = exp(khk)
        mtpqdata[k+1, 2] = exp(kh2k)
        
        ϕ = apply(l_h, ψ#=, truncate=false=#)
        canonicalform!(ϕ, ξ, χ)

        khk1 = log(inner(ψ', h, ϕ))
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
    println("\rcanonical summation finished")

    open(string(filename, "_mtpq.dat"), "w") do io
        for k in 1:kmax+1
            temp = real(n*(l - mtpqdata[k, 1])/2(k-1))
            energy = real(mtpqdata[k, 1])
            energysquare = real(mtpqdata[k, 2])
            energyvariance = energysquare - energy^2

            write(io, "$temp\t$(energy/n)\t$(energyvariance / temp^2 /n)\t$energyvariance\t$energy\t$energysquare\n")
        end
        println("output written in \"$(filename)_mtpq.dat\"")
    end

    open(string(filename, "_ctpq.dat"), "w") do io
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

            ene = ene / bet_nor
            mag = mag / bet_nor
            dev_ene = ene2 / bet_nor - ene^2
            dev_mag = mag2 / bet_nor - mag^2
            
            write(io, "$t\t$(real(ene) / n)\t$(real(mag) / n)\t$(real(dev_ene) / t^2 / n)\t$(real(dev_mag) / t / n)\n")
        end
        println("output written in \"$(filename)_ctpq.dat\"")
    end
    println("=====================================")
    println("")
end

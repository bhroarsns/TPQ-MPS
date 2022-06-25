export heisenberg, transverseising, magnetization, randomTPQMPS, canonicalform!, canonicalsummation

using ITensors
using Printf

function l_heisenberg(sites::Vector{Index{Int64}}, n::Int64, J::Float64, l::Float64)
    JJ = J / n
    mpo = OpSum()
    mpo += l, "I", 1
    for j = 2:n
        mpo += -JJ / 2.0, "S+", j, "S-", j+1
        mpo += -JJ / 2.0, "S-", j, "S+ ", j+1
        mpo += -JJ, "Sz", j, "Sz", j+1
    end
    println("MPO of the iteration operator of the Heisenberg model generated")
    return MPO(mpo, sites)
end

function l_transverseising(sites::Vector{Index{Int64}}, n::Int64, J::Float64, Γ::Float64, l::Float64)
    JJ = J / n
    ΓΓ = Γ / n
    mpo = OpSum()
    mpo += l, "I", 1
    for j = 2:n
        mpo += -JJ, "Sz", j, "Sz", j+1
        mpo += -ΓΓ, "Sx", j
    end
    mpo += ΓΓ, "Sx", n + 1
    println("MPO of the iteration operator of the transverse Ising model generated")
    return MPO(mpo, sites)
end

function heisenberg(sites::Vector{Index{Int64}}, n::Int64, J::Float64)
    JJ = J
    mpo = OpSum()
    for j = 2:n
        mpo += JJ / 2.0, "S+", j, "S-", j+1
        mpo += JJ / 2.0, "S-", j, "S+ ", j+1
        mpo += JJ, "Sz", j, "Sz", j+1
    end
    println("MPO of the energy density operator of the Heisenberg model generated")
    return MPO(mpo, sites)
end

function transverseising(sites::Vector{Index{Int64}}, n::Int64, J::Float64, Γ::Float64)
    JJ = J
    ΓΓ = Γ
    mpo = OpSum()
    for j = 2:n
        mpo += JJ, "Sz", j, "Sz", j+1
        mpo += ΓΓ, "Sx", j
    end
    mpo += ΓΓ, "Sx", n + 1
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
    sites = siteinds(sitetype, system_size)
    sizehint!(sites, system_size + 2)
    pushfirst!(sites, Index(χ, tags="Site,aux-L"))
    push!(sites, Index(χ, tags="Site,aux-R"))
    ψ = randomMPS(Complex{Float64}, sites, ξ)
    ψ[1] = delta(Complex{Float64}, sites[1], linkinds(ψ, 1))
    ψ[system_size + 2] = delta(Complex{Float64}, linkinds(ψ, system_size + 1), sites[system_size + 2])
    println("initial state successfully generated: norm = $(norm(ψ)) vs sqrt(χ) = $(sqrt(χ))")

    return (ψ, sites)
end

function canonicalform!(A::MPS, ξ::Int64, χ::Int64)
    len = length(A)
    sites = siteinds(A)

    if len != 1
        U1, s1, V1 = LinearAlgebra.svd(A[1] * A[2], sites[1])
        prev_link, _ = inds(s1)
        A[1] = U1
        A[2] = s1 * V1

        for j = 2:len-1
            Uj, sj, Vj = LinearAlgebra.svd(A[j] * A[j + 1], sites[j], prev_link)
            prev_link, _ = inds(sj)
            A[j] = Uj
            A[j + 1] = sj * Vj
        end

        Vn, sn, Un = LinearAlgebra.svd(A[len - 1] * A[len], sites[len]; lefttags = "Link,l=$(len-1)", maxdim = χ)
        prev_link, _ = inds(sn)
        A[len] = Vn
        A[len - 1] = Un * sn

        for j = len-1:-1:3
            Vj, sj, Uj = LinearAlgebra.svd(A[j - 1] * A[j], sites[j], prev_link; lefttags = "Link,l=$(j-1)", maxdim = ξ)
            prev_link, _ = inds(sj)
            A[j] = Vj
            A[j - 1] = Uj * sj
        end
    end
end

function canonicalsummation(init::MPS, h::MPO, ξ::Int64, χ::Int64, l::Float64, kmax::Int64, tempstep::Float64, numtemps::Int64, m::MPO, filename::String="output")
    println("start calculating canonical summation")

    n = length(init) - 2

    energy = zeros(Complex{Float64}, numtemps, 2kmax + 2)
    energy2 = zeros(Complex{Float64}, numtemps, 2kmax + 2)
    magnet = zeros(Complex{Float64}, numtemps, 2kmax + 2)
    magnet2 = zeros(Complex{Float64}, numtemps, 2kmax + 2)
    beta_norm = zeros(Complex{Float64}, numtemps, 2kmax + 2)
    factors = zeros(Complex{Float64}, numtemps, kmax + 2)

    # temperatures to calculate cTPQ
    temps = (numtemps * tempstep):-tempstep:tempstep
    println("temperature range: from $(tempstep) to $(numtemps * tempstep)")

    ψ = deepcopy(init)
    normalize!(ψ)

    h2 = h' * h
    m2 = m' * m

    mtpq_data = zeros(Complex{Float64}, kmax+1, 2)

    for k = 0:kmax
        @printf "\rstep %03i" k

        khk = inner(ψ', h, ψ)
        kh2k = inner(ψ'', h2, ψ)
        kmk = log(inner(ψ', m, ψ))
        km2k = log(inner(ψ'', m2, ψ))

        mtpq_data[k+1, 1] = khk
        mtpq_data[k+1, 2] = kh2k

        ϕ = l * ψ - apply(h, ψ)
        # @printf "\tmTPQ state(k = %03i) generated" k+1
        canonicalform!(ϕ, ξ, χ)
        # print("\tcanonicalization finished")

        khk1 = log(l * khk - kh2k) # ⟨k|h|k+1⟩ = ⟨k|h(l - h)|k⟩ = l⟨k|h|k⟩ - ⟨k|h²|k⟩
        kh2k = log(kh2k)
        kh2k1 = log(inner(ψ'', h2, ϕ))
        kmk1 = log(inner(ψ', m, ϕ))
        km2k1 = log(inner(ψ'', m2, ϕ))
        kk1 = log(l - khk) # ⟨k|k+1⟩ = ⟨k|(l - h)|k⟩ = l⟨k|k⟩ - ⟨k|h|k⟩ = l⟨k|k⟩ - ⟨k|h|k⟩
        khk = log(khk)

        ψ = deepcopy(ϕ)
        norm_k = []
        normalize!(ψ; (lognorm!)=norm_k)

        for (i, t) in enumerate(temps)
            factor = factors[i, k+1]
            new_factor = factor + log(n) - log(t) - log(2k + 1)
            energy[i, 2k+1] = khk + factor
            energy[i, 2k+2] = khk1 + new_factor
            energy2[i, 2k+1] = kh2k + factor
            energy2[i, 2k+2] = kh2k1 + new_factor
            magnet[i, 2k+1] = kmk + factor
            magnet[i, 2k+2] = kmk1 + new_factor
            magnet2[i, 2k+1] = km2k + factor
            magnet2[i, 2k+2] = km2k1 + new_factor
            beta_norm[i, 2k+1] = factor
            beta_norm[i, 2k+2] = kk1 + new_factor
            factors[i, k+2] = new_factor + log(n) - log(t) - log(2k + 2) + 2 * norm_k[1]
        end
    end

    open(string(filename, "_mtpq.dat"), "w") do io
        for k in 1:kmax+1
            write(io, "$(real(n*(l - mtpq_data[k, 1])/2(k-1)))\t$(real(mtpq_data[k, 1]))\t$(real(mtpq_data[k, 2] - mtpq_data[k, 1]^2))\n")
        end
    end

    println("\ncanonical summation finished")

    open(string(filename, "_ctpq.dat"), "w") do io
        for (i, t) in enumerate(temps)
            energyt = sort(energy[i, :]; by = x -> real(x))
            energy2t = sort(energy2[i, :]; by = x -> real(x))
            magnett = sort(magnet[i, :]; by = x -> real(x))
            magnet2t = sort(magnet2[i, :]; by = x -> real(x))
            beta_normt = sort(beta_norm[i, :]; by = x -> real(x))

            divider = real(last(beta_normt))
        
            ene = 0.0
            for e in energyt
                ene += exp(e - divider)
            end

            ene2 = 0.0
            for e in energy2t
                ene2 += exp(e - divider)
            end

            mag = 0.0
            for e in magnett
                mag += exp(e - divider)
            end

            mag2 = 0.0
            for e in magnet2t
                mag2 += exp(e - divider)
            end
        
            bet_nor = 0.0
            for b in beta_normt
                bet_nor += exp(b - divider)
            end

            ene = ene / bet_nor
            mag = mag / bet_nor

            dev_ene = ene2 / bet_nor - ene^2
            dev_mag = mag2 / bet_nor - mag^2
            
            write(io, "$t\t$(real(ene))\t$(real(mag))\t$(n*real(dev_ene)/t^2)\t$(n*real(dev_mag)/t)\n")
        end
    end
end
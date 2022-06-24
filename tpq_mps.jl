export heisenberg, transverseising, magnetization, randomTPQMPS, canonicalform!, canonicalsummation

using ITensors
using Printf

function heisenberg(sites, n, J)
    JJ = J / n
    mpo = OpSum()
    for j = 2:n
        mpo += JJ / 2.0, "S+", j, "S-", j+1
        mpo += JJ / 2.0, "S-", j, "S+ ", j+1
        mpo += JJ, "Sz", j, "Sz", j+1
    end
    println("MPO of the energy density operator of the Heisenberg model generated")
    return MPO(mpo, sites)
end

function transverseising(sites, n, J, Γ)
    JJ = J / n
    ΓΓ = Γ / n
    mpo = OpSum()
    for j = 2:n
        mpo += JJ, "Sz", j, "Sz", j+1
        mpo += ΓΓ, "Sx", j
    end
    mpo += ΓΓ, "Sx", n + 1
    println("MPO of the energy density operator of the transverse Ising model generated")
    return MPO(mpo, sites)
end

function magnetization(sites, n)
    mpo = OpSum()
    for j = 2:n+1
        mpo += 1/n,"Sz", j
    end
    println("MPO of magnetization generated")
    return MPO(mpo, sites)
end

function randomTPQMPS(system_size, χ, ξ=χ, sitetype="S=1/2")
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

function canonicalform!(A::MPS, ξ, χ)
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

function canonicalsummation(init::MPS, h::MPO, ξ, χ, l, kmax, tempstep, numtemps, m)
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
    println("temperature range: $(tempstep) to $(numtemps * tempstep)")
    
    ψ = deepcopy(init)
    normalize!(ψ)

    for k = 0:kmax
        @printf "\rstep %03i\t" k

        global khk = inner(ψ', h, ψ)
        global kh2k = inner(h, ψ, h, ψ)
        global kmk = log(inner(ψ', m, ψ))
        global km2k = log(inner(m, ψ, m, ψ))

        global ϕ = l * ψ - apply(h, ψ)
        @printf "mTPQ state(k = %03i) generated" k+1
        canonicalform!(ϕ, ξ, χ)

        global khk1 = log(l * khk - kh2k) # ⟨k|h|k+1⟩ = ⟨k|h(l - h)|k⟩ = l⟨k|h|k⟩ - ⟨k|h²|k⟩
        global kh2k = log(kh2k)
        global kh2k1 = log(inner(h, ψ, h, ϕ))
        global kmk1 = log(inner(ψ', m, ϕ))
        global km2k1 = log(inner(m, ψ, m, ϕ))
        global kk1 = log(l - khk) # ⟨k|k+1⟩ = ⟨k|(l - h)|k⟩ = l⟨k|k⟩ - ⟨k|h|k⟩ = l⟨k|k⟩ - ⟨k|h|k⟩
        global khk = log(khk)

        global ψ = deepcopy(ϕ)
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
            factors[i, k+2] = new_factor + log(n) - log(t) - log(2k + 2) + norm_k[1]
        end
    end
    print("\ncanonical summation finished\n")

    open("output.dat", "w") do io
        for (i, t) in enumerate(temps)
            sort!(energy[i, :]; by = x -> real(x))
            sort!(energy2[i, :]; by = x -> real(x))
            sort!(magnet[i, :]; by = x -> real(x))
            sort!(magnet2[i, :]; by = x -> real(x))
            sort!(beta_norm[i, :]; by = x -> real(x))
        
            divider = real(last(beta_norm[i, :]))
        
            ene = 0.0
            for e in energy[i, :]
                ene += exp(e - divider)
            end

            ene2 = 0.0
            for e in energy2[i, :]
                ene2 += exp(e - divider)
            end

            mag = 0.0
            for e in magnet[i, :]
                mag += exp(e - divider)
            end

            mag2 = 0.0
            for e in magnet2[i, :]
                mag2 += exp(e - divider)
            end
        
            bet_nor = 0.0
            for b in beta_norm[i, :]
                bet_nor += exp(b - divider)
            end

            ene = ene / bet_nor
            mag = mag / bet_nor

            dev_ene = ene2 / bet_nor - ene^2
            dev_mag = mag2 / bet_nor - mag^2
            
            write(io, "$t\t$(real(ene))\t$(real(mag))\t$(real(dev_ene)/t^2)\t$(n*real(dev_mag)/t)\n")
        end
    end
end
using ITensors

sites = [Index(2), Index(2)]
links = [Index(2)]

A1 = ITensor([1, 1, 1, 1], sites[1], links[1])
A2 = ITensor([1, 1, 1, 1], sites[2], links[1])
@show A1
A = MPS([A1, A2])
normalize!(A)

ampo = OpSum()
ampo += "I", 1
x = MPO(ampo, sites)
@show x[1]

norm(A)
@show A[2]

B = 2.0 * A

@show B[2]
B = contract(x, A, cutoff = -Inf)

@show B[1]
norm(B)

n = 4 # system size
χ = 2 # bond dimension

# prepare random MPS
sites = siteinds("S=1/2", n)
links = collect(Iterators.map(i -> Index(χ, "l=$i"), 1:n-1))

A1 = ITensor(collect(Iterators.map(_ -> 1, 1:2*χ)), sites[1], links[1])
A2 = ITensor(collect(Iterators.map(_ -> 1, 1:4*χ)), links[1], sites[2], links[2])
A3 = ITensor(collect(Iterators.map(_ -> 1, 1:4*χ)), links[2], sites[3], links[3])
A4 = ITensor(collect(Iterators.map(_ -> 1, 1:2*χ)), links[3], sites[4])

ψ = MPS([A1, A2, A3, A4])
@show ψ[4]
norm(ψ)

ampo = OpSum()
ampo += "I",2
x = MPO(ampo, sites)
@show x[2]

ϕ = x * ψ
norm(ϕ)
@show ϕ[4]
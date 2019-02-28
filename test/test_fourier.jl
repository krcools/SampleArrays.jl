

using SampleArrays
using Test

include("utils/APS.jl")
include("utils/subsignal.jl")

Δt = 1.0
p = 10
χ = 10.0

ωs = π / Δt
ωm = ωs / χ

ϕ = APS(stepsize=Δt, tailsteps=p, oversampling=χ)

M = 6*p
f = subsignal(SampleArrays.axis(0,Δt,M), t->sin(ωm*t), ϕ)
B = SampleArrays.sampledsignal(f, -p*Δt:Δt:(M+p)*Δt)

C = SampleArrays.fouriertransform(B)
@show SampleArrays.offset(C)
@show SampleArrays.stepsize(C)

D = SampleArrays.inversefouriertransform(C, SampleArrays.offset(B))
@test length(B) == length(D)
@test SampleArrays.samples(B) ≈ SampleArrays.samples(D)

t = SampleArrays.axis(D)
@test eltype(D) == Float64

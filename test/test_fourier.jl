

using SampledSignals
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
f = subsignal(SampledSignals.axis(0,Δt,M), t->sin(ωm*t), ϕ)
B = SampledSignals.sampledsignal(f, -p*Δt:Δt:(M+p)*Δt)

C = SampledSignals.fouriertransform(B)
@show SampledSignals.offset(C)
@show SampledSignals.stepsize(C)

D = SampledSignals.inversefouriertransform(C, SampledSignals.offset(B))
@test length(B) == length(D)
@test SampledSignals.samples(B) ≈ SampledSignals.samples(D)

t = SampledSignals.axis(D)
@test eltype(D) == Float64

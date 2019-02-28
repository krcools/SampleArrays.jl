using Test

using SampleArrays

t1 = range(-1, step=0.1, length=20)
t2 = range(-0, step=0.1, stop=last(t1))

@test last(t1) ≈ last(t2)

s1 = sampledsignal(sin, t1)
s2 = sampledsignal(cos, t2)

function truncf(f,t0,t1)
    tol = eps((t0+t1)/2) * 1e3
    x -> (t0-tol <= x <= t1+tol) ? f(x) : zero(x)
end
truncf(f,ax) = truncf(f,first(ax),last(ax))

sin1 = truncf(sin, axis(s1))
cos2 = truncf(cos, axis(s2))

y1 = sin1.(t1)
@test samples(s1) ≈ y1

y2 = cos2.(t2)
@test samples(s2) ≈ y2

SampleArrays.add!(s1,s2)

s3 = sampledsignal(x->sin1(x)+cos2(x), t1)

@test s1 ≈ s3

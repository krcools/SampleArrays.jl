using Test

using SampleArrays

ax = range(1.0, step=0.5, length=4)
sp = [1.0, 2.0, -3.0, 4.5]

s1 = sampledsignal(ax, sp)
s2 = (x -> 2x) * s1

# [2.0, 3.0, 4.0, 5.0]
@test samples(s2) ≈ [2.0, 6.0, -12.0, 22.5]

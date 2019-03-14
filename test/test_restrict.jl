using SampleArrays
using Test

h = 0.123
ax1 = range(0.0, step=h, length=100)
sg1 = SampleArrays.sampledsignal(sin, ax1)

@test last(ax1) ≈ 99*h

ax2 = range(33*h, step=h, stop=200*h)
@test length(ax2) == 168

sg2 = SampleArrays.restrict(sg1, ax2)

ax3 = SampleArrays.axis(sg2)
@test first(ax3) ≈ 33h
@test last(ax3) ≈ 99h
@test length(ax3) == 67
@test length(ax3) == length(samples(sg2))
@test samples(sg2) ≈ sin.(ax3)

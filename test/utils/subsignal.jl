struct SubSignal{U,G,T}
    samples::Vector{U}
    offset::T
    stepsize::T
    numsteps::Int
    interpolator::G
end

function support(f::SubSignal)
    tail = support(f.interpolator)
    bulk = axis(f.offset, f.stepsize, f.numsteps)
    return range(first(tail)+first(bulk), stop=last(tail)+last(bulk), length=f.numsteps+2f.interpolator.tailsteps)
    # return (first(tail)+first(bulk)) : f.stepsize : (last(tail)+last(bulk))
end


"""
    subsignal(axis, f, interpolator)

Create a signal by interpolating f in the points of axis.
"""
function subsignal(axis, f, ϕ)

    @assert length(axis) >= 2

    offset = first(axis)
    stepsize = axis[2]-axis[1]
    numsteps = length(axis)
    samples = f.(axis)

    SubSignal(samples, offset, stepsize, numsteps, ϕ)
end

function subsignal(F::SampledSignals.SampledSignal, ϕ)

    ax = axis(F)
    offset = first(ax)
    stepsize = ax[2]-ax[1]
    numsteps = length(ax)
    samp = samples(F)

    SubSignal(samp, offset, stepsize, numsteps, ϕ)
end


function subsignal(axis::AbstractVector, samples::AbstractVector, ϕ)
    SubSignal(samples, offset(axis), stepsize(axis), length(axis), ϕ)
end


function (f::SubSignal)(x)

    N = f.numsteps
    F = f.samples
    ϕ = f.interpolator

    x0 = f.offset
    Δx = f.stepsize

    s0, s1 = support(f.interpolator)
    i0 = ceil(Int, (x+s0-x0)/Δx) + 1
    i1 = floor(Int, (x+s1-x0)/Δx) + 1

    U = eltype(F)
    T = typeof(ϕ(0))
    y = zero(promote_type(U,T))
    for i in i0:i1
        (1 <= i <= N) || continue
        xi = x0 + (i-1) * Δx
        y += F[i] * ϕ(x-xi)
    end

    return y
end

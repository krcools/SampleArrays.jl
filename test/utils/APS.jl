struct APS{T}
    stepsize::T
    tailsteps::Int
    oversampling::T
end


"""
    APS(;stepsize, tailsteps, oversampling)

Returns the Approximate Prolate Spheroidal function based on the supplied required keywords.
"""
function APS(;
        stepsize=nothing,
        bandwidth=nothing,
        tailsteps=error("required"),
        oversampling=error("required"))

    xor(stepsize == nothing, bandwidth == nothing) ||
        error("supply either 'stepsize' or 'bandwidth'")

    (bandwidth == nothing) ||
        return APS(π/bandwidth, tailsteps, oversampling)

    (stepsize == nothing) ||
        return APS(stepsize, tailsteps, oversampling)

    error("neither 'stepsize' nor 'bandwidth' were supplied.")
end

function (f::APS)(x)
    Δt, p, χ = f.stepsize, f.tailsteps, f.oversampling

    abs(x/(p*Δt)) >= 1 && return zero(x)

    ωs = π/Δt
    ω0 = ωs/χ
    ω1 = (ωs+ω0)/2
    ω2 = (ωs-ω0)/2

    A = ω1/ωs
    B = (x == 0) ? one(x) : sin(ω1*x)/(ω1*x)
    s = √(1-(x/(p*Δt))^2)
    C = sinh(ω2*p*Δt*s)/(sinh(ω2*p*Δt)*s)

    return A*B*C
end


"""
    support(f::APS)

Return the numerical support of an APS function. The function is guaranteed to be negilby small outside of the returned region.
"""
support(f::APS) = -f.stepsize*f.tailsteps, f.stepsize*f.tailsteps

"""
    stepsize(f::APS)

Returns the stepsize underlying the interpolator. This size is computed as the pi over the bandwidth.
"""
stepsize(f::APS) = f.stepsize


# """
#     precision(f::APS)
#
# Return the precision of an APS function. The returned precision is in the order of magnitude of the function value on the boundary of its numerical support.
# """
# function precision(f::APS)
#
#     Δt = f.stepsize
#     χ = f.oversampling
#     p = f.tailsteps
#
#     ωs = π / Δt
#     ωm = ωs / χ
#
#     1/sinh((ωs - ωm)/2*p*Δt)
# end

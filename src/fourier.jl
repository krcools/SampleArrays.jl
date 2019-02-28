"""
    fouriertransform(array, dt, t0, dim=1)

```math
    F(\\omega) = 1/\\sqrt{2 \\pi} \\int f(t) e^{-i \\omega t} dt
```
"""
function fouriertransform(a::AbstractArray, dt, t0, dim=1)
    n = size(a,dim)
    dω = 2π / (n*dt)
    b = fftshift(fft(a, dim), dim) * dt / sqrt(2π)
    ω0 = -dω * div(n,2)
    b .*= exp.(-im*t0*range(ω0,stop=ω0+(n-1)*dω,length=n))
    b, dω, ω0
end


"""
    inversefouriertransform(array,dω,ω0,t0=zero(ω0),dim=1)
"""
function inversefouriertransform(g,dω,ω0,t0=zero(ω0),dim=1)

    n = size(g,dim)
    dt = 2π/(n*dω)

    ω = axis(ω0, dω, n)
    a = n * dω * ifft(ifftshift(g .* exp.(im*ω*t0),dim),dim) / √(2π)

    return a, dt, t0
end


function fouriertransform(signal::SampledSignal)
    Y, dω, ω0 = fouriertransform(signal.samples, stepsize(signal), offset(signal))
    return SampledSignal(Y,dω,ω0)
end

function inversefouriertransform(signal::SampledSignal, t0=zero(typeof(offset(signal))))
    A, dt, t0 = inversefouriertransform(samples(signal), stepsize(signal), offset(signal), t0)
    return SampledSignal(real.(A),dt,t0)
end

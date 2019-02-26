module SampleArrays

import LinearAlgebra
using FFTW

export sampledsignal, axis, stepsize, offset, restrict, samples
export fouriertransform, inversefouriertransform
export SignalArray

include("axis.jl")
include("signal.jl")
include("signalarray.jl")
include("fourier.jl")

end # module

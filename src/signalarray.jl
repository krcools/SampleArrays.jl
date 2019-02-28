struct SignalArray{U,N1,T} <: AbstractArray{U,N1}
    array::Array{U,N1}
    stepsize::T
    offset::T
end

function SignalArray{U}(axis; dims=dims) where {U}
    n = length(axis)
    array = fill(zero(U), dims..., n)
    SignalArray(array, stepsize(axis), offset(axis))
end

Base.IndexStyle(::Type{<:SignalArray}) = Base.IndexCartesian()

Base.size(signal::SignalArray) = size(signal.array)
Base.length(signal::SignalArray) = length(signal.array)
Base.eltype(signal::SignalArray{U,N1,T}) where {U,N1,T} = SampledSignal{U,T}

function Base.getindex(signal::SignalArray, i)
    sampledsignal(axis(signal), view(signal.array,i,1:size(signal.array,2)))
end

function Base.setindex!(signal::SignalArray, v, i::Int)
    ax1 = axis(signal)
    ax2 = axis(v)
    @assert offset(ax1) ≈ offset(ax2)
    @assert stepsize(ax1) ≈ stepsize(ax2)
    @assert length(ax1) == length(ax2)
    signal.array[i,:] = samples(v)
end

function Base.:*(b::AbstractMatrix, a::SignalArray)
    SignalArray(b*a.array, a.stepsize, a.offset)
end

stepsize(signal::SignalArray) = signal.stepsize
offset(signal::SignalArray) = signal.offset
axis(signal::SignalArray) = axis(offset(signal), stepsize(signal), size(signal.array,2))

signal(a::SignalArray, i::Int) = SampledSignal(a.array[i,:], axis(a))

struct SampleBundle{U,N1,T} <: AbstractArray{U,N1}
    array::Array{U,N1}
    stepsize::T
    offset::T
end

function SampleBundle{U}(axis; dims=dims) where {U}
    n = length(axis)
    array = fill(zero(U), dims..., n)
    SampleBundle(array, stepsize(axis), offset(axis))
end

Base.IndexStyle(::Type{<:SampleBundle}) = Base.IndexCartesian()

Base.size(signal::SampleBundle) = size(signal.array)
Base.length(signal::SampleBundle) = length(signal.array)
Base.eltype(signal::SampleBundle{U,N1,T}) where {U,N1,T} = SampledSignal{U,T}

function Base.getindex(signal::SampleBundle, i)
    sampledsignal(axis(signal), view(signal.array,i,1:size(signal.array,2)))
end

function Base.setindex!(signal::SampleBundle, v, i::Int)
    ax1 = axis(signal)
    ax2 = axis(v)
    @assert offset(ax1) ≈ offset(ax2)
    @assert stepsize(ax1) ≈ stepsize(ax2)
    @assert length(ax1) == length(ax2)
    signal.array[i,:] = samples(v)
end

function Base.:*(b::AbstractMatrix, a::SampleBundle)
    SampleBundle(b*a.array, a.stepsize, a.offset)
end

stepsize(signal::SampleBundle) = signal.stepsize
offset(signal::SampleBundle) = signal.offset
axis(signal::SampleBundle) = axis(offset(signal), stepsize(signal), size(signal.array,2))

signal(a::SampleBundle, i::Int) = SampledSignal(a.array[i,:], axis(a))

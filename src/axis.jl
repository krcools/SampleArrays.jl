export snap

"""
    axis(x0,dx,n)
"""
axis(x0,dx,n) = range(x0, stop = x0 +(n-1)*dx, length = n)
offset(ax::AbstractVector) = first(ax)
stepsize(ax::AbstractVector) = (ax[2]-ax[1])


hull(a::AbstractRange, b::AbstractRange) = min(first(a),first(b)) : max(last(a),last(b))
set_addition(a::AbstractRange, b::AbstractRange) = (first(a)+first(b)) : step(a) : (last(a)+last(b))

"""

"""
function join_snapped(a::AbstractRange, b, step)
    # x0 = first(a) + first(b)
    # x1 = last(a) + last(b)
    # n0 = floor(Int, x0/step)
    # n1 = ceil(Int, x1/step)
    # return axis(n0*step, step, n1-n0+1)
    snap(set_addition(a,b), step)
end

function snap(r::AbstractRange, step)
    i0 = floor(Int,first(r)/step)
    i1 = ceil(Int, last(r)/step)
    return axis(i0*step, step, i1-i0+1)
end

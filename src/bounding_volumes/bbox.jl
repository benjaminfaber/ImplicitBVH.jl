"""
    $(TYPEDEF)

Axis-aligned bounding box, highly optimised for computing bounding volumes for triangles and
merging into larger bounding volumes.

Can also be constructed from two spheres to e.g. allow merging [`BSphere`](@ref) leaves into
[`BBox`](@ref) nodes.

# Methods
    # Convenience constructors
    BBox(lo::NTuple{3, T}, up::NTuple{3, T}) where T
    BBox{T}(lo::AbstractVector, up::AbstractVector) where T
    BBox(lo::AbstractVector, up::AbstractVector)

    # Construct from triangle vertices
    BBox{T}(p1, p2, p3) where T
    BBox(p1, p2, p3)
    BBox{T}(vertices::AbstractMatrix) where T
    BBox(vertices::AbstractMatrix)
    BBox{T}(triangle) where T
    BBox(triangle)

    # Merging bounding boxes
    BBox{T}(a::BBox, b::BBox) where T
    BBox(a::BBox{T}, b::BBox{T}) where T
    Base.:+(a::BBox, b::BBox)

    # Merging bounding spheres
    BBox{T}(a::BSphere{T}) where T
    BBox(a::BSphere{T}) where T
    BBox{T}(a::BSphere{T}, b::BSphere{T}) where T
    BBox(a::BSphere{T}, b::BSphere{T}) where T
"""
struct BBox{N, T}
    lo::NTuple{N, T}
    up::NTuple{N, T}
end

Base.eltype(::BBox{N, T}) where {N, T} = T
Base.eltype(::Type{BBox{N, T}}) where {N, T} = T
Base.ndims(::BBox{N, T}) where {N, T} = N
Base.ndims(::Type{BBox{N, T}}) where {N, T} = N



# Convenience constructors, with and without type parameter
function BBox{T}(lo::AbstractVector, up::AbstractVector) where T
    N = length(lo)
    BBox(NTuple{N, eltype(lo)}(lo), NTuple{3, eltype(up)}(up))
end
function BBox{N, T}(lo::AbstractVector, up::AbstractVector) where {N, T}
    BBox(NTuple{N, T}(T.(lo)), NTuple{N, T}(T.(up)))
end

function BBox(lo::AbstractVector, up::AbstractVector)
    BBox{length(lo), eltype(lo)}(lo, up)
end



# Constructors from triangles
function BBox{T}(p1, p2, p3) where T

    lower = (minimum3(p1[1], p2[1], p3[1]),
             minimum3(p1[2], p2[2], p3[2]),
             minimum3(p1[3], p2[3], p3[3]))

    upper = (maximum3(p1[1], p2[1], p3[1]),
             maximum3(p1[2], p2[2], p3[2]),
             maximum3(p1[3], p2[3], p3[3]))
   
    BBox{T}(lower, upper)
end

function BBox{T}(p1, p2) where T
    lower = (minimum(p1[1], p2[1]),
             minimum(p1[2], p2[2]))
    upper = (maximum(p1[1], p2[1]),
             maximum(p1[1], p2[1]))
    BBox{T}(lower, upper)
end

# Convenience constructors, with and without explicit type parameter
function BBox(p1, p2, p3)
    BBox{eltype(p1)}(p1, p2, p3)
end

function BBox(p1, p2)
    BBox{eltype(p1)}(p1, p2)
end

function BBox{T}(triangle_or_line) where T
    # Decompose triangle into its vertices.
    # Works transparently with GeometryBasics.Triangle, GeometryBasics.Line, Vector{SVector{N, T}}, etc.
    _BBox(triange_or_line, Val(length(triangle_or_line)))
end

function BBox(triangle_or_line)
    _BBox(triange_or_line, Val(length(triangle_or_line)))
end

function _BBox(triangle, ::Val{3})
    p1, p2, p3 = triangle
    BBox{eltype(p1)}(p1, p2, p3)
end

function _BBox(line, ::Val{2})
    p1, p2 = line
    BBox{eltype(p1)}(p1, p2)
end

function BBox{T}(vertices::AbstractMatrix) where T
    _BBox(vertices, Val(size(vertices, 2)))
end

function _BBox(vertices::AbstractMatrix, ::Val{3})
    BBox{eltype(vertices)}(@view(vertices[:, 1]), @view(vertices[:, 2]), @view(vertices[:, 3]))
end

function _BBox(vertices::AbstractMatrix, ::Val{2})
    BBox{eltype(vertices)}(@view(vertices[:, 1], @view(vertices[:, 2])))
end

function BBox(vertices::AbstractMatrix)
    _BBox{eltype(vertices)}(vertices, Val(size(vertices, 2)))
end


# Overloaded center function
center(b::BBox{3, T}) where T = (T(0.5) * (b.lo[1] + b.up[1]),
                                 T(0.5) * (b.lo[2] + b.up[2]),
                                 T(0.5) * (b.lo[3] + b.up[3]))

center(b::BBox{2, T}) where T = (T(0.5) * (b.lo[1] + b.up[1]),
                                 T(0.5) * (b.lo[2] + b.up[2]))

# Merge two bounding spheres
function BSphere(a::BSphere{3, T}, b::BSphere{3, T}) where T
    length = dist3(a.x, b.x)

    # a is enclosed within b
    if length + a.r <= b.r
        return BSphere{3, T}(b.x, b.r)

    # b is enclosed within a
    elseif length + b.r <= a.r
        return BSphere{3, T}(a.x, a.r)

    # Bounding spheres are not enclosed
    else
        frac = T(0.5) * ((b.r - a.r) / length + T(1))
        centre = (a.x[1] + frac * (b.x[1] - a.x[1]),
                  a.x[2] + frac * (b.x[2] - a.x[2]),
                  a.x[3] + frac * (b.x[3] - a.x[3]))
        radius = T(0.5) * (length + a.r + b.r)
        return BSphere{3, T}(centre, radius)
    end
end

# Merge two bounding spheres
function BSphere(a::BSphere{2, T}, b::BSphere{2, T}) where T
    length = dist2(a.x, b.x)

    # a is enclosed within b
    if length + a.r <= b.r
        return BSphere{2, T}(b.x, b.r)

    # b is enclosed within a
    elseif length + b.r <= a.r
        return BSphere{2, T}(a.x, a.r)

    # Bounding spheres are not enclosed
    else
        frac = T(0.5) * ((b.r - a.r) / length + T(1))
        centre = (a.x[1] + frac * (b.x[1] - a.x[1]),
                  a.x[2] + frac * (b.x[2] - a.x[2]))
        radius = T(0.5) * (length + a.r + b.r)
        return BSphere{2, T}(centre, radius)
    end
end

BSphere(a::BSphere{N, T}, b::BSphere{N, T}) where {N, T} = BSphere{T}(a, b)
Base.:+(a::BSphere, b::BSphere) = BSphere(a, b)


# Merge two bounding boxes
function BBox(a::BBox{3, T}, b::BBox{3, T}) where T
    lower = (minimum2(a.lo[1], b.lo[1]),
             minimum2(a.lo[2], b.lo[2]),
             minimum2(a.lo[3], b.lo[3]))

    upper = (maximum2(a.up[1], b.up[1]),
             maximum2(a.up[2], b.up[2]),
             maximum2(a.up[3], b.up[3]))

    BBox{T}(lower, upper)
end

function BBox(a::BBox{2, T}, b::BBox{2, T}) where T
    lower = (minimum2(a.lo[1], b.lo[1]),
             minimum2(a.lo[2], b.lo[2]))

    upper = (minimum2(a.up[1], b.up[1]),
             minimum2(a.up[2], b.up[2]))
    BBox{2, T}(lower, upper)
end

BBox(a::BBox{N, T}, b::BBox{N, T}) where {N, T} = BBox(a, b)
Base.:+(a::BBox, b::BBox) = BBox(a, b)


# Convert BSphere to BBox
function BBox{T}(a::BSphere{3, T}) where T
    lower = (a.x[1] - a.r, a.x[2] - a.r, a.x[3] - a.r)
    upper = (a.x[1] + a.r, a.x[2] + a.r, a.x[3] + a.r)
    BBox(lower, upper)
end

function BBox{T}(a::BSphere{2, T}) where {T}
    lower = (a.x[1] - a.r, a.x[2] - a.r)
    upper = (a.x[1] + a.r, a.x[2] + a.r)
    BBox(lower, upper)
end

function BBox(a::BSphere{N, T}) where {N, T}
    BBox{T}(a)
end

# Merge two BSphere into enclosing BBox
function BBox{T}(a::BSphere{3, T}, b::BSphere{3, T}) where T
    length = dist3(a.x, b.x)

    # a is enclosed within b
    if length + a.r <= b.r
        return BBox(b)

    # b is enclosed within a
    elseif length + b.r <= a.r
        return BBox(a)

    # Bounding spheres are not enclosed
    else
        lower = (minimum2(a.x[1] - a.r, b.x[1] - b.r),
                 minimum2(a.x[2] - a.r, b.x[2] - b.r),
                 minimum2(a.x[3] - a.r, b.x[3] - b.r))

        upper = (maximum2(a.x[1] + a.r, b.x[1] + b.r),
                 maximum2(a.x[2] + a.r, b.x[2] + b.r),
                 maximum2(a.x[3] + a.r, b.x[3] + b.r))

        return BBox(lower, upper)
    end
end

function BBox{T}(a::BSphere{2, T}, b::BSphere{2, T}) where T
    length = dist2(a.x, b.x)

    # a is enclosed within b
    if length + a.r <= b.r
        return BBox(b)

    # b is enclosed within a
    elseif length + b.r <= a.r
        return BBox(a)

    # Bounding spheres are not enclosed
    else
        lower = (minimum2(a.x[1] - a.r, b.x[1] - b.r),
                 minimum2(a.x[2] - a.r, b.x[2] - b.r))

        upper = (maximum2(a.x[1] + a.r, b.x[1] + b.r),
                 maximum2(a.x[2] + a.r, b.x[2] + b.r))

        return BBox(lower, upper)
    end
end

function BBox(a::BSphere{N, T}, b::BSphere{N, T}) where {N, T}
    BBox{T}(a, b)
end

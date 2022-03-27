"""
Represents an abstract N-dimensional grid of spins with undetermined topology
"""
abstract type AbstractSpinLattice{N} <: AbstractArray{Bool, N} end

"""
Represents an abstract N-dimensional square grid of spins
"""
abstract type AbstractSquareLattice{N} <: AbstractSpinLattice{N} end

"""
Represents an abstract 2-dimensional honeycomb grid of spins.
"""
abstract type AbstractHoneycombLattice <: AbstractSpinLattice{2} end

"""
The following section is to implement a basic SquareLattice concrete subtype.
This type is essentially just an encapsulated BitArray for which we will define
more specialized functions (and forward everything else).
"""
struct SquareLattice{N} <: AbstractSquareLattice{N}
    data::BitArray{N}
    SquareLattice{N}(dims, init = (x->rand(Bool))) where N = new(init.(BitArray{N}(undef, dims)))
end

SquareLattice(dims, init = (x -> rand(Bool))) = SquareLattice{length(dims)}(dims, init)

"""
The following section is to implement a basic HoneycombLattice concrete subtype.
This type largely behaves similarly to a 2D square lattice, except you can imagine
each row of cells is "offset" by half a cell, and has 6 neighbors (rather than 4)
"""
mutable struct HoneycombLattice <: AbstractHoneycombLattice
    data::BitArray{2}
    HoneycombLattice(dims, init = (x->rand(Bool))) = new(init.(BitArray{2}(undef, dims)))
end

Base.size(A::L) where L <: AbstractSpinLattice = size(A.data)
Base.IndexStyle(::Type{<:L}) where L <: AbstractSpinLattice = IndexLinear()

Base.getindex(A::L, i) where L <: AbstractSpinLattice = getindex(A.data, i)
Base.setindex!(A::L, v, i) where L <: AbstractSpinLattice = setindex!(A.data, v, i)
Base.getindex(A::L, ind::CartesianIndex) where L <: AbstractSpinLattice = getindex(A.data, CartesianIndex(mod1.(Tuple(ind), size(A.data))))
Base.setindex!(A::L, v, ind::CartesianIndex) where L <: AbstractSpinLattice = setindex!(A.data, v, CartesianIndex(mod1.(Tuple(ind), size(A.data))))

function Base.show(io::IO, mime::MIME"text/plain", x::T) where {N, T <: AbstractSquareLattice{N}}
    if N == 1
        write(io, String((x -> x ? '\u25A1' : '\u25A0').(x.data)))
    elseif N == 2
        for i in 1:size(x.data, 2)
            write(io, String((x -> x ? '\u25A1' : '\u25A0').(x.data[:, i])) * "\n")
        end
    else
        # TODO: handle the thing here with Base.Cartesian, probably
        show(io, mime, x.data)
    end
end

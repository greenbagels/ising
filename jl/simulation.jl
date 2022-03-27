include("lattice.jl")

abstract type AbstractIsingSimulation end

"""
Represents a single instance of an Ising simulation with measured values represented
using `T`

Note that this struct demonstrates a major tradeoff (though I wouldn't go as far
as to call it a "flaw") of Julia: we cannot guarantee type invariants once mutation
is allowed. This makes protecting against accidental or unintended modification
more difficult.
"""
mutable struct IsingSimulation{T, L} <: AbstractIsingSimulation
    "Contains the actual spin lattice representation"
    lattice::L
    "Site interaction energy"
    interaction::T
    "Total lattice energy given by the hamiltonian"
    energy::T
    "Total magnetization, aka sum of spins"
    magnetization::T
    "Magnetic field strength"
    field::T
    "Inverse temperature in units of kB"
    beta::T

    function IsingSimulation{T, L}(dims::NTuple{N, Int}=(10,10), init=(x -> rand(Bool)), U::T=1., H::T=0., beta::T=2.27) where {N, T <: Real, L <: AbstractSquareLattice{N}}
        lattice = L(dims, init)

        # Calculate the magnetization here.
        # The lattice stores boolean values such that (up, down) |-> (1, 0)
        # You might be tempted to first map the 0 values to -1 and then sum over the
        # lattice. But these transformations incur, at best, an extra conditional check
        # for each lattice element. Instead, we can use some math to perform a single
        # calculation at the end:

        # ntotal - nup    = ndown
        #    nup - ntotal = -ndown
        #   2nup - ntotal = nup - ndown
        #   2nup - ntotal = m
        ntot = length(lattice)
        M = 2 * count(lattice) - ntot

        # Now, for total energy.
        # The Ising Hamiltonian is E = - U * \sum_{ij} (s_i s_j) - H \sum_i s_i
        # where the sum over i, j is understood to be a sum over neighboring pairs.
        opp_spin_pairs = 0
        for site_index in CartesianIndices(lattice)
            n = ndims(lattice)
            for k in 1:n, dx in [-1, 1]
                offset = CartesianIndex(ntuple(i -> i == k ? dx : 0, n))
                # println("$(site_index) + $(offset) = $(site_index + offset)")
                # The idea: converting from {0,1} to {-1, 1} is messy; let's just
                # do it using boolean algebra instead, and reconcile the issue
                # outside the loop.
                opp_spin_pairs += xor(lattice[site_index], lattice[site_index + offset])
            end
        end

        # So now, opp_spin_pairs stores the (double-counted) number of pairs of
        # opposing spins. The number of pairs should just be 2*length(lattice),
        # so from this we can reconstruct the number of same-sign pairs :).
        # The same-sign pairs will each add a contribution of -U to the total energy,
        # while the opposite-sign pairs add a contribution of U. So our total
        # interaction energy is just -U(same_pairs - opp_pairs)
        #
        # But this is just U(2 * opp_pairs - total_pairs), and opp_pairs was
        # already double-counted!

        int_E = U * (opp_spin_pairs - ndims(lattice) * ntot)

        # Similarly, the field E is really just -H * M, which we already calculated.

        field_E = - H * M

        new(lattice, U, (int_E + field_E), M, H, beta)
    end

    function IsingSimulation{T, L}(dims::NTuple{N, Int}=(10,10), init=(x -> rand(Bool)), U::T=1., H::T=0., beta::T=2.27) where {N, T <: Real, L <: AbstractHoneycombLattice}
        lattice = L(dims, init)

        ntot = length(lattice)
        M = 2 * count(lattice) - ntot
        opp_spin_pairs = 0
        for site_index in CartesianIndices(lattice)
            # Imagine an 8x2 honeycomb looks like this
            #  _/x\_/x\_/x\_/x\
            # /x\x/x\x/x\x/x\x/
            # \x/ \x/ \x/ \x/ \
            # / \_/ \_/ \_/ \_/
            # \_/ \_/ \_/ \_/
            #
            # The XX marked cells are represented by a single row in a BitArray{2}
            # So clearly, we still only need a BitArray{2} of shape (8,2) to store
            # this data. However, in terms of locality, notice that the "stencil"
            # for a cell's neighbors depends on the parity of its (one-indexed) x
            # coordinate; if x is odd/even, it has 3 neighbors below/above and
            # 1 neighbor above/below. It always has 2 sideways neighbors, though.
            offsets = CartesianIndex.(site_index[1] % 2 == 1 ?
                                      [ (0, -1), (-1, 0), (1, 0), (-1,  1), (0,  1), (1,  1)] :
                                      [ (0,  1), (-1, 0), (1, 0), (-1, -1), (0, -1), (1, -1)]   )
            for offset in offsets
                opp_spin_pairs += xor(lattice[site_index], lattice[site_index + offset])
            end
        end
        # This situation is identical to the square lattice up until we count
        # total_pairs. In our case, there are (6*ntot)/2 = 3 ntot pairs,
        # not 2 ntot pairs! (this corresponds to a cell interacting with 6
        # neighbors, rather than 4, a 50% increase).
        int_E = U * (opp_spin_pairs - 3 * ntot)
        field_E = - H * M

        new(lattice, U, (int_E + field_E), M, H, beta)
    end
end

function try_flip!(sim::S, site_index::CartesianIndex = CartesianIndex(ntuple(i -> rand(1:size(sim.lattice, i)), ndims(sim.lattice)))) where {T <: Real, L <: SquareLattice, S <: IsingSimulation{T, L}}
    # The idea: for a fixed i, the interaction energy with its neighbors
    # (with spin s_j) is - U * \sum_j s_i s_j = - U * s_i \sum_j s_j.
    # If we flip s_i, then it becomes U * s_i \sum_j s_j, and the
    # difference between these is just 2 U * s_i \sum_j s_j.

    # Similarly, the field energy is -H*s_i, which becomes H*s_i, a change
    # of 2H*s_i. So the total change is 2U * s_i \sum_j s_j + 2H * s_i
    #                                 = 2 s_i (U \sum_j s_j + H)

    # But \sum_j s_j is just (up_spins - down_spins)
    #                      = (up_spins - (2*ndim - up_spins))
    #                      = (2*up_spins - 2*ndim))
    #                      = 2(up_spins - ndim)
    dims = ndims(sim.lattice)
    ntot = length(sim.lattice)
    up_spins = 0
    for k in 1:dims, dx in [-1, 1]
        offset = CartesianIndex(ntuple(i -> i == k ? dx : 0, dims))
        up_spins += sim.lattice[site_index + offset]
    end
    s_i = 2 * sim.lattice[site_index] - 1
    dE = 2. * s_i * (2 * sim.interaction * (up_spins - dims) + sim.field)

    if dE <= 0. || (rand() < exp(-sim.beta * dE))
        sim.lattice[site_index] = !sim.lattice[site_index]
        sim.energy += dE
        sim.magnetization += s_i
    end
end

function try_flip!(sim::S, site_index::CartesianIndex = CartesianIndex(ntuple(i -> rand(1:size(sim.lattice, i)), ndims(sim.lattice)))) where {T <: Real, L <: HoneycombLattice, S <: IsingSimulation{T, L}}
    # The idea: for a fixed i, the interaction energy with its neighbors
    # (with spin s_j) is - U * \sum_j s_i s_j = - U * s_i \sum_j s_j.
    # If we flip s_i, then it becomes U * s_i \sum_j s_j, and the
    # difference between these is just 2 U * s_i \sum_j s_j.

    # Similarly, the field energy is -H*s_i, which becomes H*s_i, a change
    # of 2H*s_i. So the total change is 2U * s_i \sum_j s_j + 2H * s_i
    #                                 = 2 s_i (U \sum_j s_j + H)

    # But \sum_j s_j is just (up_spins - down_spins)
    #                      = (up_spins - (6 - up_spins))
    #                      = (2*up_spins - 6))
    ntot = length(sim.lattice)
    up_spins = 0
    offsets = CartesianIndex.(site_index[1] % 2 == 1 ?
                              [ (0, -1), (-1, 0), (1, 0), (-1,  1), (0,  1), (1,  1)] :
                              [ (0,  1), (-1, 0), (1, 0), (-1, -1), (0, -1), (1, -1)]   )
    for offset in offsets
        up_spins += sim.lattice[site_index + offset]
    end
    s_i = 2 * sim.lattice[site_index] - 1
    dE = 2. * s_i * (2 * sim.interaction * (up_spins - 3) + sim.field)

    if dE <= 0. || (rand() < exp(-sim.beta * dE))
        sim.lattice[site_index] = !sim.lattice[site_index]
        sim.energy += dE
        sim.magnetization += s_i
    end
end

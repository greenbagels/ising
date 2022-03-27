include("simulation.jl")

function test_lattice_jl()
    println("Testing lattice.jl:\n")
    println("\n    - Making empty 5x5 square lattice...")
    A = SquareLattice((5, 5), (x) -> false)
    display(A)
    println("\n    - Making random 5x5 square lattice...")
    A = SquareLattice((5, 5))
    display(A)
    println()

    println("\n    - Making empty 5x5 honeycomb lattice...")
    A = HoneycombLattice((5, 5), (x) -> false)
    display(A)
    println("\n    - Making random 5x5 honeycomb lattice...")
    A = HoneycombLattice((5, 5))
    display(A)
    println()
end

function test_simulation_jl()
    println("\n\n\n\nTesting simulation.jl:")
    println("    - Testing random square lattice:")
    A = IsingSimulation{Float64, SquareLattice{2}}((10, 10)#=, (x) -> false=#)
    n = length(A.lattice)
    println("    - Average magnetization (normalized): $(A.magnetization/n)")
    println("    - Average energy (normalized): $(A.energy/n)")
    println("    - `try_flip`ing 100 times...")
    for i in 1:100
        try_flip!(A)
    end
    println("    - Average magnetization (normalized): $(A.magnetization/n)")
    println("    - Average energy (normalized): $(A.energy/n)\n\n\n")

    println("    - Testing random honeycomb lattice:")
    A = IsingSimulation{Float64, HoneycombLattice}((10, 10)#=, (x) -> false=#)
    n = length(A.lattice)
    println("    - Average magnetization (normalized): $(A.magnetization/n)")
    println("    - Average energy (normalized): $(A.energy/n)")
    println("    - `try_flip`ing 100 times...")
    for i in 1:100
        try_flip!(A)
    end
    println("    - Average magnetization (normalized): $(A.magnetization/n)")
    println("    - Average energy (normalized): $(A.energy/n)")
end

test_lattice_jl()
test_simulation_jl()

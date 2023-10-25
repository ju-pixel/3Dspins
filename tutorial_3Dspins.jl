# using Pkg
# Pkg.add("GLMakie")
# Pkg.add("GeometryBasics")
# Pkg.add("Colors")

using GLMakie, GeometryBasics, Colors

GLMakie.activate!()

struct Par
    N::Int64      # number of cells across the lattice
    L::Float64    # length of the elements
    J::Float64   # NN distance (=1 here)
    β::Float64    # angle spins make with the x axis 
end

#____________________________________________________________________


function get_spin_centers(C::Par)

    # generate a periodic lattice of points

    # side length l in terms of J spacing
    # l is the distance from outermost spin centers on each side
    
    l = (2C.N - 1) * C.J

    spin_centers = Vector{Float64}[]
    
    for x = -l/2 : C.J : l/2 
        for y = l/2 : -C.J : -l/2
            push!(spin_centers, Float32[x,y,0.0])
        end
    end

    return spin_centers
end

#____________________________________________________________________


function create_T2_config(C::Par)

    # individual spin vectors grouped in "cells" of 4
    spin_1 = [cos(C.β), -sin(C.β), 0]
    spin_2 = [cos(C.β), cos(C.β), 0]
    spin_3, spin_4 = spin_2, spin_1

    # number of columns in the lattice (it is also a # of spins per side)
    columns = 2C.N

    μ_hat = Vector{Float64}[]

    for col = 1:columns
        spins = isodd(col) ? [spin_1, spin_2] : [spin_3, spin_4]
        append!(μ_hat, repeat(spins, C.N))
    end

    return μ_hat
end


#____________________________________________________________________


function create_random_config(C::Par)

    # individual spin vectors grouped in "cells" of 4
    spin_1 = [cos(C.β), -sin(C.β), 0]
    spin_2 = [cos(C.β), cos(C.β), 0]
    spin_3, spin_4 = spin_2, spin_1

    # number of columns in the lattice (it is also a # of spins per side)
    columns = 2C.N
    
    μ_hat_rand = Vector{Float64}[]

    for col = 1:columns
        spins = isodd(col) ? [spin_1 *rand([1,-1]) , spin_2 *rand([1,-1])] : [spin_3 *rand([1,-1]), spin_4 *rand([1,-1])]
        append!(μ_hat_rand, repeat(spins, C.N))
    end

    return μ_hat_rand
end


#____________________________________________________________________


function plot_spins(N, L, J, β)

    C = Par(N, L, J, β)

    spin_centers = get_spin_centers(C)
    #μ_hat = create_T2_config(C)
    μ_hat_rand = create_random_config(C)

    spins_per_side = 2C.N
    num_spins = spins_per_side^2

    spin_center = [Point3f(i) for i in spin_centers] 
    #spin_conf = [(C.L/2) * Vec3f(i) for i in μ_hat]
    spin_conf = [(C.L/2) * Vec3f(i) for i in μ_hat_rand]

    fig = Figure(resolution=(720, 720), dpi=600)
    ax = Axis3(fig[1,1], aspect=:data, perspectiveness=0.6)

    arrows!(ax, spin_center, spin_conf,          # axes, positions, directions
            linecolor=:darkgrey, 
            arrowcolor=RGB(0.43, 0.43, 0.43), 
            quality=32,                          # Sets the quaity of the arrow.
            arrowsize=Vec3f0(0.3, 0.3, C.L/2),   # Scales the size of the arrow head. 
            linewidth=0.15, 
            align=:center,)                      # Sets how arrows are positioned.

    hidedecorations!(ax)
    hidespines!(ax)

    #save("$(num_spins)_orderedConfig_spins3D.png", fig)
    save("$(num_spins)_randConfig_spins3D.png", fig)

    return fig
end

#plot_spins(N,  L,  J,  β)
plot_spins(10, 1.0, 1.0, π/4)



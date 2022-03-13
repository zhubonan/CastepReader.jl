#=
Code for writing for CHGCAR file containing 3D grid data
=#


module Chgcar

using Printf
using LinearAlgebra

"""Enumerate fortran float where a leading zero is written"""
function ffloat(f)
    string = @sprintf "%.10E" f
    if f > 0
        return @sprintf "0.%s%sE%+03d" string[1:1] string[3:12] Int(ceil(log10(f)))
    end
    return @sprintf "-0.%s%sE%+03d" string[2:2] string[4:13]  Int(ceil(log10(-f)))
end

"""
Write out a CHGCAR file with 
"""
function write_chgcar(io::IO, datasets; lattice_col=Diagonal([10, 10, 10]), 
                      species=(:Si, ), coords=[0., 0., 0.][:, :], comment="VolumetricData")
    write(io, "$(comment)\n")
    write(io, "   1.00000000000000\n") 
    write(io, @sprintf(" %12.6f%12.6f%12.6f\n", lattice_col[:, 1]...))
    write(io, @sprintf(" %12.6f%12.6f%12.6f\n", lattice_col[:, 2]...))
    write(io, @sprintf(" %12.6f%12.6f%12.6f\n", lattice_col[:, 3]...))
    # Symbols and counts
    usymbols = unique(species)
    counts = zeros(Int, length(usymbols))
    for (i, symbol) in enumerate(unique(species))
        write(io, @sprintf "%5s" symbol)
        counts[i] = count(x -> x == symbol, species)
    end
    write(io, "\n")
    for c in counts
        write(io, @sprintf "%6d" c)
    end
    write(io, "\n")
    # Write the coordinates
    write(io, "Direct\n")
    for coord in eachcol(coords)
        @sprintf "%10.6f%10.6f%10.6f\n" coord[1] coord[2] coord[3]
    end
    write(io, " \n") 
    # Now starts the volumetric data set
    for data in datasets
        write_data(io, data)
    end
end

"""
Write a CHGCAR file
"""

function write_chgcar(fname::AbstractString, datasets; lattice_col=Diagonal([10, 10, 10]), 
                      species=(:Si, ), coords=[0., 0., 0.][:, :], comment="VolumetricData")
    open(fname, "w") do io
        write_chgcar(io, datasets; lattice_col, species, coords, comment)
    end
end

"""
Write the volumetric data block of the CHGCAR file
"""
function write_data(io::IO, data)

    ngx, ngy, ngz = size(data)
    write(io, "   $ngx   $ngy   $ngz\n")
    count = 1
    for k=1:ngz, j=1:ngy, i=1:ngx
        # Start of a new line
        count % 5 == 1 && write(io, " ")
        # Write data
        if count % 5 == 0
            write(io, ffloat(data[i, j, k]), "\n")
        else
            write(io, ffloat(data[i, j, k]), " ")
        end
        count += 1
    end
    # Write a final line break
    if count % 5 != 1
        write(io, "\n")
    end

end

export write_chgcar
end # module chgcar

import .Chgcar: write_chgcar

"""
"""
function write_chgcar(fname, density::ChargeDensity;sum_and_spin=true)
    if sum_and_spin
        total = (density.density[:, :, :, 1] .+ density[:, :, :, 2]) ./ 2
        diff = (density.density[:, :, :, 1] .- density[:, :, :, 2]) ./ 2
        datasets = (total, diff)
    else
        datasets = (density[:, :, :, 1], density[:, :, :, 2])
    end
    write_chgcar(fname, datasets;lattice_col=density.basis, species=(:Si, ), coords=[0., 0., 0.][:, :], )
end
#=
Wavefunction related routines

The goal is to recover the full wavefunction grid from the planewave coefficients
=#

"""
    coeff_to_recip(coeff_array, nwaves_at_kp, grid_coords, ngx, ngy, ngz)

Convert planewave coefficients in the reduced representation to grids
"""
function coeff_to_recip(coeff_array, nwaves_at_kp, grid_coords, ngx, ngy, ngz)
    npw, nspinor, band_max, nkpts, nspins = size(coeff_array)
    grid_size = Int[ngx, ngy, ngz]
    indices = coord_to_indices(grid_coords, grid_size)
    # Grid data index with ngx, ngy, ngz, nspinor, band_max, nkpts, nspins
    # Note that the kpoints order may be different from those read from other sections
    # so be careful with that.....
    grid = zeros(ComplexF64, ngx, ngy, ngz, nspinor, band_max, nkpts, nspins)

    for is in 1:nspins
        for ik in 1:nkpts
            for ib in 1:band_max, ispinor in 1:nspinor
                for ipw in 1:nwaves_at_kp[ik] 
                    grid[indices[1, ipw, ik], indices[2, ipw, ik],indices[3, ipw, ik], 
                         ispinor, ib, ik, is] = coeff_array[ipw, ispinor, ib, ik, is]
                end
            end
        end
    end
    grid
end

"""
Convert G-vector coordinates (unit of reciprocal lattice basis) to 1-based grid indices
"""
function coord_to_indices(grid_coords, grid_size)
    indices = copy(grid_coords)
    _, npw, nkpts = size(grid_coords)
    # The grid_coords have sizes in -n/2 to n/2 - we move it to 1 to n
    for k in 1:nkpts
        for j in 1:npw
            for i in 1:3
                if indices[i, j, k] < 0
                    indices[i, j, k] += grid_size[i]
                end
            end
        end
    end
    # Switch to 1 based indexing
    indices .+= 1
end

"""
Load wavefunction using data read from the checkpoint file.

Returns the reciprocal representation on the grids. 
Needs to be FFTed to get to the real space representation.
"""
function wavefunction(data_dict::Dict)
    wave = data_dict[:WAVEFUNCTION_DATA]
    coeff_to_recip(wave.coeffs, wave.nwaves_at_kp, wave.pw_grid_coord, wave.ngx, wave.ngy, wave.ngz)
end
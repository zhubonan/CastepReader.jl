#=
Wavefunction related routines

The goal is to recover the full wavefunction grid from the planewave coefficients
=#
import FFTW
using FFTW
using LinearAlgebra

export WaveFunction, ChargeDensity, chargedensity, sliceband, slicekpoint, slicespin, slicespinor


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
    coeff_to_recip(coeff_array, data_dict::Dict)

Convert plane wave coefficients to a reciprocal space grid.
"""
function coeff_to_recip(coeff_array, data_dict::Dict)
    wave = data_dict[:WAVEFUNCTION_DATA]
    coeff_to_recip(coeff_array, wave.nwaves_at_kp, wave.pw_grid_coord, wave.ngx, wave.ngy, wave.ngz)
end

"""
Convert full grid representation back to reciprocal coefficients
The original specifications of the wave function needs to be reused
"""
function recip_to_coeff(grid, coeff_array_orig, nwaves_at_kp, grid_coords, ngx, ngy, ngz)
    npw, nspinor, band_max, nkpts, nspins = size(coeff_array_orig)
    grid_size = Int[ngx, ngy, ngz]
    indices = coord_to_indices(grid_coords, grid_size)
    # Grid data index with ngx, ngy, ngz, nspinor, band_max, nkpts, nspins
    # Note that the kpoints order may be different from those read from other sections
    # so be careful with that.....
    # Allocate the coefficients
    coeff_array = similar(coeff_array_orig)
    fill!(coeff_array, 0.)

    for is in 1:nspins
        for ik in 1:nkpts
            for ib in 1:band_max, ispinor in 1:nspinor
                for ipw in 1:nwaves_at_kp[ik] 
                    coeff_array[ipw, ispinor, ib, ik, is] = grid[indices[1, ipw, ik], indices[2, ipw, ik],indices[3, ipw, ik], 
                         ispinor, ib, ik, is] 
                end
            end
        end
    end
    coeff_array
end

function recip_to_coeff(grid, data_dict)
    wave = data_dict[:WAVEFUNCTION_DATA]
    recip_to_coeff(grid, wave.coeffs, wave.nwaves_at_kp, wave.pw_grid_coord, wave.ngx, wave.ngy, wave.ngz)
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


"""
Type for storing a wavefunction
"""
mutable struct WaveFunction
    wave::Array{ComplexF64, 7}
    real_space::Bool
    ngx::Int
    ngy::Int
    ngz::Int
    nspinor::Int
    nbands::Int
    nkpts::Int
    nspins::Int
    basis::Matrix{Float64}
    "Occupations in (ib, ik, is)"
    occupations::Array{Float64, 3}
    "Eigenvalues in (ib, ik, is)"
    eigenvalues::Array{Float64, 3}
    "Weight of the kpoints"
    kpoint_weights::Array{Float64, 1}
    "Coordinates of the kpoints"
    kpoints::Array{Float64, 2}
end

struct ChargeDensity{T}
    density::Array{T, 4}
    ngx::Int
    ngy::Int
    ngz::Int
    nspins::Int
    basis::Matrix{Float64}
end

"""
    ChargeDensity(density)

Construction from a density array.
"""
function ChargeDensity(density, basis=diagm([10., 10., 10.]))
    ChargeDensity(density, size(density)..., basis)
end

"""
    ChargeDensity(ngx, ngy, ngz)

Construction an empty ChargeDensity.
"""
function ChargeDensity{T}(ngx::Int, ngy::Int, ngz::Int) where {T}
    ChargeDensity(zeros(T, ngx, ngy, ngz))
end

"""
    WaveFunction(data_dict::Dict)

Construct the WaveFunction from a Dict of read data from castep_bin/check file
"""
function WaveFunction(data_dict::Dict)
    wave = wavefunction(data_dict)
    ngx, ngy, ngz, nspinor, nbands, nkpts, nspins = size(wave)
    basis = data_dict[:REAL_LATTICE] .* Bohr2Ang
    eigenvalues = data_dict[:EIGENVALUES] .* Ha2eV
    occupations = data_dict[:OCCUPATIONS]
    
    kpoints_weights = data_dict[:KPOINT_WEIGHTS]
    kpoints = data_dict[:KPOINTS]
    WaveFunction(wave, false, ngx, ngy, ngz, nspinor, nbands, nkpts, nspins, 
                 basis, occupations, eigenvalues, kpoints_weights, kpoints)
end

slicespinor(wave::WaveFunction, ispinor) = selectdim(wave.wave, 4, ispinor)

sliceband(wave::WaveFunction, ib) = selectdim(wave.wave, 5, ib)

slicekpoint(wave::WaveFunction, ik) = selectdim(wave.wave, 6, ik)

slicespin(wave::WaveFunction, is) = selectdim(wave.wave, 7, is)

"""
    FFTW.fft!(wavef::WaveFunction)
In place FFT to transfer the wavefunction to the frequency space

The in-place transformation is apply slice by slice
"""
function FFTW.fft!(wavef::WaveFunction)
    @assert wavef.real_space "WaveFunction is already in the frequency space"
    xyz_slice = zeros(ComplexF64, wavef.ngx, wavef.ngy, wavef.ngz) 
    for i=1:wavef.nspins,j=1:wavef.nkpts,k=1:wavef.nbands,l=1:wavef.nspinor
        xyz_slice .= wavef.wave[:, :, :, l, k, j, i]
        wavef.wave[:, :, :, l, k, j, i] = FFTW.fft(xyz_slice)
    end
    wavef.real_space = false
    wavef
end

"""
    fft(wavef::WaveFunction)
FFT to transfer the wavefunction to the frequency space
"""
function FFTW.fft(wavef::WaveFunction)
    @assert wavef.real_space "WaveFunction is already in the frequency space"
    new_wave = deepcopy(wavef)
    FFTW.fft!(new_wave)
end


"""
In place FFT to transfer the wavefunction to the real space

The in-place transformation is apply slice by slice
"""
function FFTW.ifft!(wavef::WaveFunction)
    @assert ~wavef.real_space "WaveFunction is already in the real space"
    xyz_slice = zeros(ComplexF64, wavef.ngx, wavef.ngy, wavef.ngz) 
    for i=1:wavef.nspins,j=1:wavef.nkpts,k=1:wavef.nbands,l=1:wavef.nspinor
        xyz_slice .= wavef.wave[:, :, :, l, k, j, i]
        wavef.wave[:, :, :, l, k, j, i] = FFTW.ifft(xyz_slice)
    end
    wavef.real_space = true
    wavef
end

"""
    fft(wavef::WaveFunction)
FFT to transfer the wavefunction to the frequency space
"""
function FFTW.ifft(wavef::WaveFunction)
    @assert ~wavef.real_space "WaveFunction is already in the real space"
    new_wave = deepcopy(wavef)
    FFTW.ifft!(new_wave)
end

ensure_recip!(wavef::WaveFunction) = wavef.real_space && fft!(wavef)
ensure_real!(wavef::WaveFunction) = wavef.real_space || ifft!(wavef)

"""
    chargedensity(wavef::WaveFunction, kidx, bidx, sidx)

Construct charge density for given kpoints, band, and spin indices
"""
function chargedensity(wavef::WaveFunction, bidx, kidx, sidx; weighted=true, ivbm=-1)
    orig_state = wavef.real_space
    ensure_real!(wavef)
    density = zeros(wavef.ngx, wavef.ngy, wavef.ngz, wavef.nspins)
    # Force occupations
    if ivbm > 1
        occ = zeros(size(wavef.occupations))
        for ib =1:ivbm
            occ[ib, : , :] .= 1.0
        end
    else
        occ = wavef.occupations
    end
    for is in sidx,ik in kidx,ib in bidx,ispinor=1:wavef.nspinor
        for z=1:wavef.ngz, y=1:wavef.ngy, x=1:wavef.ngx
            val = norm(wavef.wave[x, y, z, ispinor, ib, ik, is])
            # Do we weight by kpoint weights and occupations? 
            if weighted
                val = val * occ[ib, ik, is] * wavef.kpoint_weights[ik]
            end
            density[x, y, z, is] += val 
        end
    end
    orig_state == false && ensure_recip!(wavef)
    ChargeDensity(density, wavef.basis)
end

"""
    chargedensity(wavef::WaveFunction)

Compute charge density from the wavefunction

Note that the density does not include the augmentation charges involved in ultrasoft pseudopotentials
"""
function chargedensity(wavef::WaveFunction; ivbm=-1)
    @assert ivbm > 1 || sum(wavef.occupations[:, 1, 1]) > 0. "VBM index must be passed for unweighted WaveFunction"
    chargedensity(wavef, 1:wavef.nbands, 1:wavef.nkpts, 1:wavef.nspins;ivbm=ivbm)
end

"""
    chargedensity(wavef::WaveFunction, ib, ik, is)

Get the charge density for a single slice of band, kpoint, spin channel
"""
function chargedensity(wavef::WaveFunction, ib::Int, ik::Int, is::Int)
    chargedensity(wavef, ib:ib, ik:ik, is:is;weighted=false)
end


"""
Interpolate charge density from `c1` into `c2` with a finner grid
"""
function fine_grid_interpolate(c1::ChargeDensity{T}, c2::ChargeDensity{T}) where {T}
    # Zero the grid to be interpolated into
    fill!(0., c2.density)
    ngx, ngy, ngz = c1.ngx, c1.ngy, c1.ngz
    ngxf, ngyf, ngzf = c2.ngx, c2.ngy, c2.ngz
    @assert (ngxf > ngx) && (ngyf > ngy) && (ngzf > ngz)

    workspace = zeros(ComplexF64, ngx, ngy, ngz)
    workspace_fine = zeros(ComplexF64, ngxf, ngyf, ngzf)
    for is in 1:c1.nspins
        workspace .= c1.density[:, :, :, is]
        # To reciprocal space
        fft!(workspace)
        for z = 1:ngz, y=1:ngy, x=1:ngx
            # Compute the index in the interpolated grid in the reciprocal space
            x < ngx/2+1 ? fx = x : fx = ngxf - ngx + x
            y < ngy/2+1 ? fy = y : fy = ngyf - ngy + y
            z < ngz/2+1 ? fz = z : fz = ngzf - ngz + z
            work_space_fine[fx, fy, fz] = workspace[x, y, z]
        end
        # Back to real space
        ifft!(workspace_fine)
        # Function barrier for copying back the density
        copy_density(c2, workspace_fine, is)
    end
end

"""
Copy density from a temporary array in to the spin channel
"""
function copy_density(density::ChargeDensity{ComplexF64}, tmp, spin) 
    density.density[:, :, :, spin] .= tmp
end

"""
Copy density from a temporary array in to the spin channel
"""
function copy_density(density::ChargeDensity{Float64}, tmp, spin) 
    density.density[:, :, :, spin] .= real.(tmp)
end


"""
Computed weighted density from some given weights by bands
"""
function weighted_density_by_bands(fname, weights)
    data = read_castep_check(fname)
    nbands = data[:NBANDS]
    nspins = data[:NSPINS]
    nkpts = data[:NKPTS]
    @assert size(weights, 1) nbands
    wavef = WaveFunction(data)
    ifft!(wavef)
    chargedensity(wavef, ensure_three_dims(weights, nkpts, nspins), data[:KPOINTS_WEIGHTS])
end


ensure_three_dims(weights::AbstractArray{T, 1}, d1, d2) where {T} = repeat(weights[:, :, :], d1, d2)
ensure_three_dims(weights::AbstractArray{T, 2}, d1, d2) where {T} = repeat(weights[:, :, :], d2)
ensure_three_dims(weights::AbstractArray{T, 3}, d1, d2) where {T} = weights

"""
Update wavefunction of an existing file
"""
function CastepBin.update_wavefunction(src, dst, wavef::WaveFunction, data::Dict=read_castep_check(src))
    @assert wavef.real_space == false "Wave must be in the reciprocal space"
    new_coeff = recip_to_coeff(wavef.wave, data)
    CastepBin.update_wavefunction(src, dst, new_coeff)
end

"""
Update wavefunction of an existing file
"""
function CastepBin.update_wavefunction!(src, wavef::WaveFunction, data::Dict=read_castep_check(src))
    @assert wavef.real_space == false "Wave must be in the reciprocal space"
    new_coeff = recip_to_coeff(wavef.wave, data)
    CastepBin.update_wavefunction!(src, new_coeff)
end
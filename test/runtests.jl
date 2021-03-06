using Test
using CastepReader
using Unitful
using UnitfulAtomic
using FFTW

const cr=CastepReader
const testdata = joinpath(pkgdir(cr), "test/testdata")

gaas = joinpath(testdata, "GaAs")

# Test file name functions
@test cr.suffix("a.b") == "b"
@test cr.stem("x/a.b") == "x/a"
@test cr.filename("c/x/a.b") == "a.b"
@test cr.stem(cr.filename("c/x/a.b")) == "a"

@testset "Bands" begin
    kpts, kw, kidx, eigen, ne, efermi, cell_file = cr.read_bands_castep_raw(joinpath(gaas, "GaAs.bands"))
    @test size(kpts) == (3, 10)
    @test size(kw)[1] == 10
    @test size(eigen) == (1, 10, 40)
    @test kpts[:, 8] == [0., 0., 0.]
    @test kw[8] ≈ 0.008
    @test kidx[8] == 8
    @test ne == 28
    @test efermi[1] == uconvert(u"eV", 0.113556u"Eh_au")
    @test cell_file[1] == uconvert(u"Å", 5.432963u"a0_au")

    # Read the optical matrix
    ns, nk, nb = size(eigen)
    bin = cr.read_ome_bin(joinpath(gaas, "GaAs.ome_bin"), nb, nk, ns)
    @test size(bin) == (nb, nb, 3, nk, ns)
    @test bin[1] ≈ (0.00022207179206183356 - 0.0im) * cr.HaBohr

    # TODO - add test for cst_ome
    cst = cr.read_cst_ome(joinpath(gaas, "GaAs.cst_ome"), nb, nk, ns)

    @test size(cst) == (nb, nb, 3, nk, ns)
    @test cst[1] ≈ (0.0004183035272252899 - 0.0im) * cr.HaBohr2perAng

end



# Test reading the cell file
@testset "read cell" begin
    cell_file = cr.read_cell_castep(joinpath(gaas, "GaAs.cell"))
    @test :symmetry_ops in keys(cell_file.blocks)
    @test begin
        :kpoints_mp_spacing in keys(cell_file.kws) && cell_file.kws[:kpoints_mp_spacing] == "0.07"
    end
    cellmat = cr.read_lattice(cell_file)
    @test cellmat[1] == 2.875u"Å"

    symm, disp = cr.read_symm_disp(cell_file)
    @test begin
        symm[1, 1, 2] == 0. && symm[2, 1, 2] == -1. && all(disp .== 0.)
    end
end

@testset "castep_bin" begin
    bin = joinpath(gaas, "GaAs.castep_bin")    
    output = read_castep_check(bin)
    @test :FORCES in keys(output)
    @test output[:VERSION_NUMBER] ≈ 20.1100006103
    @test output[:NUM_IONS]  == 2
    @test output[:NUM_SPECIES]  == 2
    @test output[:MAX_IONS_IN_SPECIES]  == 1

    @test size(output[:FORCES]) == (3,1, 2)
end

@testset "wave" begin
    checkfile = joinpath(gaas, "GaAs.check")    
    output = read_castep_check(checkfile)
    wavef = CastepReader.WaveFunction(output)
    ngx = 25
    @test wavef.ngx == ngx
    out = CastepReader.sliceband(wavef, 1) 
    @test size(out, 1) == ngx
    ifft!(wavef)
    @test wavef.real_space == true
    fft!(wavef)
    @test wavef.real_space == false
    # Charge density - with out interpolation
    ifft!(wavef)
    out = CastepReader.chargedensity(wavef, ones(wavef.nbands, wavef.nkpts, wavef.nspins), ones(wavef.nkpts))
    @assert out.ngx == ngx

    # For test updating wave coefficients
    @testset "update" begin
        dst = joinpath(mktempdir(), "GaAs.check")
        cp(checkfile, dst)
        data = read_castep_check(dst)
        getcoeff = data -> data[:WAVEFUNCTION_DATA].coeffs 
        getcoeff(data)[1011] = 0.
        # Update the data
        CastepReader.update_wavefunction(dst, getcoeff(data))
        new_data = read_castep_check(dst)
        @test getcoeff(new_data)[1011] == 0.

        # Cannot update castep_bin files - should throw assertion error 
        binfile = joinpath(gaas, "GaAs.castep_bin")    
        dst = joinpath(mktempdir(), "GaAs.castep_bin")
        cp(binfile, dst)
        @test_throws AssertionError  CastepReader.update_wavefunction(dst, getcoeff(data))

        # Update with given wavefunction type
        wf = CastepReader.WaveFunction(data)
        dst = joinpath(mktempdir(), "GaAs.check")
        CastepReader.update_wavefunction(checkfile, dst, wf, data)
        new_data = read_castep_check(dst)
        @test getcoeff(new_data)[1011] == 0.
    end
end


@testset "pdos_bin" begin
    bin = joinpath(gaas, "GaAs.pdos_bin")
    output = read_pdos_bin(bin)
    @test size(output.kpoints)[2] == 10
    @test size(output.pdos_weights) == (20, 40, 10, 1)
    @test maximum(output.species) == 2

    m_mask = cr.m_channel_mask(output)
    @test m_mask[1] == Orbitals.dz2
    @test m_mask[6] == Orbitals.s
    @test m_mask[10] == Orbitals.s

    indices = cr.pdos_index_by_site(output)
    @test indices[1][Orbitals.dz2] == [1] 
    @test indices[1][Orbitals.s] == [6, 10] 
end

@testset "om" begin
    seed = joinpath(gaas, "GaAs")
    @testset "ome_bin" begin
        output = cr.read_ome_bin(seed)
        @test size(output)[4] == 10
    end


    @testset "cst_ome" begin
        output = cr.read_cst_ome(seed)
        @test size(output)[4] == 10
    end

    output = cr.read_om_castep(seed)
    @test size(output)[4] == 10
    
end

@testset "writers" begin
    @testset "chgcar" begin
        buffer = IOBuffer()
        CastepReader.write_chgcar(buffer, [rand(10, 10, 10)])
        seek(buffer, 0)
        content = take!(buffer)
        @test endswith(String(content), "\n")
    end
end

using FortranFiles
using Dates
using Unitful
using UnitfulAtomic

"Suffix of the file"
suffix(fname::AbstractString) = fname[findlast(isequal('.'), fname)+1:end]

"Part of the file exclude the suffix"
stem(fname::AbstractString) = fname[1:findlast(isequal('.'), fname)-1]

"Name of the file"
filename(fpath::AbstractString) = fpath[findlast(isequal('/'), fpath) + 1:end]

#const HaBohr2perAng = u"Eh_au" * u"a0_au" * u"a0_au" / u"Å" / u"ħ" * u"me_au"
#const HaBohr2perAng = u"Eh_au" * u"Å" * u"Å" / u"a0_au" / u"ħ" * u"me_au"

#Note that is an extra factor (ħ-1 me) attached to the matrix eleemnts (to make it momentum)
# For simplicity we carry it forward into the c0 term
const HaBohr2perAng = u"Eh_au" * u"a0_au" * u"a0_au" / u"Å" 
const HaBohr = u"Eh_au" * u"a0_au" 

"""
Read the cst_ome file with given number of bands, kpoints, and spins
"""
function read_cst_ome(fname, nb, nk, ns; endian="big-endian")
    om = Array{typeof((1.0 + im) * HaBohr2perAng)}(undef, nb, nb, 3, nk, ns)
    ffile = FortranFile(fname, convert=endian)
    for ik = 1:nk, is = 1:ns, c=1:3, ib = 1:nb, jb = 1:nb
        # Unit is Ha Bohr^2 / Ang
        @inbounds om[ib, jb, c, ik, is] = read(ffile, ComplexF64) * HaBohr2perAng
    end
    close(ffile)
    uconvert.(HaBohr, om)
end


function write_cst_ome(fname::AbstractString, om::Array{T}; endian="big-endian") where {T<:ComplexF64}
    nb, _, _, nk, ns = size(om)
    ffile = FortranFile(fname, "w", convert=endian)
    for ik = 1:nk, is = 1:ns, c=1:3, ib = 1:nb, jb = 1:nb
        # Unit is Ha Bohr^2 / Ang
        write(ffile, om[ib, jb, c, ik, is])
    end
    close(ffile)
end

function write_cst_ome(fname::AbstractString, om::Array{T}; endian="big-endian") where {T<:Quantity}
    nb, _, _, nk, ns = size(om)
    ffile = FortranFile(fname, "w", convert=endian)
    for ik = 1:nk, is = 1:ns, c=1:3, ib = 1:nb, jb = 1:nb
        # Unit is Ha Bohr^2 / Ang
        write(ffile, ustrip.(HaBohr2perAng, om[ib, jb, c, ik, is]))
    end
    close(ffile)
end

"""
    read_cst_ome(seed;)

Read the `<seed>.cst_ome` file, requires `<seed>.bands` file to be present.
"""
function read_cst_ome(seed)
    output = read_bands_castep_raw("$(seed).bands")
    ns, nk, nb = size(output.bands)
    ome = read_cst_ome("$(seed).cst_ome", nb, nk, ns)
    return ome
end

"""
    read_ome_bin(seed; legacy_format=false)

Read the `<seed>.ome_bin` file, requires `<seed>.bands` file to be present.
"""
function read_ome_bin(seed; legacy_format=false)
    output = read_bands_castep_raw("$(seed).bands")
    ns, nk, nb = size(output.bands)
    ome = read_ome_bin("$(seed).ome_bin", nb, nk, ns, legacy_format=legacy_format)
    return ome
end


"""
Read the ome_bin file
"""
function read_ome_bin(fname, nb, nk, ns; endian="big-endian", legacy_format=false)
    ffile = FortranFile(fname, convert=endian)
    version = read(ffile, Float64)
    supported = 1.0
    @assert version - supported < 0.001 "Unsupported .ome_bin file version"
    header = read(ffile, FString{80})

    # The ome_bin is written in the strange legacy format like the cst_ome file...
    if legacy_format == true
        om = Array{typeof((1.0 + im) * HaBohr2perAng)}(undef, nb, nb, 3, nk, ns)
        for ik = 1:nk, is = 1:ns, c = 1:3, ib = 1:nb, jb = 1:nb
            @inbounds om[ib, jb, c, ik, is] = read(ffile, ComplexF64) * HaBohr2perAng
        end
    else
        om = Array{typeof((1.0 + im) * HaBohr)}(undef, nb, nb, 3, nk, ns)
        for ik = 1:nk
            for is = 1:ns
                @inbounds om[:, :, :, ik, is] = read(ffile, (ComplexF64, nb, nb, 3)) * HaBohr
            end
        end
    end
    close(ffile)
    om
end


"""
    read_om_castep(ome_file, nb, nk, ns;endian="big-endian", legacy_format=false)

Read optical matrix element for CASTEP
"""
function read_om_castep(ome_file, nb, nk, ns;endian="big-endian", legacy_format=false)
    if endswith(ome_file, "cst_ome")
        out = read_cst_ome(ome_file, nb, nk, ns;endian)
    elseif endswith(ome_file, "ome_bin")
        out = read_ome_bin(ome_file, nb, nk, ns;endian, legacy_format)
    else
        throw(ErrorException("Unsupported file: $ome_file"))
    end
    out
end

"""
    read_om_castep(seed; endian="big-endian", legacy_format=false)

Read optical matrix element given a CASTEP seed name.
"""
function read_om_castep(seed; endian="big-endian", legacy_format=false)
    output = read_bands_castep_raw("$(seed).bands")
    if isfile("$(seed).ome_bin")
        ome_file = "$(seed).ome_bin"
    else
        ome_file = "$(seed).cst_ome"
    end
    ns, nk, nb = size(output.bands)
    read_om_castep(ome_file, nb, nk, ns;endian, legacy_format)
end

function write_ome_bin(fname::AbstractString, om::Array{T}; endian="big-endian", version=1.0) where {T<:ComplexF64}
    nb, _, _, nk, ns = size(om)
    ffile = FortranFile(fname, "w", convert=endian)
    write(ffile, version)
    write(ffile, FString(80, "Written by NLTools at $(now())"))
    tmp = om[:, :, :, 1, 1]
    for ik = 1:nk
        for is = 1:ns
            tmp .= om[:, :, :, ik, is]
            write(ffile, tmp)
        end
    end
    close(ffile)
end

function write_ome_bin(fname::AbstractString, om::Array{T}; endian="big-endian", version=1.0) where {T<:Quantity}
    nb, _, _, nk, ns = size(om)
    ffile = FortranFile(fname, "w", convert=endian)
    write(ffile, version)
    write(ffile, FString(80, "Written by NLTools at $(now())"))
    for ik = 1:nk
        for is = 1:ns
            # Ensure the unit written is HaBohr
            tmp = ustrip.(HaBohr, om[:, :, :, ik, is])
            write(ffile, tmp)
        end
    end
    close(ffile)
end

precompile(read_ome_bin, (String,))
precompile(read_cst_ome, (String,))

module CastepBin
using FortranFiles
CastepFortranFile(io) = FortranFile(io, convert="big-endian")
const IntF = Int32

"""
Macro for compact representation of how data should read from the binary file

For example, the following definition:
```julia
@readtag f outdict tags begin

    CELL_VERSION_NUMBER
    VERSION_NUMBER::Float64
    skip
    X::(Float64, :B)

end
```
expands into to:
```julia
if "CELL_VERSION_NUMBER" in keys(tags)
    gototag("CELL_VERSION_NUMBER", f, tags)
    output[:VERSION_NUMBER] = read(f, Float64)
    skiptag(f)
    output[:X] = read(f, outdict[:B])
end

Note that tags that are not present will simply be skipped.
```

"""
macro readtag(f, outdict, tags, expr)
    @assert expr.head == :block "Expect a block definition"
    eout = []
    thistag = Pair(nothing, [])

    function flush_if_block!(eout, thistag)
        if !isnothing(thistag.first)
            push!(eout, Expr(:if, 
            quote
                $(thistag.first) in keys($tags)
            end, 
            Expr(:block, thistag.second...))
            )
        end
    end

    for line in expr.args
        isa(line, LineNumberNode) && continue
        # This defines a new tags 
        # TODO - add condition to check if the tag exists or not 
        # Do nothing is the tags is not there!!!
        # Symbol as tag name
        if isa(line, Symbol)  
            if line == :skip || line == :SKIP
                push!(thistag.second, :(skiptag($f)))
                continue
            else
                line = String(line)
            end
        end
        # String - treat as the name of the required tag
        if isa(line, String)
            # Create the condition - only parse if the tag exists
            if !isnothing(thistag.first)
                flush_if_block!(eout, thistag)
            end
            thistag = Pair(line, [])
            push!(thistag.second, quote
                gototag($line, $f, $tags)     
            end)
            continue
        end
        # This is a read line
        if line.head == :(::)
            var = line.args[1]
            type = line.args[2]
            # Parse the existing parameters
            # Check if any parameters refer to exiting ones
            # X::(Float64, :A)
            # Expands to read(f, (Float64, outdict[:A]))
            if isa(type, Expr) && type.head == :tuple
                type_args = []
                # Construct the type field
                for arg in  type.args
                    if isa(arg, QuoteNode)
                        push!(type_args, quote
                            $outdict[$arg]    
                        end
                        )
                    else
                        push!(type_args, arg)
                    end
                end
                type = Expr(:tuple, type_args...)
            end
            push!(thistag.second, quote
                $outdict[$(Meta.quot(var))] = read($f, $type)
            end)
        end
    end

    # Flust out the last entry
    if !isnothing(thistag.first)
        flush_if_block!(eout, thistag)
    end

    esc(Expr(:block, eout...))
end

"""
Scan the binary file for tags
"""
function locate_tags(f::FortranFile;max_length=256)
    REC_LENGTH=4
    out = Dict{String, Int}()
    while true
        # All patterns has been found? just break the loop
        local rec
        try
            rec = FortranFiles.Record(f)
        catch y
            if isa(y, EOFError)
                break
            end
        end
        rec.subreclen <= max_length && record_tag!(out, rec;offset=REC_LENGTH)
        close(rec)
    end
    out
end

function locate_tags(name::AbstractString)
    open(name) do f
        locate_tags(CastepFortranFile(f))
    end
end

"""
Record the positions of eligible records
"""
function record_tag!(out::Dict{String,Int}, rec;offset::Int)
    loc = position(rec.io)
    fs = read(rec, FString{rec.subreclen})  # Read as string
    any(x -> x < 0, fs.data) && return nothing
    str = strip(convert(String, fs))
    m = match(r"[A-Z_0-9]*$", str)
    isnothing(m) && return
    if str in keys(out)
        str = str * "_SECOND"
    end
    out[str] = loc - offset
end

"Read in the unit cell information"
function read_cell!(output, f::FortranFile, tags::Dict{String,Int})

    @readtag f output tags begin
        BEGIN_UNIT_CELL_SECOND
        skip
        REAL_LATTICE::(Float64, 3, 3)
        skip
        RECIP_LATTICE::(Float64, 3, 3)
        skip
        VOLUME::Float64
        skip
        NUM_SPECIES::IntF
        skip
        NUM_IONS::IntF
        skip
        MAX_IONS_IN_SPECIES::IntF

        CELL_VERSION_NUMBER_SECOND
        VERSION_NUMBER::Float64

        "CELL%NUM_IONS_IN_SPECIES_SECOND"

        NUM_IONS_IN_SPECIES::(IntF, :NUM_SPECIES)

        "CELL%IONIC_POSITIONS_SECOND"

        IONIC_POSITIONS::(Float64, 3, :MAX_IONS_IN_SPECIES, :NUM_SPECIES)

        "CELL%IONIC_VELOCITIES_SECOND"

        IONIC_VELOCITIES::(Float64, 3, :MAX_IONS_IN_SPECIES, :NUM_SPECIES)

        "CELL%SPECIES_SYMBOL_SECOND"

        SPECIES_SYMBOL::(FString{8}, :NUM_SPECIES)
    end

    @readtag f output tags begin
        NKPTS_SECOND
        NKPTS::IntF

        KPOINTS_SECOND
        KPOINTS::(Float64, 3, :NKPTS)

        KPOINTS_WEIGHTS_SECOND
        KPOINTS_WEIGHTS::(Float64, :NKPTS)

    end

    output
end

function read_grid_properties!(output, f::FortranFile, tags::Dict{String, Int}, with_wavefunction=false)

    @readtag f output tags begin
        END_CELL_GLOBAL_SECOND
        FOUND_GROUND_STATE_WAVEFUNCTION::IntF
        FOUND_GROUND_STATE_DENSITY::IntF
        TOTAL_ENERGY::Float64
        FERMI_ENERGY::Float64
    end
    nbands, nspins = read(f, IntF, IntF)
    output[:NBANDS] = nbands
    output[:NSPINS] = nspins
    if with_wavefunction
        read_wavefunction_complex!(output, f, tags)
    end
    read_eigenvalue_and_occ!(output, f, tags)
    skiptag(f)
    ngxf, ngyf, ngzf = read(f, IntF, IntF, IntF)
    output[:NGX_FINE] = ngxf
    output[:NGY_FINE] = ngyf
    output[:NGZ_FINE] = ngzf
    read_charge_density!(output, f, tags)
end

"""
Read the grid properties but update the plane wave coefficients as given
"""
function update_wavefunction(output, f::FortranFile, tags::Dict{String, Int}, coeffs)

    @readtag f output tags begin
        END_CELL_GLOBAL_SECOND
        FOUND_GROUND_STATE_WAVEFUNCTION::IntF
        FOUND_GROUND_STATE_DENSITY::IntF
        TOTAL_ENERGY::Float64
        FERMI_ENERGY::Float64
    end
    nbands, nspins = read(f, IntF, IntF)
    output[:NBANDS] = nbands
    output[:NSPINS] = nspins
    # Write the wave function
    write_wavefunction_complex!(output, f, tags, coeffs)
end


"""
Read eigenvalues and occupations
"""
function read_eigenvalue_and_occ!(output, f, tags)
    nkpts::Int = output[:NKPTS]
    nspin::Int = output[:NSPINS]
    nbands::Int = output[:NBANDS]

    kpoints = zeros(3, nkpts)
    occ = zeros(nbands, nkpts, nspin)
    eignvalues = zeros(nbands, nkpts, nspin)
    for ik in 1:nkpts
        kpoints[:, ik] = read(f, (Float64, 3))
        for is in 1:nspin
            occ[:, ik, is] = read(f, (Float64, nbands))
            eignvalues[:, ik, is] = read(f, (Float64, nbands))
        end
    end
    output[:OCCUPATIONS] = occ
    output[:EIGENVALUES] = eignvalues
    output[:KPOINTS_OF_EIGENVALUES] = kpoints
end

"""
Read the charge density section from the binary file

NOTE: only support linear spin-polarisation for now
"""
function read_charge_density!(output, f, tags)
    ngxf::Int = output[:NGX_FINE]
    ngyf::Int = output[:NGY_FINE]
    ngzf::Int = output[:NGZ_FINE]
    nspin::Int = output[:NSPINS]
    # NOTE - assume not NCM
    spin_density = zeros(ComplexF64, ngxf, ngyf, ngzf)
    charge_density = zeros(ComplexF64, ngxf, ngyf, ngzf)

    for _ in 1:ngxf
        for _ in 1:ngyf
            if nspin == 2
                nx, ny, zcol, spincol = read(f, IntF, IntF, (ComplexF64, ngzf), (ComplexF64, ngzf))
                charge_density[nx, ny, :] = zcol
                spin_density[nx, ny, :] = spincol
            else
                nx, ny, zcol = read(f, IntF, IntF, (ComplexF64, ngzf))
                charge_density[nx, ny, :] = zcol
            end
        end
    end
    if nspin == 2
        output[:SPIN_DENSITY] = spin_density
    end
    output[:CHARGE_DENSITY] = charge_density
end

function read_wavefunction_complex!(output, f, tags)
    header = read(f, FString{4})
    ngx, ngy, ngz = read(f, IntF, IntF, IntF)
    coeff_size_1, spinorcomps, nbands_max, nkpts, nspins = read(f, IntF, IntF, IntF, IntF, IntF) 
    coeffs = zeros(ComplexF64, coeff_size_1, spinorcomps, nbands_max, nkpts, nspins)
    nwaves_at_kp = zeros(Int, nkpts)
    kpts = zeros(3, nkpts)
    pw_grid_coord = zeros(IntF, 3, coeff_size_1, nkpts)

    for is in 1:nspins
        for ik in 1:nkpts
            kpts[:, ik], nwaves = read(f, (Float64, 3), IntF)
            nwaves_at_kp[ik] = nwaves

            # Grid coordinates
            coords_x = read(f, (IntF, nwaves))
            coords_y = read(f, (IntF, nwaves))
            coords_z = read(f, (IntF, nwaves))
            pw_grid_coord[1, 1:nwaves, ik] = coords_x 
            pw_grid_coord[2, 1:nwaves, ik] = coords_y 
            pw_grid_coord[3, 1:nwaves, ik] = coords_z 

            for ib in 1:nbands_max
                for ispinor in 1:spinorcomps
                    coeff = read(f, (ComplexF64, nwaves))
                    coeffs[1:nwaves, ispinor, ib, ik, is] = coeff
                end
            end
        end
    end
    data = (
        ngx=ngx,
        ngy=ngy,
        ngz=ngz,
        pw_grid_coord=pw_grid_coord,
        coeffs=coeffs,
        kpts=kpts,
        nwaves_at_kp=nwaves_at_kp,
    )
    output[:WAVEFUNCTION_DATA] = data
end

"""
Write the wavefunction coefficients back to the file
"""
function write_wavefunction_complex!(output, f, tags, coeffs)
    header = read(f, FString{4})
    ngx, ngy, ngz = read(f, IntF, IntF, IntF)
    coeff_size_1, spinorcomps, nbands_max, nkpts, nspins = read(f, IntF, IntF, IntF, IntF, IntF) 
    nwaves_at_kp = zeros(Int, nkpts)
    kpts = zeros(3, nkpts)
    pw_grid_coord = zeros(IntF, 3, coeff_size_1, nkpts)

    for is in 1:nspins
        for ik in 1:nkpts
            kpts[:, ik], nwaves = read(f, (Float64, 3), IntF)
            nwaves_at_kp[ik] = nwaves

            # Grid coordinates
            coords_x = read(f, (IntF, nwaves))
            coords_y = read(f, (IntF, nwaves))
            coords_z = read(f, (IntF, nwaves))
            pw_grid_coord[1, 1:nwaves, ik] = coords_x 
            pw_grid_coord[2, 1:nwaves, ik] = coords_y 
            pw_grid_coord[3, 1:nwaves, ik] = coords_z 

            for ib in 1:nbands_max
                for ispinor in 1:spinorcomps
                    write(f, coeffs[1:nwaves, ispinor, ib, ik, is])
                end
            end
        end
    end
    data = (
        ngx=ngx,
        ngy=ngy,
        ngz=ngz,
        pw_grid_coord=pw_grid_coord,
        coeffs=coeffs,
        kpts=kpts,
        nwaves_at_kp=nwaves_at_kp,
    )
    output[:WAVEFUNCTION_DATA] = data
end

"""
Read the forces
"""
function read_forces!(output, f::FortranFile, tags::Dict{String,Int})
    @readtag f output tags begin
        FORCES
        FORCES::(Float64, 3, :MAX_IONS_IN_SPECIES, :NUM_SPECIES)
    end
    output
end


function read_stress!(output, f::FortranFile, tags::Dict{String,Int})
    @readtag f output tags begin
        STRESS
        STRESS::(Float64, 6)
        STRAIN::(Float64, 3,3)
    end
end


"""
Seek to a tag and be read to read in the data following it
"""
function gototag(tag, f::FortranFile, tags::Dict{String,Int})
    seek(f.io, tags[string(tag)])
    rec = FortranFiles.Record(f)
    close(rec)
end

"""Skip to the next tag"""
function skiptag(f)
    rec = FortranFiles.Record(f)
    close(rec)
end

"""
READ RAW data from CASTEP binary file

Expect everyting to be in the atomic units
Atoms are indexed by (species, index in specie) 
"""
function read_castep_check(fname)
    output = Dict{Symbol, Any}()
    open(fname) do fhandle
        f = CastepFortranFile(fhandle)
        tags = locate_tags(f)
        read_cell!(output, f, tags)
        read_forces!(output, f, tags)
        read_stress!(output, f, tags)
        has_wave = !("CASTEP_BIN" in keys(tags))
        read_grid_properties!(output, f, tags, has_wave)
    end
    output
end


"""
    update_wavefunction(fname, coeffs::AbstractArray{T, 5}) where {T}

Update the wavefunction (coefficients) stored in the check file. 
The coefficients array passed should be modified based on that was originally
read from the file.

* **NOTE** this function will overwrite target file.
"""
function update_wavefunction(fname, coeffs::AbstractArray{T, 5}) where {T}
    output = Dict{Symbol, Any}()
    open(fname, "r+") do fhandle
        f = CastepFortranFile(fhandle)
        tags = locate_tags(f)
        has_wave = !("CASTEP_BIN" in keys(tags))
        @assert has_wave "The file does not contain any wavefunction data"
        update_wavefunction(output, f, tags, coeffs)
    end
end


"""
    update_wavefunction(src, dst, coeffs::AbstractArray{T, 5}) where {T}

Update the wavefunction (coefficients) stored in the check file. 
The coefficients array passed should be modified based on that was originally
read from the file.
"""
function update_wavefunction(src, dst, coeffs::AbstractArray{T, 5}) where {T}
    cp(src, dst)
    update_wavefunction(dst, coeffs)
end

precompile(read_castep_check, (String, ))
export read_castep_check, update_wavefunction

end

using .CastepBin
import .CastepBin
using FortranFiles

"""
Locate the tag and read in data

```julia
@readtag f outdict tags begin

    CELL_VERSION_NUMBER
    VERSION_NUMBER::Float64
    X::(Float64, :B)

end
```

translates to
```julia
gototag("CELL_VERSION_NUMBER", f, tags)
output[:VERSION_NUMBER] = read(f, Float64)
output[:X] = read(f, outdict[:B])
```
"""
macro readtag(f, outdict, tags, expr)
    @assert expr.head == :block "Expect a block definition"
    eout = []
    for line in expr.args
        isa(line, LineNumberNode) && continue
        # This defines a new tags 
        # TODO - add condition to check if the tag exists or not 
        # Do nothing is the tags is not there!!!
        if isa(line, Symbol)  
            if line == :skip
                push!(eout, :(skiptag($f)))
            else
                push!(eout, quote
                    gototag($(Meta.quot(line)), $f, $tags)     
                end
                )
            end
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
            push!(eout, quote
                $outdict[$(Meta.quot(var))] = read($f, $type)
            end)
        end
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
function record_tag!(out, rec;offset)
    loc = position(rec.io)
    fs = read(rec, FString{rec.subreclen})  # Read as string
    any(x -> x < 0, fs.data) && return nothing
    str = strip(convert(String, fs))
    m = match(r"[A-Z_0-9]*$", str)
    isnothing(m) || (out[str] = loc - offset)
end

"Read in the unit cell information"
function read_cell!(output, f::FortranFile, tags)
    seek(f.io, tags["BEGIN_UNIT_CELL"])
    read(f, FString{256})
    read(f, FString{256})
    output[:REAL_LATTICE] = read(f, (Float64, 3, 3))
    read(f, FString{256})
    output[:RECIP_LATTICE] = read(f, (Float64, 3, 3))
    read(f, FString{256})
    output[:VOLUME] = read(f, Float64)

    read(f, FString{256})
    nspecies = read(f, IntF)
    output[:NUM_SPCIES] = nspecies

    read(f, FString{256})
    output[:NUM_IONS] = read(f, IntF)
    
    read(f, FString{256})
    max_ions = read(f, IntF)
    output[:MAX_IONS_IN_SPECIES] = max_ions

    # The following sections may not be ordered or may not always be present
    # So we read by direct seeking
    gototag("CELL_VERSION_NUMBER", f, tags)
    output[:VERSION_NUMBER] = read(f, Float64)
    gototag("CELL%NUM_IONS_IN_SPECIES", f, tags)
    output[:NUM_IONS_IN_SPECIES] = read(f, (IntF, nspecies))
    gototag("CELL%IONIC_POSITIONS", f, tags)
    output[:IONIC_POSITIONS] = read(f, (Float64, 3, max_ions, nspecies))
    gototag("CELL%IONIC_VELOCITIES", f, tags)
    output[:IONIC_VELOCITIES] = read(f, (Float64, 3, max_ions, nspecies))
    output
end

"""
Read the forces
"""
function read_forces!_(output, f::FortranFile, tags)
    gototag("FORCES", f, tags)
    nspecies = output[:NUM_SPCIES]
    nmax = output[:MAX_IONS_IN_SPECIES]
    output[:FORCES] = read(f, (Float64, 3, nmax, nspecies))
    output
end

"""
Read the forces
"""
function read_forces!(output, f::FortranFile, tags)

    @readtag f output tags begin
        FORCES
        FORCES::(Float64, 3, :MAX_IONS_IN_SPECIES, :NUM_SPCIES)
    end
    output
end


function read_stress!(output, f::FortranFile, tags)
    gototag("STRESS", f, tags)
    output[:STRESS] = read(f, (Float64, 6))   # Components of the stress tensor (Stress)
    output[:STRAIN] = read(f, (Float64, 3,3))  # Components of the stress tensor (External)
    output
end


"""
Seek to a tag and be read to read in the data following it
"""
function gototag(tag, f::FortranFile, tags)
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
    end
    output
end

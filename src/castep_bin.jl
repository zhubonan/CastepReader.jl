
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
    isnothing(m) || (out[str] = loc - offset)
end

"Read in the unit cell information"
function read_cell!(output, f::FortranFile, tags::Dict{String,Int})

    @readtag f output tags begin
        BEGIN_UNIT_CELL
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

        CELL_VERSION_NUMBER
        VERSION_NUMBER::Float64

        "CELL%NUM_IONS_IN_SPECIES"
        NUM_IONS_IN_SPECIES::(Float64, :NUM_SPECIES)

        "CELL%IONIC_POSITIONS"
        NUM_IONS_IN_SPECIES::(Float64, 3, :MAX_IONS_IN_SPECIES, :NUM_SPECIES)

        "CELL%IONIC_VELOCITIES"
        IONIC_VELOCITIES::(Float64, 3, :MAX_IONS_IN_SPECIES, :NUM_SPECIES)
    end
    output
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
    end
    output
end

precompile(read_castep_check, (String, ))
export read_castep_check

end

using .CastepBin
#=
Reader for cell/param files
=#
using Unitful
using UnitfulAtomic

"Clean lines by removing new line symbols and comments"
function clean_cell_lines(lines_in::Vector{String})
    lines_out = String[]
    for line in lines_in
        tmp = strip(line)
        length(tmp) == 0 && continue
        # Location of the comment
        matched = match(r"[#!]", tmp, 1)
        # No comment - keep the whole line
        if matched === nothing
            push!(lines_out, string(tmp))
        elseif matched.offset == 1
            continue
        else
            # Remove the trailing part
            push!(lines_out, string(strip(tmp[1:matched.offset-1])))
        end
    end
    lines_out
end

"Read in CASTEP's cell/param files"
function read_cell_castep(fname)
    f = open(fname)
    lines = clean_cell_lines(readlines(f))
    close(f)
    kws = Dict{Symbol, String}()
    blocks = Dict{Symbol, Vector{String}}()
    inblock = false
    block_name = ""
    current_block = String[]
    for il in eachindex(lines)
        line = lowercase(lines[il])
        # Empty line - do nothing
        length(il) == 0 && continue

        # Are we hitting a block
        mblock = match(r"%block (\w+)", line)
        if mblock !== nothing
            block_name = mblock.captures[1]
            inblock = true
            continue
        end

        mblock = match(r"%endblock (\w+)", line)
        if mblock !== nothing
            block_name = mblock.captures[1]
            inblock = false
            # Record the block
            blocks[Symbol(block_name)] = current_block
            # Reset current block
            current_block = String[]
            block_name = ""
            continue
        end

        # In inside the block - push the lines
        if inblock == true
            push!(current_block, string(strip(line)))
            continue
        end

        # Process keyword argumetns
        tokens = split(line, r"[:= ]+", limit=2)  
        if length(tokens) == 1
            kws[Symbol(tokens[1])] = ""
        else
            kws[Symbol(tokens[1])] = string(tokens[2])
        end
    end
    CastepInFile(kws, blocks)
end

"Type for CASTEP input file"
struct CastepInFile
    kws::Dict
    blocks::Dict
end


function Base.show(io::IO, c::CastepInFile)
    println("CastepInFile with keys: ", join(keys(c.kws), ", "))
    println("Blocks: ", join(keys(c.blocks), ", "))
end


"Allow access fields by c[:xx] notation"
function Base.getindex(c::CastepInFile, s::Symbol)
    if (s in keys(c.kws)) & (s in keys(c.blocks))
        throw(ErrorException("Key $(s) exists in both key and blocks"))
    end
    if s in keys(c.kws)
        return c.kws[s]
    end
    if s in keys(c.kws)
        return c.blocks[s]
    end
    throw(KeyError(s))
end

"""Read out the lattice from the input file"""
function read_lattice(ci::CastepInFile)
    cellmat = zeros(Float64, 3, 3)
    unit = ""
    if :lattice_cart in keys(ci.blocks)
        i = 1
        for line in ci.blocks[:lattice_cart]
            tokens = split(line)
            if length(tokens) == 1
                unit = tokens[1]
                continue
            end
            cellmat[:, i] = [parse(Float64, s) for s in tokens]
            i += 1
        end
    elseif :lattice_abc in keys(ci.blocks)
        a, b, c = [parse(Float64, x) for x in split(ci.blocks[:lattice_abc][1])]
        α, β, γ = [parse(Float64, x) for x in split(ci.blocks[:lattice_abc][2])]
        cellmat = cellpar2mat(a, b, c, α, β, γ, degree=true)
    else
        throw(ErrorException("No lattice related block is found"))
    end
    cellmat * u"Å"
end


"""
Read in the symmetry operations

Returns matrix of the rotations and that for the displacements
"""
function read_symm_disp(ci::CastepInFile)
    lines = ci.blocks[:symmetry_ops]
    vecs = []
    for line in lines
        push!(vecs, [parse(Float64, x) for x in split(line)])
    end
    nvec = length(vecs)
    data = reshape(hcat(vecs...), 3, 4, div(nvec, 4))
    # Separate rotational symmetries and the displacement
    data[:, 1:3, :], data[:, 4, :]
end
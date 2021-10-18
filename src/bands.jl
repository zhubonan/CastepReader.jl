using Unitful
using UnitfulAtomic

"""
Read the bands file

The bands file contains the eigenvalues of the 

# Arguments

- reorder: If true, recover the orignial kpoint orders as defined in the cell file. Note that the 
order of kpoint in other quantities, such as the optical matrix, are the same as they appear in the
bands file rather than following the internal "counter". 

"""
function read_bands_castep_raw(fname::AbstractString)
    nk, ns, nelect, neign, fermi_engs, cell_mat = read_bands_header(fname)
    # Allocate the ararys
    bands = zeros(Float64, ns, nk, neign)
    kweight = zeros(nk)
    kpoints = zeros(3, nk)

    bfile = open(fname)
    ecount = 1
    is = 0
    ik = 0
    kcount = 1
    # Record the "original" kpoint index, as they defined in the input cell file
    kidx = zeros(Int64, nk)
    for (i, line) in enumerate(eachline(bfile))
        # Skip to line 11
        i < 10 && continue
        if contains(line, "K-point") 
            tokens = split(line)
            ik = parse(Int, tokens[2])
            kidx[kcount] = ik
            kvec = [parse(Float64, token) for token in tokens[3:5]]
            weight = parse(Float64, tokens[6])
            kpoints[:, kcount] = kvec
            kweight[kcount] = weight
            kcount += 1
            continue
        end
        if contains(line, "Spin component")
            is = parse(Int, split(line)[3])  # Covnert unit to eV
            ecount = 1
            continue
        end
        # Read the eigen values for the corresponding kpoint
        # Because kcount was already incremeted - the current kpoint is there kcount - 1
        bands[is, kcount - 1, ecount] = parse(Float64, strip(line))
        ecount += 1
    end

    # Convert unit to eV
    bands_with_unit = uconvert.(u"eV", bands .* u"Eh_au")


    return (kpoints=kpoints, 
            kweights=kweight, 
            kidx=kidx, 
            bands=bands_with_unit, 
            nelect=nelect, 
            fermi_engs=fermi_engs, 
            cell_mat=cell_mat)
end


"""
Read the header for the bands file
"""
function read_bands_header(fname::AbstractString)
    bfile = open(fname)
    nk = 0
    ns = 0
    nelect = 0.
    neigen = 0
    efermi = Float64[]
    for ln in eachline(bfile)
        if contains(ln, "Number of k-points")
            nk = parse(Int, split(ln)[end])
            continue
        end

        if contains(ln, "Number of spin components")
            ns = parse(Int, split(ln)[end])
            continue
        end

        if contains(ln, "Number of electrons")
            nelect = parse(Float64, split(ln)[end])
            continue
        end

        if contains(ln, "Number of eigenvalues")
            neigen = parse(Int, split(ln)[end])
            continue
        end
        if contains(ln, "Fermi ")
            tokens = split(ln)
            efermi = [uconvert(u"eV", parse(Float64, tokens[idx]) * u"Eh_au") for idx = 6:length(tokens)]
            continue
        end
        if contains(ln, "Unit cell vectors")
            break
        end
    end
    # Read the cell vectors
    cell_vec = zeros(Float64, 3, 3)
    cell_vec[:, 1] = Float64[parse(Float64, token) for token in split(readline(bfile))]
    cell_vec[:, 2] = Float64[parse(Float64, token) for token in split(readline(bfile))]
    cell_vec[:, 3] = Float64[parse(Float64, token) for token in split(readline(bfile))]
    cell_vec_with_unit = uconvert.(u"Å", cell_vec * u"a0_au")

    close(bfile)
    @assert nk > 0
    @assert ns > 0
    @assert nelect > 0
    @assert neigen > 0
    @assert length(efermi) > 0
    return nk, ns, nelect, neigen, efermi, cell_vec_with_unit
end


precompile(read_bands_castep_raw, (String,))
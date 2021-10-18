using Parameters


"""
    read_pdos_bin(f::FortranFile)

Read the `pdos_bin` file.
"""
function read_pdos_bin(f::FortranFile)

    fversion = read(f, Float64)
    fheader = strip(convert(String, read(f, FString{80})))

    nk = read(f, IntF)
    ns = read(f, IntF)
    num_popn_orb = read(f, IntF)
    max_eignenv = read(f, IntF)

    # Now start the read the actual data
    species = read(f, (IntF, num_popn_orb))
    ion = read(f, (IntF, num_popn_orb))
    am_channel = read(f, (IntF, num_popn_orb))

    # Allocate storage space
    pdos_weights = zeros(Float64, (num_popn_orb, max_eignenv, nk, ns))
    kpos = zeros(Float64, (3, nk))
    num_eignvalues = zeros(Int, ns)

    for ik in 1:nk
        _, kpos[:, ik] = read(f, IntF, (Float64, 3)) 
        for is in 1:ns
            read(f, IntF)
            num_eignvalues[is] = read(f, IntF)
            for ib in 1:num_eignvalues[is]
                pdos_weights[:, ib, ik, is] = read(f, (Float64, num_popn_orb))
            end
        end
    end
    return (
        fversion=fversion,
        fheader=fheader,
        nk=nk,
        ns=ns,
        norb=num_popn_orb,
        max_eignenv=max_eignenv,
        species=species,
        ion=ion,
        am_channel=am_channel,
        pdos_weights=pdos_weights,
        kpoints=kpos,
        num_eignvalues=num_eignvalues,
    )
end

"""
    read_pdos_bin(fname::AbstractString)

Read the `pdos_bin` file.
"""
function read_pdos_bin(fname::AbstractString)
    open(fname) do fh
        f = CastepFortranFile(fh)
        read_pdos_bin(f)
    end
end

precompile(read_pdos_bin, (String,))
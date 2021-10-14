using Parameters

"""
Type storing the information read from the PDOS file
"""
@with_kw struct PDOS
    fversion::Float64
    fheader::String
    nk::Int
    ns::Int
    norb::Int
    max_eignenv::Int
    "Index array for species"
    species::Vector{Int}
    "Index array for ion number in the species"
    ion::Vector{Int}
    "Index array for the angumar momentum channel"
    am_channel::Vector{Int}
    pdos_weights::Array{Float64, 4}
    kpoints::Array{Float64, 2}
    num_eignvalues::Int
end

"""
Read the pdos_bin file
"""
function read_pdos_bin(f::FortranFile)

    fversion = read(f, Float64)
    fheader = strip(string(read(f, FString{80})))

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
        kpos[:, ik] = read(f, IntF, (Float64, 3)) 
        for is in 1:ns
            read(f, IntF)
            num_eignvalues[is] = read(f, IntF)
            for ib in 1:num_eignvalues[is]
                pdos_weights[:, ib, ik, is] = read(f, (Float64, num_popn_orb))
            end
        end
    end
    PDOS(
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

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
        species=species,   # Indexing array for species
        ion=ion,   # Indexing array for the ion numbers
        am_channel=am_channel,  # Indexing array for the angular momentum channels (l)
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

"""
Reorder the PDOS information

Return a dictionary which can be used as

```julia
output = pdos_index_by_site(pdos_data)
idx = output[1][Orbitals.s]   # Index of the s orbital in the pdos_data.weights
weights = pdos_data.pdos_weights[idx, :, :, :]  # Slice the weight array to obtain the weigths for this site/am combination
```
"""
function pdos_index_by_site(parsed_items)
    mapping = Orbitals.castep_orbital_order
    unique_species = sort(unique(parsed_items.species))
    site_index = 1
    output_data = Dict{Int, Dict{Orbital, Vector{Int}}}()
    for ns in unique_species
        specie_mask = parsed_items.species .== ns
        total_ions = maximum(parsed_items.ion[specie_mask])  # Total number of ions for this species
        # Treat each ion for this speice
        for nion in 1:total_ions
            ion_mask = (parsed_items.ion .== nion) .& specie_mask
            max_am = maximum(parsed_items.am_channel[ion_mask])
            orbs = Dict{Orbital, Vector{Int}}()
            # Angular momentum channel for each site
            for am in 0:max_am
                ion_am_mask = (parsed_items.am_channel .== am) .& ion_mask
                ion_am_idx = findall(ion_am_mask)
                for (iam, iloc) in enumerate(ion_am_idx)
                    #iloc is the index of the orbital
                    #iam is the m channel
                    this_orb = mapping[am+1][((iam - 1) % (2 * am + 1)) + 1]  # m number of this orbital
                    if this_orb in keys(orbs)
                        push!(orbs[this_orb], iloc)
                    else
                        orbs[this_orb] = [iloc]
                end
            end
        end
        output_data[site_index] = orbs 
        site_index += 1
        end
    end
    output_data
end

"""
Get an masking for the m channels inferred from the repeated appearance of l channels

CASTEP does not write the m channel labels explicity, but they can be inferred from repeated l channels.
This assumes that l channels of the same atom is written consecutively, which seem to be the case.
"""
function m_channel_mask(pdos_data)
    mapping = Orbitals.castep_orbital_order
    norb = pdos_data.norb
    last_am = pdos_data.am_channel[1]
    current_m = 0
    m_channels = Array{Orbital, 1}(undef, norb)
    m_channels[1] = mapping[last_am + 1][1]
    for i in 2:norb
        current_am = pdos_data.am_channel[i]
        if current_am == last_am
            current_m += 1
        else
            current_m = 0
        end
        m_channels[i] = mapping[current_am+1][(current_m % (2 * current_am + 1)) + 1]
        last_am = current_am
    end
    m_channels
end

precompile(read_pdos_bin, (String,))

export m_channel_mask, pdos_index_by_site
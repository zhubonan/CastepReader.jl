
module CastepReader



    include("castep_bin.jl")
    include("optical_matrix.jl")
    include("pdos_bin.jl")
    include("bands.jl")
    include("cell.jl")

    export read_castep_check, read_cst_ome, read_ome_bin, read_pdos_bin
end
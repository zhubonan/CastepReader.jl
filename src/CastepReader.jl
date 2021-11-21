
module CastepReader

    CastepFortranFile(io) = FortranFile(io, convert="big-endian")

    const IntF = Int32

    """
    Module for electronic Enum types
    """
    module Orbitals
        @enum Spin spinup=1 spindown=2
        @enum Orbital begin
            s=0
            px=10
            py=11
            pz=12
            dxy=20
            dyz=21
            dz2=22
            dxz=23
            dx2=24
            f_xxx=30
            f_yyy=31
            f_zzz=32
            f_xyz=33
            f_z_xx_yy=34
            f_y_zz_xx=35
            f_x_yy_zz=36
        end
        function l(x::Orbital)
            div(Integer(x), 10) 
        end
        function m(x::Orbital)
            Integer(x) % 10
        end

        const _castep_orb_names = Dict{Orbital, String}(
            s => "S",
            px => "Px",
            py => "Py",
            pz => "Pz",
            dxy => "Dxy",
            dyz => "Dzy",
            dz2 => "Dzz",
            dxz => "Dzx",
            dx2 => "Dxx-yy",
            f_xxx => "Fxxx",
            f_yyy => "Fyyy",
            f_zzz => "Fzzz",
            f_xyz => "Fxyz",
            f_z_xx_yy => "Fz(xx-yy)",
            f_y_zz_xx => "Fy(zz-xx)",
            f_x_yy_zz => "Fx(yy-zz)",
        )
        function castep_name(orb::Orbital)
            _castep_orb_names[orb]
        end
        const castep_orbital_order = [
            [Orbitals.s],
            [Orbitals.px, Orbitals.py, Orbitals.pz],
            [Orbitals.dz2, Orbitals.dyz, Orbitals.dxz, Orbitals.dx2, Orbitals.dxy],
            [Orbitals.f_xxx, Orbitals.f_yyy, Orbitals.f_zzz, Orbitals.f_xyz, Orbitals.f_z_xx_yy, Orbitals.f_y_zz_xx, Orbitals.f_x_yy_zz]
        ]
        export Spin, Orbital, l, m, castep_name
    end

    using .Orbitals
    export Orbitals

    include("castep_bin.jl")
    include("wave.jl")
    include("optical_matrix.jl")
    include("pdos_bin.jl")
    include("bands.jl")
    include("cell.jl")

    export read_castep_check, read_cst_ome, read_ome_bin, read_pdos_bin
end
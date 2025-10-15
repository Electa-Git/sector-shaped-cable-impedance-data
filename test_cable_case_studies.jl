
#NB:: For this file to work as of now you need to dev a fork of LineCableModels where 
# the Sector-Shaped data models and the OpenDSS formulations are implemented 
# dev ->  https://github.com/MohamedNumair/LineCableModels.jl/tree/integration/sector-fem-dss
# make sure to checkout the branch `integration/sector-fem-dss`

using LineCableModels
using LineCableModels.DataModel
using LineCableModels.Engine.FEM

using DataFrames
import LineCableModels.BackendHandler: renderfig #hide
fullfile(filename) = joinpath(@__DIR__, filename); #hide
set_verbosity!(3); #hide

using FileIO
using JLD2
using GLMakie

function compute_dss_formulation!(problem, results_df, internal_impedance, earth_impedance)
    dss_formulation = FormulationSet(:DSS,
        internal_impedance=internal_impedance,
        earth_impedance=earth_impedance,
        options=LineCableModels.Engine.DSSOptions()
    )

    time_taken = @elapsed begin
        _, line_params = compute!(problem, dss_formulation)
    end

    Internal_form = split(split(string(dss_formulation.internal_impedance),".")[end],"(")[1]
    Earth_form = split(split(string(dss_formulation.earth_impedance),".")[end],"(")[1]
    formulation_name = Earth_form == Internal_form ? Internal_form : "$(Internal_form)_$(Earth_form)"
    
    @info "Computed line parameters using DSS with $(Internal_form) internal and $(Earth_form) earth formulations."
    push!(results_df, (problem.system.system_id, formulation_name, time_taken, line_params.Z[:,:,1].*1000, line_params.Y[:,:,1].*1000, "", ""))
end

function compute_fem_formulation!(problem, results_df, earth, materials)
    f = problem.frequencies
    system = problem.system
    domain_radius = calc_domain_size(earth, f)
    mesh_transition = MeshTransition(system,
        [1],
        r_min=0.09,
        r_length=0.25,
        mesh_factor_min=0.01 / (domain_radius / 5),
        mesh_factor_max=0.25 / (domain_radius / 5),
        n_regions=3)

    opts = (force_remesh=true,
        force_overwrite=true,
        plot_field_maps=false,
        mesh_only=false, save_path="C:\\LCMFEMM\\fem_output_$(system.system_id)",
        keep_run_files=true,
        verbosity=0)

    fem_formulation = FormulationSet(:FEM,
        impedance=Darwin(),
        admittance=Electrodynamics(),
        domain_radius=domain_radius,
        domain_radius_inf=domain_radius * 1.25,
        elements_per_length_conductor=3,
        elements_per_length_insulator=3,
        elements_per_length_semicon=1,
        elements_per_length_interfaces=3,
        points_per_circumference=10,
        mesh_size_min=1e-6,
        mesh_size_max=domain_radius / 5,
        mesh_transitions=[mesh_transition],
        mesh_size_default=domain_radius / 10,
        mesh_algorithm=5,
        mesh_max_retries=20,
        materials=materials,
        options=opts,
    )

    time_taken = @elapsed begin
        _, line_params = compute!(problem, fem_formulation)
    end
    @info "Computed line parameters using FEM with Darwin internal and Electrodynamics earth formulations."
    push!(results_df, (system.system_id, "FEM", time_taken, line_params.Z[:,:,1].*1000, line_params.Y[:,:,1].*1000, "", ""))
end


function get_krn(zhat::Matrix{ComplexF64}; nphase::Int=3, ncct::Int=4)

    # Partition the matrix based on Kersting's equation 4.46
    zij = zhat[1:nphase, 1:nphase]
    zin = zhat[1:nphase, nphase+1:ncct]
    znj = zhat[nphase+1:ncct, 1:nphase]
    znn = zhat[nphase+1:ncct, nphase+1:ncct]

    # Perform Kron reduction as per Kersting's equation 4.55
    zabc = zij - zin * inv(znn) * znj

    # Calculate the neutral transformation matrix from Kersting's equation 4.52
    #tn = -inv(znn) * znj

    return zabc
end

# Schur Complement https://en.wikipedia.org/wiki/Schur_complement


function get_z012(m_phase::Matrix)
    # Check if the matrix is 3x3
    if size(m_phase) != (3, 3)
        throw(DimensionMismatch("Input matrix must be 3x3."))
    end

    # Define the complex operator 'a' (1 angle 120 degrees)
    α = exp(im * 2 * pi / 3) # cis(120°) or -0.5 + im*sqrt(3)/2

    # Define the Fortescue transformation matrix T
        T = ComplexF64[
        1 1 1;
        1 α^2 α;
        1 α α^2
    ]

    # Define the inverse of the Fortescue transformation matrix T_inv
    # T_inv = (1/3) * T' conjugate transpose (Hermitian transpose)
    # Or explicitly:
    T_inv = (1 / 3) * ComplexF64[
        1 1 1;
        1 α α^2;
        1 α^2 α
    ]

    # Calculate the sequence matrix: M_seq = T_inv * M_phase * T
    m_seq = T_inv * m_phase * T

    # Extract the diagonal components
    zero_sequence = m_seq[1, 1]
    positive_sequence = m_seq[2, 2]
    negative_sequence = m_seq[3, 3]

    #return m_seq, zero_sequence, positive_sequence, negative_sequence
    println("===================================")
    println("Zero Sequence z⁰: $zero_sequence")
    println("Positive Sequence z⁺: $positive_sequence")
    println("Negative Sequence z⁻: $negative_sequence")
    println("===================================")
    return m_seq
end


####################################################
################## Materials  ########################
####################################################
materials = MaterialsLibrary(add_defaults=true)
pvc = Material(Inf, 8.0, 1.0, 20.0, 0.1) # simple PVC
add!(materials, "pvc", pvc)

delete!(materials, "aluminum")
delete!(materials, "copper")
R_dc_nom = 0.32 # ohm/km
L_m = 1000.0 # m
A_sector = 92.74893091060221 * 1e-6 # m²
rho_from_Rdc = (R_dc_nom / L_m) * A_sector
aluminum = Material(rho_from_Rdc, 1.0, 1.0, 20.0, 0.00403) # Aluminum with corrected resistivity from IEC 60287-3-2
add!(materials, "aluminum", aluminum)
r_strand = 0.79e-3
rho_neutral_from_Rdc = (0.32 / 1000) * (30 * (pi * r_strand^2))
copper = Material(rho_neutral_from_Rdc, 1.0, 1.0, 20.0, 0.00393) # Copper with corrected resistivity from IEC 60287-3-2
add!(materials, "copper", copper)
materials

FileIO.save("materials.jld2","materials",materials)



## PREPARE DATAFRAME TO STORE RESULTS

results_df = DataFrame(CableID=String[], Formulation=String[], time_taken = Float64[], line_params_Z = Any[], line_params_Y = Any[], Jresiduals=Any[], SE_kpi=Any[])





####################################################
################## NAYY 95  ########################
####################################################
##
    n_sectors = 3
    r_back_mm = 10.24       # sector radius b
    d_sector_mm = 9.14        # sector depth s
    r_corner_mm = 1.02      # corner radius c
    theta_cond_deg = 119.0    # sector angle φ
    theta_cond_deg = 120.0    # sector angle φ
    ins_thick = 1.1e-3        # core insulation thickness

    sector_params = SectorParams(
        n_sectors,
        r_back_mm / 1000,
        d_sector_mm / 1000,
        r_corner_mm / 1000,
        theta_cond_deg,
        ins_thick
    )


    rot_angles = (0.0, 120.0, 240.0)
    sectors = [Sector(sector_params, ang, aluminum) for ang in rot_angles]
    insulators = [SectorInsulator(sectors[i], ins_thick, pvc) for i in 1:3]

    components = [
        CableComponent("core1", ConductorGroup(sectors[1]), InsulatorGroup(insulators[1])),
        CableComponent("core2", ConductorGroup(sectors[2]), InsulatorGroup(insulators[2])),
        CableComponent("core3", ConductorGroup(sectors[3]), InsulatorGroup(insulators[3]))
    ]


    design = CableDesign("NAYY_O_3x95_30x2_5", components[1])
    add!(design, components[2])
    add!(design, components[3])

    FileIO.save("$(design.cable_id).jld2",design.cable_id, design)


    println("Previewing cable design...")
    plt, _ = preview(design)
    plt 

    design

    cablepos = CablePosition(
        design,
        0.0, -1.0,
        Dict("core1" => 1, "core2" => 2, "core3" => 3)
    )


    system = LineCableSystem("NAYY_O", 1.0, cablepos)
    f = [50.0]
    earth = EarthModel(f, 100.0, 1.0, 1.0) # 100 Ω·m resistivity




    # We create the `LineParametersProblem`.
    problem = LineParametersProblem(system;
        temperature=20.0,
        earth_props=earth,
        frequencies=f
    )


    ## DIFFERENT FORMULATIONS ON SAME CABLE

    # FEM
    compute_fem_formulation!(problem, results_df, earth, materials)

    # ANALYTICAL
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.SimpleCarson(), LineCableModels.Engine.SimpleCarson())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.FullCarson(), LineCableModels.Engine.FullCarson())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.DeriModel(), LineCableModels.Engine.FullCarson())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.DeriModel(), LineCableModels.Engine.DeriModel())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.DeriModel(), LineCableModels.Engine.Saad())

##
####################################################





####################################################
################## NAYCWY 95  ######################
####################################################
##
    n_sectors = 3
    r_back_mm = 10.24       # sector radius b
    d_sector_mm = 9.14        # sector depth s
    r_corner_mm = 1.02      # corner radius c
    theta_cond_deg = 119.0    # sector angle φ
    theta_cond_deg = 120.0    # sector angle φ
    ins_thick = 1.1e-3        # core insulation thickness

    sector_params = SectorParams(
        n_sectors,
        r_back_mm / 1000,
        d_sector_mm / 1000,
        r_corner_mm / 1000,
        theta_cond_deg,
        ins_thick
    )


    rot_angles = (0.0, 120.0, 240.0)
    sectors = [Sector(sector_params, ang, aluminum) for ang in rot_angles]
    insulators = [SectorInsulator(sectors[i], ins_thick, pvc) for i in 1:3]

    components = [
        CableComponent("core1", ConductorGroup(sectors[1]), InsulatorGroup(insulators[1])),
        CableComponent("core2", ConductorGroup(sectors[2]), InsulatorGroup(insulators[2])),
        CableComponent("core3", ConductorGroup(sectors[3]), InsulatorGroup(insulators[3]))
    ]

    # === Neutral Conductor ===
    n_neutral = 30
    d_neutral_strand_mm = 1.58
    r_strand = (d_neutral_strand_mm / 2) * 1e-3
    d_outer_mm = 38.3
    t_neutral_sheath_mm = 2.2

    R_O = (d_outer_mm / 2) * 1e-3
    outer_jacket_thickness = t_neutral_sheath_mm * 1e-3

    # From the geometry, the outer radius of the cable is the sum of the radius to the center of the neutral wires, 
    # the radius of a neutral strand, and the outer jacket thickness.
    # R_O = R_N + r_strand + outer_jacket_thickness
    R_N = R_O - r_strand - outer_jacket_thickness
    inner_radius_neutral = R_N - r_strand

    neutral_wires = WireArray(
        inner_radius_neutral,
        Diameter(2 * r_strand),
        n_neutral,
        0.0, # lay ratio
        copper
    )

    neutral_jacket = Insulator(neutral_wires, Thickness(outer_jacket_thickness), pvc)
    neutral_component = CableComponent("neutral", ConductorGroup(neutral_wires), InsulatorGroup(neutral_jacket))



    design = CableDesign("NAYCWY_O_3x95_30x2_5", components[1])
    add!(design, components[2])
    add!(design, components[3])
    add!(design, neutral_component)

    design_naycwy = design

    FileIO.save("$(design_naycwy.cable_id).jld2",design_naycwy.cable_id, design_naycwy)

    # loaded = FileIO.load("NAYCWYO.jld2","design_naycwy")



    println("Previewing cable design...")
    plt, _ = preview(design)
    plt 

    design

    cablepos = CablePosition(
        design,
        0.0, -1.0,
        Dict("core1" => 1, "core2" => 2, "core3" => 3, "neutral" => 4)
    )


    system = LineCableSystem("NAYCWY_O", 1000.0, cablepos)
    f = [50.0]
    earth = EarthModel(f, 100.0, 1.0, 1.0) # 100 Ω·m resistivity




    # We create the `LineParametersProblem`.
    problem = LineParametersProblem(system;
        temperature=20.0,
        earth_props=earth,
        frequencies=f
    )


    ## DIFFERENT FORMULATIONS ON SAME CABLE

    # FEM
    compute_fem_formulation!(problem, results_df, earth, materials)

    # ANALYTICAL
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.SimpleCarson(), LineCableModels.Engine.SimpleCarson())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.FullCarson(), LineCableModels.Engine.FullCarson())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.DeriModel(), LineCableModels.Engine.FullCarson())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.DeriModel(), LineCableModels.Engine.DeriModel())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.DeriModel(), LineCableModels.Engine.Saad())
##
####################################################





####################################################
################## 4x95 NAYY  ########################
####################################################
##
    materials = MaterialsLibrary(add_defaults=true)
    pvc = Material(1e15, 8.0, 1.0, 20.0, 0.1)             # simple PVC
    add!(materials, "pvc", pvc)
    copper = get(materials, "copper")

    R_dc_nom = 0.32 # ohm/km
    L_m = 1000.0 # m
    A_sector = 92.14e-6 # m²
    rho_from_Rdc = (R_dc_nom / L_m) * A_sector
    aluminum_corrected = Material(rho_from_Rdc, 1.0, 1.0, 20.0, 0.00403) # Aluminum with corrected resistivity from IEC 60287-3-2
    add!(materials, "aluminum_corrected", aluminum_corrected)


    n_sectors = 4
    r_back_mm = 12.00
    d_sector_mm = 10.35
    r_corner_mm = 1.20
    theta_cond_deg = 89.0
    theta_cond_deg = 90.0
    ins_thick = 1.1e-3

    sector_params = SectorParams(
        n_sectors,
        r_back_mm / 1000,
        d_sector_mm / 1000,
        r_corner_mm / 1000,
        theta_cond_deg,
        ins_thick
    )
        
    rot_angles = (0.0, 90.0, 180.0,270.0)
    sectors = [Sector(sector_params, ang, aluminum_corrected) for ang in rot_angles]
    insulators = [SectorInsulator(sectors[i], ins_thick, pvc) for i in 1:4]

    components = [
        CableComponent("core1", ConductorGroup(sectors[1]), InsulatorGroup(insulators[1])),
        CableComponent("core2", ConductorGroup(sectors[2]), InsulatorGroup(insulators[2])),
        CableComponent("core3", ConductorGroup(sectors[3]), InsulatorGroup(insulators[3])),
        CableComponent("core4", ConductorGroup(sectors[4]), InsulatorGroup(insulators[4]))
    ]


    # === Assemble cable design ===
    design = CableDesign("4_core_95NAYY", components[1])
    add!(design, components[2])
    add!(design, components[3])
    add!(design, components[4])



    println("Previewing cable design...")
    plt, _ = preview(design)
    plt 


    cablepos = CablePosition(
        design,
        0.0, -1.0,
        Dict("core1" => 1, "core2" => 2, "core3" => 3, "core4" => 4)
    )

    system = LineCableSystem("4c_NAYY", 1000.0, cablepos)

    f = [50.0]
    earth = EarthModel(f, 100.0, 1e-16, 1.0)
    problem = LineParametersProblem(system;
        temperature=20.0,
        earth_props=earth,
        frequencies=f
    )


    ## DIFFERENT FORMULATIONS ON SAME CABLE

    # FEM
    compute_fem_formulation!(problem, results_df, earth, materials)

    # ANALYTICAL
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.SimpleCarson(), LineCableModels.Engine.SimpleCarson())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.FullCarson(), LineCableModels.Engine.FullCarson())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.DeriModel(), LineCableModels.Engine.FullCarson())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.DeriModel(), LineCableModels.Engine.DeriModel())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.DeriModel(), LineCableModels.Engine.Saad())

##
####################################################


FileIO.save("cables_impedance_data.jld2","results_df",results_df)

results_df











































####################################################
################## Materials  ########################
####################################################
materials = MaterialsLibrary(add_defaults=true)
pvc = Material(Inf, 8.0, 1.0, 20.0, 0.1) # simple PVC
add!(materials, "pvc", pvc)

delete!(materials, "aluminum")
delete!(materials, "copper")
R_dc_nom = 0.32 # ohm/km
L_m = 1000.0 # m
A_sector = 0.88*95 * 1e-6 # m²
rho_from_Rdc = (R_dc_nom / L_m) * A_sector
aluminum = Material(rho_from_Rdc, 1.0, 1.0, 20.0, 0.00403) # Aluminum with corrected resistivity from IEC 60287-3-2
add!(materials, "aluminum", aluminum)
r_strand = 0.79e-3
rho_neutral_from_Rdc = (0.32 / 1000) * (30 * (pi * r_strand^2))
copper = Material(rho_neutral_from_Rdc, 1.0, 1.0, 20.0, 0.00393) # Copper with corrected resistivity from IEC 60287-3-2
add!(materials, "copper", copper)
materials

# FileIO.save("materials.jld2","materials",materials)



# ## PREPARE DATAFRAME TO STORE RESULTS

# results_df = DataFrame(CableID=String[], Formulation=String[], time_taken = Float64[], line_params_Z = Any[], line_params_Y = Any[], Jresiduals=Any[], SE_kpi=Any[])





####################################################
################## NAYY 95  ########################
####################################################
##
    n_sectors = 3
    r_back_mm = 10.24       # sector radius b
    d_sector_mm = 9.14        # sector depth s
    r_corner_mm = 1.02      # corner radius c
    theta_cond_deg = 119.0    # sector angle φ
    theta_cond_deg = 120.0    # sector angle φ
    ins_thick = 1.1e-3        # core insulation thickness

    sector_params = SectorParams(
        n_sectors,
        r_back_mm / 1000,
        d_sector_mm / 1000,
        r_corner_mm / 1000,
        theta_cond_deg,
        ins_thick
    )


    rot_angles = (0.0, 120.0, 240.0)
    sectors = [Sector(sector_params, ang, aluminum) for ang in rot_angles]
    insulators = [SectorInsulator(sectors[i], ins_thick, pvc) for i in 1:3]

    components = [
        CableComponent("core1", ConductorGroup(sectors[1]), InsulatorGroup(insulators[1])),
        CableComponent("core2", ConductorGroup(sectors[2]), InsulatorGroup(insulators[2])),
        CableComponent("core3", ConductorGroup(sectors[3]), InsulatorGroup(insulators[3]))
    ]


    design = CableDesign("A88n_NAYY_O_3x95_30x2_5", components[1])
    add!(design, components[2])
    add!(design, components[3])

    FileIO.save("$(design.cable_id).jld2",design.cable_id, design)


    println("Previewing cable design...")
    plt, _ = preview(design)
    plt 

    design

    cablepos = CablePosition(
        design,
        0.0, -1.0,
        Dict("core1" => 1, "core2" => 2, "core3" => 3)
    )


    system = LineCableSystem("A88n_NAYY_O", 1.0, cablepos)
    f = [50.0]
    earth = EarthModel(f, 100.0, 1.0, 1.0) # 100 Ω·m resistivity




    # We create the `LineParametersProblem`.
    problem = LineParametersProblem(system;
        temperature=20.0,
        earth_props=earth,
        frequencies=f
    )


    ## DIFFERENT FORMULATIONS ON SAME CABLE

    # FEM
    compute_fem_formulation!(problem, results_df, earth, materials)

    # ANALYTICAL
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.SimpleCarson(), LineCableModels.Engine.SimpleCarson())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.FullCarson(), LineCableModels.Engine.FullCarson())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.DeriModel(), LineCableModels.Engine.FullCarson())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.DeriModel(), LineCableModels.Engine.DeriModel())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.DeriModel(), LineCableModels.Engine.Saad())

##
####################################################





####################################################
################## NAYCWY 95  ######################
####################################################
##
    n_sectors = 3
    r_back_mm = 10.24       # sector radius b
    d_sector_mm = 9.14        # sector depth s
    r_corner_mm = 1.02      # corner radius c
    theta_cond_deg = 119.0    # sector angle φ
    theta_cond_deg = 120.0    # sector angle φ
    ins_thick = 1.1e-3        # core insulation thickness

    sector_params = SectorParams(
        n_sectors,
        r_back_mm / 1000,
        d_sector_mm / 1000,
        r_corner_mm / 1000,
        theta_cond_deg,
        ins_thick
    )


    rot_angles = (0.0, 120.0, 240.0)
    sectors = [Sector(sector_params, ang, aluminum) for ang in rot_angles]
    insulators = [SectorInsulator(sectors[i], ins_thick, pvc) for i in 1:3]

    components = [
        CableComponent("core1", ConductorGroup(sectors[1]), InsulatorGroup(insulators[1])),
        CableComponent("core2", ConductorGroup(sectors[2]), InsulatorGroup(insulators[2])),
        CableComponent("core3", ConductorGroup(sectors[3]), InsulatorGroup(insulators[3]))
    ]

    # === Neutral Conductor ===
    n_neutral = 30
    d_neutral_strand_mm = 1.58
    r_strand = (d_neutral_strand_mm / 2) * 1e-3
    d_outer_mm = 38.3
    t_neutral_sheath_mm = 2.2

    R_O = (d_outer_mm / 2) * 1e-3
    outer_jacket_thickness = t_neutral_sheath_mm * 1e-3

    # From the geometry, the outer radius of the cable is the sum of the radius to the center of the neutral wires, 
    # the radius of a neutral strand, and the outer jacket thickness.
    # R_O = R_N + r_strand + outer_jacket_thickness
    R_N = R_O - r_strand - outer_jacket_thickness
    inner_radius_neutral = R_N - r_strand

    neutral_wires = WireArray(
        inner_radius_neutral,
        Diameter(2 * r_strand),
        n_neutral,
        0.0, # lay ratio
        copper
    )

    neutral_jacket = Insulator(neutral_wires, Thickness(outer_jacket_thickness), pvc)
    neutral_component = CableComponent("neutral", ConductorGroup(neutral_wires), InsulatorGroup(neutral_jacket))



    design = CableDesign("A88n_NAYCWY_O_3x95_30x2_5", components[1])
    add!(design, components[2])
    add!(design, components[3])
    add!(design, neutral_component)

    design_naycwy = design

    FileIO.save("$(design_naycwy.cable_id).jld2",design_naycwy.cable_id, design_naycwy)

    # loaded = FileIO.load("NAYCWYO.jld2","design_naycwy")



    println("Previewing cable design...")
    plt, _ = preview(design)
    plt 

    design

    cablepos = CablePosition(
        design,
        0.0, -1.0,
        Dict("core1" => 1, "core2" => 2, "core3" => 3, "neutral" => 4)
    )


    system = LineCableSystem("A88n_NAYCWY_O", 1000.0, cablepos)
    f = [50.0]
    earth = EarthModel(f, 100.0, 1.0, 1.0) # 100 Ω·m resistivity




    # We create the `LineParametersProblem`.
    problem = LineParametersProblem(system;
        temperature=20.0,
        earth_props=earth,
        frequencies=f
    )


    ## DIFFERENT FORMULATIONS ON SAME CABLE

    # FEM
    compute_fem_formulation!(problem, results_df, earth, materials)

    # ANALYTICAL
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.SimpleCarson(), LineCableModels.Engine.SimpleCarson())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.FullCarson(), LineCableModels.Engine.FullCarson())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.DeriModel(), LineCableModels.Engine.FullCarson())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.DeriModel(), LineCableModels.Engine.DeriModel())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.DeriModel(), LineCableModels.Engine.Saad())
##
####################################################





####################################################
################## 4x95 NAYY  ########################
####################################################
##


    n_sectors = 4
    r_back_mm = 12.00
    d_sector_mm = 10.35
    r_corner_mm = 1.20
    theta_cond_deg = 89.0
    theta_cond_deg = 90.0
    ins_thick = 1.1e-3

    sector_params = SectorParams(
        n_sectors,
        r_back_mm / 1000,
        d_sector_mm / 1000,
        r_corner_mm / 1000,
        theta_cond_deg,
        ins_thick
    )
        
    rot_angles = (0.0, 90.0, 180.0,270.0)
    sectors = [Sector(sector_params, ang, aluminum_corrected) for ang in rot_angles]
    insulators = [SectorInsulator(sectors[i], ins_thick, pvc) for i in 1:4]

    components = [
        CableComponent("core1", ConductorGroup(sectors[1]), InsulatorGroup(insulators[1])),
        CableComponent("core2", ConductorGroup(sectors[2]), InsulatorGroup(insulators[2])),
        CableComponent("core3", ConductorGroup(sectors[3]), InsulatorGroup(insulators[3])),
        CableComponent("core4", ConductorGroup(sectors[4]), InsulatorGroup(insulators[4]))
    ]


    # === Assemble cable design ===
    design = CableDesign("A88n_4_core_95NAYY", components[1])
    add!(design, components[2])
    add!(design, components[3])
    add!(design, components[4])



    println("Previewing cable design...")
    plt, _ = preview(design)
    plt 


    cablepos = CablePosition(
        design,
        0.0, -1.0,
        Dict("core1" => 1, "core2" => 2, "core3" => 3, "core4" => 4)
    )

    system = LineCableSystem("A88n_4c_NAYY", 1000.0, cablepos)

    f = [50.0]
    earth = EarthModel(f, 100.0, 1e-16, 1.0)
    problem = LineParametersProblem(system;
        temperature=20.0,
        earth_props=earth,
        frequencies=f
    )


    ## DIFFERENT FORMULATIONS ON SAME CABLE

    # FEM
    compute_fem_formulation!(problem, results_df, earth, materials)

    # ANALYTICAL
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.SimpleCarson(), LineCableModels.Engine.SimpleCarson())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.FullCarson(), LineCableModels.Engine.FullCarson())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.DeriModel(), LineCableModels.Engine.FullCarson())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.DeriModel(), LineCableModels.Engine.DeriModel())
    compute_dss_formulation!(problem, results_df, LineCableModels.Engine.DeriModel(), LineCableModels.Engine.Saad())

##
####################################################




FileIO.save("cables_impedance_data_Wo_An_A88n.jld2","results_df",results_df)

results_df




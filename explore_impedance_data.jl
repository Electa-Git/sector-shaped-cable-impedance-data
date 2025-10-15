# This script shows a simple way to explore the impedance data stored in a .jld2 
# The cable_impedance_data.jld2 file or more will be attached for more cable designs and configurations
using Pkg
Pkg.activate(".") 
using FileIO
using JLD2
using DataFrames
using CSV

impedance_df = load("cables_impedance_data.jld2", "results_df")  # this .jld2 file or more will be attached for more cable designs and configurations
# or 
#impedance_df = load("cables_impedance_data_Wo_An_A88n.jld2", "results_df")  # this .jld2 file or more will be attached for more cable designs and configurations


for row in eachrow(impedance_df)
    println("Cable ID: ", row.CableID)
    println("Formulation: ", row.Formulation)
    println("Computation Time [s]: ", row.time_taken)
    println("Line Parameters Z: [Ω/km]")
    display(row.line_params_Z)
    println("Line Parameters Y: [S/km]")
    display(row.line_params_Y)

    z_df = DataFrame(row.line_params_Z, :auto)
    y_df = DataFrame(row.line_params_Y, :auto)

    CSV.write("Impedances_csv_files\\$(row.CableID)_$(row.Formulation)_Z.csv", z_df, header=false)
    CSV.write("Impedances_csv_files\\$(row.CableID)_$(row.Formulation)_Y.csv", y_df, header=false)
    @info "Saved series impedance CSV for $(row.CableID) - $(row.Formulation) at Impedances_csv_file\\\$(row.CableID)_$(row.Formulation)_Z.csv"
    @info "Saved shunt admittance CSV for $(row.CableID) - $(row.Formulation) at Impedances_csv_file\\\$(row.CableID)_$(row.Formulation)_Y.csv"

    println("--------------------------------------------------")
end
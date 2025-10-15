# Impedance Data for Sector-Shaped LV Cables

This repository contains the data and models associated with the paper: *Impact of LV Cable Impedance Model Fidelity on Distribution System State Estimation*.

It provides a high-quality, open-access dataset of series impedance parameters for common Low Voltage (LV) underground cables with sector-shaped conductors. The goal of this dataset is to support reproducible research in LV network analysis and to provide a benchmark for future studies.

## Repository Content

This repository contains:

**Impedance Matrices `/Impedances_csv_files`:** A comprehensive dataset of calculated series impedance matrices for various LV cable types. The data can be exported to a structured format for easy use in power system simulation tools.
*   **Calculation Code `test_cable_case_studies.jl` :**  This code computes impedance values are provided from multiple modeling approaches discussed in the paper:
    *   **Finite Element (FE) Models:** High-fidelity benchmark results obtained from 2D FE simulations. These models serve as the ground truth for our comparisons.
    *   **Analytical Models:** Results from various analytical formulations, including:
        *   Classical analytical equations (e.g., based on Carson's work).
        *   Analytical models with our proposed corrections for sector-shaped conductor geometry.
        *   Models showing the impact of common simplifications (e.g., neglecting earth return, ignoring sheath shielding effects).

    > [!NOTE]
    > The calculation code requires using a development fork of `LineCableModels.jl`. Please clone the fork from [https://github.com/MohamedNumair/LineCableModels.jl/tree/integration/sector-fem-dss](https://github.com/MohamedNumair/LineCableModels.jl/tree/integration/sector-fem-dss) and checkout the `integration/sector-fem-dss` branch that includes sector-shaped cable support and OpenDSS analytical methods. This requirement will remain until the relevant pull requests are merged into the main `LineCableModels.jl` package.

*   **FE Models `FEM-files`:** The source files for the Finite Element models used to generate the benchmark impedance data. This allows for verification, extension, and further analysis by other researchers.

## Relation to the Paper

The data in this repository forms the basis for the analysis presented in our paper, "Impact of LV Cable Impedance Model Fidelity on Distribution System State Estimation". The paper uses this dataset to:

1.  Validate analytical impedance models with geometric corrections against high-fidelity FE benchmarks.
2.  Conduct a sensitivity analysis on the effects of earth return and cable shielding.
3.  Quantify how different impedance modeling assumptions impact the accuracy of Distribution System State Estimation (DSSE).

We encourage users of this dataset to cite our paper.

## Future Updates

This repository is a living project. We plan to update it with:

*   Impedance data for a wider range of cable types and configurations.
*   Analysis of frequency-dependent effects for harmonic and transient studies.
*   Models for cables with different material properties (e.g., steel-armored cables).
*   Shunt Admittance Data, the shunt admittance matrices (`G` and `B` matrices) derived from the FE models are to be included for completeness.

All code and datasets are provided to ensure the reproducibility of our results and to foster further research into high-fidelity modeling of modern LV distribution networks.

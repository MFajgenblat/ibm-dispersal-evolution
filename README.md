# Code for the individual-based model in Fajgenblat, De Meester & Urban (2024)
This repository contains the code for the individual-based model presented in:

Fajgenblat, M., De Meester, L. & Urban, M.C. (2024). Dispersal evolution alters evolution-mediated priority effects in a metacommunity. Proceedings of the Philosophical Transactions B.

`model.R` contains the code for the individual-based model, defined as a function with several arguments that pertain to model parameters.

`species_pool_scenarios.csv` comprises the 16 species pool scenarios that are evaluated in the paper.

`run_scenarios.R` actually runs the individual-based model for all scenarios (full factorial combination of model parameters) and saves their combined output.

`visualization.R` turns the combined output into the visualizations shown in the paper.

`demo_monopolization.R` features a slightly modified version of the individual-based model to visualize the process of evolution-mediated priority effects.

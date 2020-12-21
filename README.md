# CF-Oxygen

Spatial and non-spatial agent based oxygen Cystic Fibrosis models in Julia. 

Modeling climax (C) and attack (F) communities, the climax community has an oxygen dependent growth rate r_c(w). The non-spatial (spatially homogeneous) model has a global oxygen amount and cells can reproduce anywhere in the grid, oxygen is governed by a mass-balance ODE and is consumed by C. The spatial model has an n x n system of ODE's that determine the oxygen amount in each space. Non-spatial version uses defualt Julia ODE solver, spatial solves the n x n system with RK4.

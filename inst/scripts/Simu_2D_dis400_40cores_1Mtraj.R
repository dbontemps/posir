
## general parameters

dimension = 2
sim_path="Simulations"
maxbatchsize = 4*10**6

## parameters for trajectories simulation

Ndiscretisation = 400
deltagrid = seq(100,1,-1)/100 # decreasing order

## Parameters for simulation quantiles estimation

Ntraj_batch = 500
alphagrid = seq(500,1,-1)/1000

## Parameters of Monte-carlo estimations

Ncores = 40
Ntrajmin = 10^6

## File names parameters

Namedeltagrid = "1to.01by.01"

# library loading

library(posir)
library(future)
logger::log_threshold(logger::DEBUG, namespace = "posir")

# Exec

plan(multisession, workers = Ncores)
compute_quantiles(Ntrajmin, Ntraj_batch, Ndiscretisation, deltagrid,
                  Namedeltagrid, alphagrid, maxbatchsize, sim_path,
                  d=dimension)
plan(sequential)

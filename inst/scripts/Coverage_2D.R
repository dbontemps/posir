
## general parameters

dimension = 2
sim_path="Simulations"
maxbatchsize = 4*10**6

## parameters for trajectories simulation

Ndiscretisation = 400
deltagrid = seq(100,1,-1)/100 # decreasing order
positions = seq(1,91,10) # only delta values for which n*delta is an integer
aux=deltagrid[positions]*n
if(sum(abs(round(aux)-aux)<.001)<length(positions)) {
  print("ERROR : incorect delta grid")
  stop()
}

## simulation parameters

Ntraj_batch = 500

## Parameters for Monte-carlo estimations

Ncores = 40
Ntrajmin = 10^6

## File names parameters

Namedeltagridinit = "1to.01by.01"
Namedeltagridfinal = "1to.1by.1"

# library loading

library(posir)
library(future)
logger::log_threshold(logger::DEBUG, namespace = "posir")

# Exec

plan(multisession, workers = Ncores)
for (n in c(50, 100, 200, 400)) {
  compute_error_levels(Ntrajmin, Ntraj_batch, n, deltagrid, positions,
                       Namedeltagridfinal, maxbatchsize, sim_path,
                       rdistrib = rCenteredPareto, NameDis = "Pareto3",
                       d=dimension)
}
plan(sequential)

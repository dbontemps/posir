
## general parameters

dimension = 1
sim_path = "Simulations"
maxbatchsize = 5*10**6

## parameters for trajectories simulation

Ndiscretisation = 20000
deltagrid = seq(200,1,-1)/200 # decreasing order
Ntraj_batch = 500

## Parameters of Monte-carlo estimations

Ncores = 4
Ntrajmin = 10^6

## File names parameters

Namedeltagridinit = "1to.005by.005"

# fonctions definalitions

library(posir)
library(EnvStats)
library(future)
logger::log_threshold(logger::DEBUG, namespace = "posir")

myfun = function(n, positions, Namedeltagridfinal,
                 rdistrib = rCenteredPareto, NameDis = "Pareto3") {
  aux=deltagrid[positions]*n
  if(sum(abs(round(aux)-aux)<.001)<length(positions)) {
    print("ERREUR : grille pour delta incorrecte")
    stop()
  }
  compute_error_levels(Ntrajmin, Ntraj_batch, n, deltagrid, positions,
                       Namedeltagridfinal, maxbatchsize, sim_path,
                       rdistrib = rdistrib, NameDis = NameDis, d=dimension)
}

# Exec

plan(multisession, workers = Ncores)

# several discretizations with Pareto 3 distribution
for(j in c(30, 50, 100, 400, 1000, 5000)) {
  myfun(n = j, positions = seq(1,181,20),
        Namedeltagridfinal = "1to.1by.1",
        rdistrib = rCenteredPareto, NameDis = "Pareto3")
}
plan(sequential)

# for(j in c(30, 50, 100, 400, 1000, 5000, 10000, 20000)) {
for(j in c(30, 50, 100, 400, 1000, 5000)) {
    myfun(n = j, positions = seq(1,181,20),
        Namedeltagridfinal = "1to.1by.1",
        rdistrib = function(k) {
          return(EnvStats::rpareto(k, location=1, shape=2.1) - (2.1/(2.1-1)))
        },
        NameDis = "Pareto2.1")
}

# Idem with Laplace distribution
rlaplace = function(j) {
  return(rexp(j)*(2*rbiName(j,size=1,prob=.5)-1))
}
for(j in c(30, 50, 100, 400)) {
  myfun(n = j, positions = seq(1,181,20),
        Namedeltagridfinal = "1to.1by.1",
        rdistrib = rlaplace, NameDis = "Laplace")
}

# Idem with normal distribution
for(j in c(30, 50, 100, 400)) {
  myfun(n = j, positions = seq(1,181,20),
        Namedeltagridfinal = "1to.1by.1",
        rdistrib = rnorm, NameDis = "Gauss")
}


# Pareto symétrisé avec discrétisation 100
# xm=1
# shapepareto = 2.1
# NameLoicourt = "Pareto2.1Sym"
# LoiErreurs = function(j) {
#   return(EnvStats::rpareto(j, location=xm, shape=shapepareto)
#          -(shapepareto*xm/(shapepareto-1))*(2*rbiName(j,size=1,prob=.5)-1)
#   )
# }
# myfun(n = 100, positions = seq(1,181,20),
#       Namedeltagridfinal = "1to.1by.1",
#       rdistrib = LoiErreurs, NameDis = NameLoicourt)

# Several Pareto shapes with discretization 100
for(shapepareto in c(2.4, 4, 10)) { #shapepareto=2.1 et 3 already done
  myfun(n = 100, positions = seq(1,181,20),
        Namedeltagridfinal = "1to.1by.1",
        rdistrib = function(j) {
          return(EnvStats::rpareto(j, location=1, shape=shapepareto)
                 -(shapepareto/(shapepareto-1)) # *(2*rbiName(j,size=1,prob=.5)-1)
          )
        }, NameDis = paste("Pareto",toString(shapepareto),sep=""))
}

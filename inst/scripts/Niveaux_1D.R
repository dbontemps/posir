
# Paramètres


## Paramètres généraux

dimension = 1
sim_path="Simulations"
Taillemaxbatch = 5*10**6

## Paramètres de simulations des trajectoires pour les quantiles

Ndiscretisation = 20000 # pour les quantiles
grilledelta = seq(200,1,-1)/200 # en ordre décroissant, valeur 1 permise
Ntraj_batch = 500
grillealpha = seq(500,1,-1)/1000

## Paramètres des estimations Monte-carlo

Ncores = 40
Ntrajmin = 10^6

## Paramètres de nommage des fichiers

NomfichierQuantiles = "Table_quantiles_1D.txt"
Nomgrilledeltainit = "1to.005by.005"

# Chargement des fonctions

library(posir)
library(EnvStats)
library(future)
logger::log_threshold(logger::DEBUG, namespace = "posir")

myfun = function(n, positions, Nomgrilledeltafin,
                 rdistrib = rCenteredPareto, NameDis = "Pareto3") {
  aux=grilledelta[positions]*n
  if(sum(abs(round(aux)-aux)<.001)<length(positions)) {
    print("ERREUR : grille pour delta incorrecte")
    stop()
  }
  compute_error_levels(Ntrajmin, Ntraj_batch, Ndiscretisation, n, grilledelta,
                       positions, Nomgrilledeltafin, Nomgrilledeltainit,
                       Taillemaxbatch, sim_path, NameFQ = NomfichierQuantiles,
                       rdistrib = rdistrib, NameDis = NameDis,
                       d=dimension)
}

# Exécutions

plan(multisession, workers = Ncores)

# différentes discrétisations pour Pareto3
for(j in c(30, 50, 100, 400, 1000, 5000)) {
  myfun(n = j, positions = seq(1,181,20),
        Nomgrilledeltafin = "1to.1by.1",
        rdistrib = rCenteredPareto, NameDis = "Pareto3")
}
plan(sequential)

for(j in c(30, 50, 100, 400, 1000, 5000, 10000, 20000)) {
  myfun(n = j, positions = seq(1,181,20),
        Nomgrilledeltafin = "1to.1by.1",
        rdistrib = function(k) {
          return(EnvStats::rpareto(k, location=1, shape=2.1) - (2.1/(2.1-1)))
        },
        NameDis = "Pareto2.1")
}

# Idem avec une distribution de Laplace
rlaplace = function(j) {
  return(rexp(j)*(2*rbinom(j,size=1,prob=.5)-1))
}
for(j in c(30, 50, 100, 400)) {
  myfun(n = j, positions = seq(1,181,20),
        Nomgrilledeltafin = "1to.1by.1",
        rdistrib = rlaplace, NameDis = "Laplace")
}

# Pareto symétrisé avec discrétisation 100
# xm=1
# shapepareto = 2.1
# NomLoicourt = "Pareto2.1Sym"
# LoiErreurs = function(j) {
#   return(EnvStats::rpareto(j, location=xm, shape=shapepareto)
#          -(shapepareto*xm/(shapepareto-1))*(2*rbinom(j,size=1,prob=.5)-1)
#   )
# }
# myfun(n = 100, positions = seq(1,181,20),
#       Nomgrilledeltafin = "1to.1by.1",
#       rdistrib = LoiErreurs, NameDis = NomLoicourt)

# Pareto recentré avec discrétisation 100, différentes shapes
for(shapepareto in c(2.4, 4, 10)) { #shapepareto=2.1 et 3 déjà faits ailleurs
  myfun(n = 100, positions = seq(1,181,20),
        Nomgrilledeltafin = "1to.1by.1",
        rdistrib = function(j) {
          return(EnvStats::rpareto(j, location=1, shape=shapepareto)
                 -(shapepareto/(shapepareto-1)) # *(2*rbinom(j,size=1,prob=.5)-1)
          )
        }, NameDis = paste("Pareto",toString(shapepareto),sep=""))
}

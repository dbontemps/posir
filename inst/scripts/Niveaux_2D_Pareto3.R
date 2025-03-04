
# Paramètres


## Paramètres généraux

dimension = 2
sim_path="Simulations"
Taillemaxbatch = 4*10**6 # max 4*10^6 sur mon PC portable ? 10^8 chez moi (avec 24Go de RAM)

## Paramètres de simulations des trajectoires

Ndiscretisation = 400 # pour les quantiles
n = 50 #, 100, 200, 400 # pour les niveaux effectifs
grilledelta = seq(100,1,-1)/100 # en ordre décroissant, valeur 1 permise
positions = seq(1,91,10) # donc delta décroit de 1 à .1 par .1
aux=grilledelta[positions]*n
if(sum(abs(round(aux)-aux)<.001)<length(positions)) {
  print("ERREUR : grille pour delta incorrecte")
  stop()
}

## Paramètres de simulation

Ntraj_batch = 500
grillealpha = seq(500,1,-1)/1000

## Paramètres des estimations Monte-carlo

Ncores = 40
Ntrajmin = 10^6

## Paramètres de nommage des fichiers

#NomfichierQuantiles = "Table_quantiles_2D.txt"
Nomgrilledeltainit = "1to.01by.01"
Nomgrilledeltafin = "1to.1by.1"

# Chargement des fonctions

library(posir)
library(future)
logger::log_threshold(logger::DEBUG, namespace = "posir")

# Exécutions

plan(multisession, workers = Ncores)
for (n in c(50, 100, 200, 400)) {
  compute_error_levels(Ntrajmin, Ntraj_batch, n, grilledelta, positions,
                       Nomgrilledeltafin, Taillemaxbatch, sim_path,
                       rdistrib = rCenteredPareto, NameDis = "Pareto3",
                       d=dimension)
}
plan(sequential)

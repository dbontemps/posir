
# Paramètres

## Paramètres généraux

dimension = 1
sim_path="Simulations"
Taillemaxbatch = 10**7 # max 10^8 chez moi (avec 24Go de RAM)

## Paramètres de simulations des trajectoires

Ndiscretisation = 20000 # max 2000 sur mon PC portable en 2D, 10^4 sur mon PC fixe perso
grilledelta = seq(200,1,-1)/200 # en ordre décroissant, valeur 1 permise

## Paramètres de simulation des quantiles

Ntraj_batch = 500
grillealpha = seq(500,1,-1)/1000

## Paramètres des estimations Monte-carlo

Ncores = 40
Ntrajmin = 10^6

## Paramètres de nommage des fichiers

Nomgrilledelta = "1to.005by.005"
#BaseNomQuantile = "Table_quantiles_1Mtraj.txt"

## Fichier quantile utilisé pour le calcul des niveaux effectifs, si différent

# Chargement des fonctions

library(posir)
library(future)
logger::log_threshold(logger::DEBUG, namespace = "posir")

# Exécutions

plan(multisession, workers = Ncores)
compute_quantiles(Ntrajmin, Ntraj_batch, Ndiscretisation, grilledelta,
                  Nomgrilledelta, grillealpha, Taillemaxbatch, sim_path,
                  d=dimension)
plan(sequential)

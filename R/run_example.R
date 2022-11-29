graphics.off()

#################################################
#################################################
### Definitive example script
#################################################
#################################################

set.seed(17)

### Imports
library("usethis")
library("roxygen2")
library("devtools")
library(Rcpp)
library(RcppEigen)
library(methods)
library(lattice)
library(ggplot2)
library(dplyr)
library(purrr)
library(Matrix)
library("viridis")
library("optAM")
library(hrbrthemes)
library(grid)
library(raster)
source("R/building_example_functions.R")
source("R/plot_exemple_functions.R")
### Global Details
col = c('green3','green2','blue1','blue3')
col = terrain.colors(20)
hea = (viridis(10))

#################################################
## Initialisation & Variables
#################################################

# Map - size
nrow = 50
ncol = 50
nr = seq(0,1,length=nrow)
nc = seq(0,1,length=ncol)

# Map - time caracteristics
N_cycles = 50
Tadd_cyc = 1.5

# Map - space caracteristics
Tmin_alt_0 = 10
Tmax_alt_0 = 15
Tadd_alt_1 = -8
Trange = Tmax_alt_0 - (Tmin_alt_0 + Tadd_alt_1)
Tmarge = c(Tmin_alt_0 + Tadd_alt_1, Tmax_alt_0+Tadd_cyc)

# Algorithm - genetic algo caracteristics
npop = 50
nsur = 10
ngen = 10

# Algorithm - optimum condition values
threshold = 50
confidence = 0.8

# Species - caracteristics
Trange_spe = c(8.8,9.1,11.6,11.9)
migr_spe = array(0.012, c(3, 3))
migr_spe[2,2] = 0.99

# Species - presence - area 1
nb_cell_occuped1 = 10
xmin = 0.15
xmax = 0.32
ymin = 0.2
ymax = 0.32
lim1 = c(xmin,xmax,ymin,ymax)

# Species - presence - area 2
nb_cell_occuped2 = 10
xmin = 0.75
xmax = 0.9
ymin = 0.75
ymax = 0.9
lim2 = c(xmin,xmax,ymin,ymax)



#################################################
## Maps construction
#################################################

# Common for all the figures
height_map = height_map(nrow,ncol)
climate_maps = climat_maps_maker(height_map,Tadd_cyc,N_cycles)
suitability_maps = suitability_maps(climate_maps,Trange_spe)
presence_1st_area = presence_map(nrow,ncol,list_suit,height_map,lim1,nb_cell_occuped1)
presence_2nd_area = presence_map(nrow,ncol,list_suit,height_map,lim2,nb_cell_occuped2)
presence_map = presence_1st_area + presence_2nd_area

presence_map[is.na(presence_map)] = 0
presence_map = Matrix(presence_map, sparse = T)
presence_map = as(presence_map,"dgCMatrix")

cost1 = cost_map(nrow,ncol,height_map,cost_mapping1)
cost2 = cost_map(nrow,ncol,height_map,cost_mapping2)
cost3 = cost_map(nrow,ncol,height_map,cost_mapping3)

#################################################
## Maps plot
#################################################

do_plot_fig2(nr, nc, height_map, climate_maps, N_cycles)
do_plot_fig3(nr, nc, cost1,cost2, cost3)

#################################################
## Results
#################################################

gss = rcpp_global_suitable_sites(suitability_maps)
gsc = rcpp_global_suitable_coordinates(gss)
ltm = rcpp_local_transition_matrix(gss,gsc,migr_spe)
tm = rcpp_transition_matrices(list_suit,ltm,gsc)
cm = rcpp_colonisation_matrices(tm)

vs = rcpp_viable_sites(cm)
vt = rcpp_viable_triplets(vs,cm,gsc,gss,cost)
vv = rcpp_viable_values(vt,vs,gss,cm)
#cv = rcpp_get_current_vector(pres,cm,gss)

ph = rcpp_pheromons(vt)
ecp = rcpp_eval_current_prob(threshold,presence_map[nrow:1,],cm,gss)
ntp = threshold - which(cumsum(ecp)>0.95)[1] 

po = rcpp_generate_population(ph,gss,npop,ntp)
resultat1 = rcpp_algorithm_opt(ph,vt,po,cost1,presence_map[nrow:1,],cm,gss,vv,threshold,confidence,npop,nsur,ngen,ntp)
resultat2 = rcpp_algorithm_opt(ph,vt,po,cost2,presence_map[nrow:1,],cm,gss,vv,threshold,confidence,npop,nsur,ngen,ntp)
resultat3 = rcpp_algorithm_opt(ph,vt,po,cost3,presence_map[nrow:1,],cm,gss,vv,threshold,confidence,npop,nsur,ngen,ntp)
choix1 = rcpp_result_to_choice(resultat1,vt)
choix2 = rcpp_result_to_choice(resultat2,vt)
choix3 = rcpp_result_to_choice(resultat3,vt)


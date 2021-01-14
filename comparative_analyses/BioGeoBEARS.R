library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))

wd = np("~/Dropbox/Encelia/analysis/phylogeny/BioGeoBEARS/")
setwd(wd)

t = read.tree("~/Dropbox/Encelia/analysis/phylogeny/astral/astral.miss0.6.brlen.dated.tre")
td_nom = drop.tip(td, c("Encelia_californica1", "Encelia_virginensis2",
                        "Encelia_frutescens_glandulosa", "Encelia_farinosa_phenicodonta",
                        "Xylorhiza_tortifolia", "Enceliopsis_covillei"))
td_nom$tip.label = gsub("\\d", "", td_nom$tip.label)
write.tree(td_nom, "nominal.tre")

trfn = "nominal.tre"
geogfn = "regions.txt"
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)

numstates_from_numareas(numareas=6, maxareas=4, include_null_range=FALSE)
numstates_from_numareas(numareas=6, maxareas=3, include_null_range=TRUE)
numstates_from_numareas(numareas=6, maxareas=2, include_null_range=TRUE)
max_range_size = 4

#######################################################
# Run DEC
#######################################################

# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

BioGeoBEARS_run_object$on_NaN_error = -1e50    
# returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE
# shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"   
# if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 1
# to avoid pathology
BioGeoBEARS_run_object$force_sparse = FALSE 

BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

runslow = FALSE
resfn = "encelia_DEC.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}

#######################################################
# Run DEC + J
#######################################################

# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "encelia_DEC_J.Rdata"
runslow = FALSE
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resDECj = res
} else {
  # Loads to "res"
  load(resfn)
  resDECj = res
}



# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
class(td_nom) = "phylo"
res2 = plot_BioGeoBEARS_results(results_object,
                                plotwhat="text", label.offset=0.45,
                                tipcex=0.7, statecex=0.7, 
                                splitcex=0.6, titlecex=0.8, 
                                plotsplits=TRUE, 
                                cornercoords_loc=scriptdir, 
                                include_null_range=TRUE, tr=td_nom, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object,plotwhat="pie", 
                         label.offset=0.45, tipcex=0.7, 
                         statecex=0.7, splitcex=0.6, titlecex=0.8, 
                         plotsplits=TRUE, cornercoords_loc=scriptdir, 
                         include_null_range=TRUE, tr= td_nom, 
                         tipranges=tipranges)


root = results_object$ML_marginal_prob_each_state_at_branch_top_AT_node[13, ]
names(root) = seq(1, length(root))
sum(root[root > 0.05])
states = root[root > 0.05]
results_object$inputs$all_geog_states_list_usually_inferred_from_areas_maxareas[as.numeric(names(states))]


############################
#  statistical testing
############################

LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)

numparams2 = 2
numparams1 = 3
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

# DEC, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# DEC+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

rbind(res2, res1)
tmp_tests = conditional_format_table(stats)
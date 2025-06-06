# Numerical experiments of cknockoff

## Description and structure

This repo hosts the code used to generate figures in the cKnockoff paper.

`R/`: host R code for the numerical experiments

`data/`: store the numerical results

`figs/`: `R/plot.R` will read the numerical results in `data/` and save the figures in `figs/`

`jobs/`: Slurm job files to be submitted to cluster for large scale simulation

`log/`: store the Slurm logs.

### Dependency

Required packages:
- dbh: `devtools::install_github("lihualei71/dbh")`
- knockoff
- cknockoff: `devtools::install_github("yixiangluo/cknockoff")`
- glmnet
- KernSmooth
- here
- abind
- foreach
- doParallel

## Usage

### FDR and Power

#### local experiments

1. Create a `my_expr_setting.R` file in `R/settings`, in which specify the experiment settings. Existing setting files in `R/settings` can be used to reproduce the numerical experiment figures in the paper.
2. Run `Rscript R/do_expr.R my_expr_setting`.
3. The real time progress are shown in `data/temp/progress-my_expr_setting.txt`, the results are stored in `data/my_expr_setting.RData`, and the produced figure is `figs/simu-my_expr_setting.pdf`.

#### cluster experiments

1. Echo the names of the experiments (e.g. `my_expr_setting`) in `jobs/expr_names.sh`
2. Specify the Slurm arguments in `jobs/expr_template.sh` according your cluster settings. e.g. change the notification email at line 11. But please do NOT add/delete lines or change the order of lines, since this is a template file the content will be read and modified by `R/cluster/set_expr.R` according to which line the content is in.
3. Run `jobs/set_exprs.sh` (e.g. `sbatch jobs/set_exprs.sh`). This will prepare the shared data and structure used by each task in cluster computing.
4. Submit the jobs to the cluster by `sh jobs/do_exprs.sh`.
5. Check the progress at `sh R/cluster/show_progress.R`.
6. After all the experiments are done, run `sh jobs/post_process.sh`. This will aggregate the numerical results from all tasks and draw figures.
7. If there is an error, check the log in `log/`.

### Scalability

Modify the experiment settings in `R/illustrative/scale_m_expr.R` or `R/illustrative/scale_n_expr.R` and run them.

### HIV drug resistance

Run the files in `R/HIV_drug_expr` sequentially: `HIV_preprocess.R`, `HIV_expr.R`, `postprocess_HIV.R`, and `plot_HIV_expr.R`.

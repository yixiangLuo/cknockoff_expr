# Numerical experiments of cknockoff

## Usage

### FDR and Power

Create a `expr_setting.R` file in `R/settings` (see existed examples). Run `Rscript R/do_expr.R expr_settings.R`. The real time progress are shown in `data/temp/progress-expr_setting.txt`, the results are stored in `data/expr_setting.RData`, and the produced figure is `figs/simu-expr_setting.pdf`.

### Scalability

Modify and run `R/scale_expr.R`. Figure in `/figs`.

### HIV drug resistance

Run the R files in `R/HIV_drug_expr`. Figure in `/figs`.

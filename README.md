# EczemaTreat

This repository contains the code developed for the paper "Data-driven personalised recommendations for eczema treatment using a Bayesian model of severity dynamics" (submitted to publication, [preprint here](https://doi.org/10.1101/2024.01.21.24301575)).

The code is written in the R language for statistical computing and the probabilistic programming language [Stan](https://mc-stan.org/) for the models.

## File structure

This repository is organized as a [research compendium](https://doi.org/10.1080/00031305.2017.1375986), with a similar structure as R packages.
Nevertheless, this project is not a literal package with a DESCRIPTION file but prefers to use [renv](https://rstudio.github.io/renv/index.html) to manage package dependencies (see [Reproduciblity section](#Reproduciblity)) and git tags for version control.

Functions specific to this project are located in the [`R/`](R/) directory.

### Stan models

All Stan models are implemented using a single Stan file, [FullModel.stan](models/FullModel.stan), with optional parameters that can be switched on and off for evaluating the contribution of the different model components.
The model is manipulated using `ScoradPred` objects (inherinting from the class `EczemaModel` defined in [EczemaPred](https://ghurault.github.io/EczemaPred/)).

### Analysis scripts

The analysis code is located in the [`analysis/`](analysis/) directory:

- [`01_check_models.R`](analysis/01_check_models.R): conduct prior predictive checks and fake data checks.
- [`02_run_fit.R`](analysis/02_run_fit.R): fit the model to data.
- [`03_check_fit.Rmd`](analysis/03_check_fit.Rmd): diagnose fit, inspect posterior, posterior predictive checks.
- [`04a_run_validation.R`](analysis/04a_run_validation.R): run the validation process (forward chaining).
- [`04b_run_validation_reference.R`](analysis/04b_run_validation_reference.R): run the validation process (forward chaining) for the reference (univariate) models.
- [`05_check_performance.Rmd`](analysis/05_check_performance.Rmd): analyse validation results, performance.
- [`06_analyse_recommendations.Rmd`](analysis/06_analyse_recommendations.Rmd): generate and analyse treatment recommendations.

- Scripts to generate figures for the paper:
  - [`07_plot_data.R`](analysis/07_plot_data.R)
  - [`07_plot_fit.R`](analysis/07_plot_fit.R)
  - [`07_plot_performance.R`](analysis/07_plot_performance.R)
  - [`07_plot_powerprior.R`](analysis/07_plot_powerprior.R)

In addition, [`generate_reports.R`](analysis/generate_reports.R) renders reports from [`03_check_fit.Rmd`](analysis/03_check_fit.Rmd) and [`05_check_performance.Rmd`](analysis/05_check_performance.Rmd) for all models and severity items/scores.
[`view_reports.R`](analysis/view_reports.R) creates an HTML document to easily browse these reports.

### Note on the terminology

- "ScoradPred" refers to the base model (independent state-space models for all severity items, defined by an ordered logistic measurement distribution and latent random walk dynamic).
- Modifications/Improvements of the "ScoradPred" model as referred to as "ScoradPred+improvement".
For example, the model consisting of ordered logistic measurement distribution for all severity items and a latent multivariate random walk (i.e. modelling the correlations between changes of latent severity) is denoted as "ScoradPred+corr".
- The different flavours of the base ScoradPred model are implemented in a model class also called ScoradPred.
- The Stan file implementing all of these models is called "FullModel.stan".

## Reproducibility

This project uses [renv](https://rstudio.github.io/renv/index.html) to manage R package dependencies.
The details of the packages needed to reproduce the analysis is stored in [`renv.lock`](renv.lock) and configuration files and the project library (ignored by git) is stored in [`renv/`](renv/).
After installing `renv` itself (`install.packages("renv")`), the project library can be restored by calling `renv::restore(exclude = "TanakaData")`.
Note that this command explicitly avoid installing the `TanakaData` package, a proprietary (unavailable) package containing the data used in this project.
The data is not available at the time of writing.

In addition, we provide a [Dockerfile](Dockerfile) to fully reproduce the computational environment with [Docker](https://www.docker.com/):

- building the image: `docker build . -t eczematreat`.
In addition to installing R packages using renv, the Docker image will also install the correct version of R and system dependencies required to use Stan.
- running the container `docker run -d --rm -p 1212:8787 -e ROOT=TRUE -e DISABLE_AUTH=true -v ${PWD}:/home/rstudio/EczemaTreat -v /home/rstudio/EczemaTreat/renv eczematreat`.
This commands launch an RStudio server session (without authentication, giving the user access to root) accessible at [`http://localhost:1212/`](http://localhost:1212/), while mounting the current directory into the container.

After that, to reproduce the analysis, we suggest to open the RStudio project (`.Rproj` file) and runs the [analysis scripts](analysis/) in the order indicated by their prefix.
Intermediate and output files are saved to a `results/` directory.

NB: this project relies on [EczemaPred version v0.3.0](https://github.com/ghurault/EczemaPred/releases/tag/v0.3.0).

## License

This open source version of this project is licensed under the GPLv3 license, which can be seen in the [LICENSE](LICENSE.md) file.

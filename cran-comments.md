# Second submission

## Responses to CRAN reviewer comments

**1. Examples wrapped in `\dontrun{}`**

> Please unwrap the examples if they are executable in < 5 sec, or create
> additionally small toy examples to allow automatic testing. You could also
> replace `\dontrun{}` with `\donttest{}` if it takes longer than 5 sec.

The `BRRAT()` function calls `rstan::sampling()` (MCMC), which cannot
complete in under 5 seconds. All `\dontrun{}` blocks have been replaced with
`\donttest{}`. The package also includes a `tests/testthat/` suite covering
input validation and basic model fitting that runs under `R CMD check`.

**2. Writing to the user's home filespace**

> Please ensure that your functions do not write by default or in your
> examples/vignettes/tests in the user's home filespace. Please omit any
> default path in writing functions. In your examples/vignettes/tests you can
> write to `tempdir()`.
> -> R/table_helpers.R

The default path arguments (`file_global = 'global.csv'`, `file_per_id =
'per_id.csv'`) have been removed from `export_BRRAT_csv()`. Users must now
supply explicit paths. The example has been updated to use `tempfile()`.

**3. Setting a specific seed within a function**

> Please do not set a seed to a specific number within a function.
> -> R/plotting_helpers.R

Two `set.seed(1)` calls used for downsampling large scatter plots have been
removed from `plot_BRRAT_population()` and `plot_BRRAT_by_id()`.

# Test environments

* Local Windows 11, R version 4.4.3 x86_64-w64-mingw32 (64-bit)
* github actions MacOS, windows, ubuntu
* win-builder (devel and release)

## Github actions

Status: OK

## Win-builder release (devel returned 550 error)

Status: 1 NOTE

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Matt Gibbs <matt.gibbs@csiro.au>'

## Local Windows 11

0 errors | 0 warnings | 0 notes

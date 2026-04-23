# BRRAT 0.0.2

* Replaced `\dontrun{}` with `\donttest{}` in all examples, as the Stan-based
  fitting cannot complete in under 5 seconds.
* Fixed typo in `summarise_BRRAT()` example (`BRRAAT` -> `BRRAT`).
* Removed default file paths from `export_BRRAT_csv()` to avoid writing to
  the user's working directory. Users must now supply explicit paths.
* Updated `export_BRRAT_csv()` example to write to `tempfile()`.
* Removed `set.seed()` calls from `plot_BRRAT_population()` and
  `plot_BRRAT_by_id()`.

# BRRAT 0.0.1

* Initial CRAN submission.

## Submission

This is a new submission.

## Test environments

* local macOS (aarch64-apple-darwin20), R 4.5.1: 0 errors | 0 warnings | 0 notes
  (apart from "New submission" with incoming checks, and a local-only HTML Tidy
  note: "Skipping checking HTML validation: 'tidy' doesn't look like recent
  enough HTML Tidy.")
* win-builder, R-devel (2026-11-06 r90845 ucrt): 1 NOTE (see below)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission (maintainer: Samuele Aglieco).
* "Possibly misspelled words in DESCRIPTION": CWM (acronym for
  cluster-weighted model) and the surnames of cited authors
  (Ingrassia, Minotti, Vittadini, Rigby, Stasinopoulos). All are
  intentional.

## Data provenance

* `ENRON` is derived from the 'Rsafd' package (MIT licence); see `inst/COPYRIGHTS`.
* `ELECTION` was provided by Eric Sanders (Radboud University), with TweetGenie age estimates provided by Dong Nguyen (Utrecht University); both gave permission (personal communication) to redistribute it under CC BY 4.0; see `inst/COPYRIGHTS`.

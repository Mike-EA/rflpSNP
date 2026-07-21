# This silences a well-known R CMD check false positive: R CMD check's
# static analysis cannot always tell that a bare symbol used as a function's
# default argument value (here, `enzymes = restriction_enzymes` in
# scan_pira_candidates()) refers to valid lazy-loaded package data declared
# via LazyData. Without this, `checking R code for possible problems`
# reports a spurious "no visible binding for global variable" NOTE.
utils::globalVariables("restriction_enzymes")

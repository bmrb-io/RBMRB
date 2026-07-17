#' @keywords internal
"_PACKAGE"

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "RBMRB v", utils::packageVersion(pkgname),
    " -- R interface to BioMagResBank (https://bmrb.io)\n",
    "Run bmrb_test() to verify API connectivity.")
}

utils::globalVariables(".")

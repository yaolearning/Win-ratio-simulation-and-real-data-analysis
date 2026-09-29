required_packages <- c("Rcpp", "survival", "devtools")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install these packages first: ",
    paste(missing_packages, collapse = ", ")
  )
}

if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the MaxWin project root.")
}

desc <- read.dcf("DESCRIPTION")
if (!identical(unname(desc[1, "Package"]), "MaxWin")) {
  stop("DESCRIPTION Package field must be 'MaxWin'.")
}

message("1/2 Generating Rcpp interfaces...")
Rcpp::compileAttributes(".")

message("2/2 Loading MaxWin...")
devtools::load_all(".", reset = TRUE)

message("")
message("MaxWin loaded successfully.")
message("Generated files should now include:")
message("  R/RcppExports.R")
message("  src/RcppExports.cpp")

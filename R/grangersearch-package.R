#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom generics tidy
#' @importFrom generics glance
#' @importFrom rlang enquo eval_tidy as_label
#' @importFrom stats complete.cases
#' @importFrom tibble tibble
## usethis namespace: end
NULL

# Re-export tidy and glance generics
#' @export
generics::tidy

#' @export
generics::glance

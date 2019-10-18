#' Ownership matrix in BLP's car example.
#'
#' @format Dummy variables.
#' \describe{
#'   \item{column i}{1, if product in row j is produced by firm i, 0 otherwise}
#' }
#' @source \url{https://dataverse.harvard.edu/file.xhtml?persistentId=doi:10.7910/DVN/26803/SOF9FW&version=1.0}
"dummies_cars"


#' Product data of BLP's car example.
#'
#' @format A data frame with product data of 2217 cars in 20 markets.
#' \describe{
#'   \item{share}{car market share,}
#'   \item{price}{car price,}
#'   \item{hpwt}{horsepower-weight ratio,}
#'   \item{air}{1, if car has air conditioning, 0 otherwise,}
#'   \item{mpg}{market identifier,}
#'   \item{space}{length times width of the car,}
#'   \item{const}{constant,}
#'   \item{id}{uniquely identifies a car,}
#'   \item{cdid}{uniquely identifies the market of a product,}
#'   \item{firmid}{uniquely identifies the firm of a product (corresponds to column number in the ownership matrix).}
#'    }
#' @source \url{https://dataverse.harvard.edu/file.xhtml?persistentId=doi:10.7910/DVN/26803/SOF9FW&version=1.0}
"productData_cars"

#' Draws for observed heterogeneity in Nevo's cereal example.
#'
#' @format Draws for observed heterogeneity for each demographic.
#' \describe{
#'   \item{cdid}{market identifier,}
#'   \item{draws_}{20 draws differing across markets.}
#' }
#' @source \url{https://dataverse.harvard.edu/file.xhtml?persistentId=doi:10.7910/DVN/26803/SOF9FW&version=1.0}
"demographicData_cereal"


#' Draws for unobserved heterogeneity in Nevo's cereal example.
#'
#' @format Each list entry contains draws (unobserved heterogeneity) for a random coefficient.
#' \describe{
#'   \item{cdid}{market identifier,}
#'   \item{draws_}{20 draws differing across markets.}
#' }
#' @source \url{https://dataverse.harvard.edu/file.xhtml?persistentId=doi:10.7910/DVN/26803/SOF9FW&version=1.0}
"originalDraws_cereal"


#' Product data of Nevo's cereal example.
#'
#' @format A data frame with product data of 24 cereals in each of 94 markets.
#' \describe{
#'   \item{share}{cereals market share,}
#'   \item{price}{cereals price,}
#'   \item{const}{constant,}
#'   \item{sugar}{cereals sugar,}
#'   \item{mushy}{cereals mushy,}
#'   \item{cdid}{market identifier,}
#'   \item{product_id}{uniquely identifies a product in a market,}
#'   \item{productdummy}{uniquely identifies a product in a market,}
#'   \item{IV1}{1. instrument,}
#'   \item{IV2}{2. instrument,}
#'  \item{IV3}{3. instrument,}
#'  \item{IV4}{4. instrument,}
#'  \item{IV5}{5. instrument,}
#'  \item{IV6}{6. instrument,}
#'  \item{IV7}{7. instrument,}
#'  \item{IV8}{8. instrument,}
#'  \item{IV9}{9. instrument,}
#'  \item{IV10}{10. instrument,}
#'  \item{IV11}{11. instrument,}
#'  \item{IV12}{12. instrument,}
#'  \item{IV13}{13. instrument,}
#'  \item{IV14}{14. instrument,}
#'  \item{IV15}{15. instrument,}
#'  \item{IV16}{16. instrument,}
#'  \item{IV17}{17. instrument,}
#'  \item{IV18}{18. instrument,}
#'  \item{IV19}{19. instrument,}
#'  \item{IV20}{20. instrument}
#'    }
#' @source \url{https://dataverse.harvard.edu/file.xhtml?persistentId=doi:10.7910/DVN/26803/SOF9FW&version=1.0}
"productData_cereal"



#' Parameter starting guesses for Nevo's cereal example.
#'
#' @format A matrix with 4 random coefficients (rows) and columns for 4 demographics and one unobserved heterogeneity column (5 cols in total).
#' @source \url{https://dataverse.harvard.edu/file.xhtml?persistentId=doi:10.7910/DVN/26803/SOF9FW&version=1.0}
"theta_guesses_cereal"


#' Mean utility starting guesses for Nevo's cereal example.
#'
#' @format A numeric vector of 2256 values.
#' @source \url{https://dataverse.harvard.edu/file.xhtml?persistentId=doi:10.7910/DVN/26803/SOF9FW&version=1.0}
"w_guesses_cereal"

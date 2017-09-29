#' Expression data set of grafting.
#'
#' A dataset containing the length normed count data from intact, grafted
#' and separated samples. 
#' 
#' @format A data frame with 27416 rows and 82 variables:
#' \describe{
#'   \item{0A Intact}{length normed counts}
#'   \item{0B Intact}{length normed counts}
#'   ...
#'   \item{240B Grafted top}{length normed counts}
#' }
#' 
"exp.data"

#' List of foldchanges 
#'
#' A dataset containing a list of foldchanges for each time point
#' of grafted and separarted samples compared to intact samples
#' 
#' @format A list with 4 elements
#' \describe{
#'   \item{fc.list$`Grafted bottom`}{matrix of foldchanges}
#'   \item{fc.list$`Grafted top`}{matrix of foldchanges}
#'   \item{fc.list$`Seperated bottom`}{matrix of foldchanges}
#'   \item{fc.list$`Seperated top`}{matrix of foldchanges}
#' }
#' 
"fc.list"

#' List of marginal likelihoods of upregulated genes compared to intact
#'
#' A dataset containing a list of marginal likelihoods for each time point
#' of grafted and separarted samples compared to intact samples
#' 
#' @format A list with 4 elements
#' \describe{
#'   \item{ml.list.up$`Grafted bottom`}{matrix of foldchanges}
#'   \item{ml.list.up$`Grafted top`}{matrix of foldchanges}
#'   \item{ml.list.up$`Seperated bottom`}{matrix of foldchanges}
#'   \item{ml.list.up$`Seperated top`}{matrix of foldchanges}
#' }
#' 
"ml.list.up"

#' List of marginal likelihoods of downregulated genes compared to intact
#'
#' A dataset containing a list of marginal likelihoods for each time point
#' of grafted and separarted samples compared to intact samples
#' 
#' @format A list with 4 elements
#' \describe{
#'   \item{ml.list.down$`Grafted bottom`}{matrix of foldchanges}
#'   \item{ml.list.down$`Grafted top`}{matrix of foldchanges}
#'   \item{ml.list.down$`Seperated bottom`}{matrix of foldchanges}
#'   \item{ml.list.down$`Seperated top`}{matrix of foldchanges}
#' }
#' 
"ml.list.down"

#' Vector containing groups for summarized samples
#'
#' 
#' @format A numeric vector indicating the relation to a certain group
#' 
#' 
"groupS"
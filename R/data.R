##' ENRON and utility sector indices
##'
##' Daily values of the ENRON stock index and of an aggregate utility sector
##' index over the 260 trading days of 1993 (4 January to 31 December).
##'
##' The data are the columns \code{ENRON.index} and \code{UTILITY.index}
##' of the data frame \code{UTIL.index} distributed with the \pkg{Rsafd}
##' package, which accompanies Carmona (2014). The two columns have been
##' renamed and the row names (dates) dropped; the third column of
##' \code{UTIL.index}, \code{DUKE.index}, is not included.
##'
##' @format A data frame with 260 rows and 2 numeric variables:
##' \describe{
##'   \item{ENRON}{Daily value of the ENRON stock index.}
##'   \item{UTILITY}{Daily value of an aggregate utility sector index.}
##' }
##' @source Carmona, R. (2025). \pkg{Rsafd}: Statistical Analysis of Financial
##'   Data in R. R package version 1.0.1,
##'   \url{https://github.com/PrincetonUniversity/Rsafd}, data frame
##'   \code{UTIL.index}. Distributed under the MIT license,
##'   Copyright (c) 2025 PrincetonUniversity; see \file{inst/COPYRIGHTS}.
##' @references Carmona, R. (2014). \emph{Statistical Analysis of Financial
##'   Data in R}, 2nd edition. Springer, New York.
##'   \doi{10.1007/978-1-4614-8788-3}
##' @keywords datasets
"ENRON"


##' Annotated age and TweetGenie age estimate for Twitter accounts
##'
##' Anonymised records of Twitter accounts from the study of Sanders et al.
##' (2016), pairing the age assigned by a human annotator with the age
##' estimated by the TweetGenie tool.
##'
##' Each row is one annotation record, not one account: accounts annotated by
##' more than one annotator appear in more than one row (Sanders et al. (2016)
##' analyse 1,815 accounts, whereas this data frame has 1,956 rows). The
##' values of \code{age} are human judgements, not verified ages.
##'
##' @format A data frame with 1956 rows and 2 variables:
##' \describe{
##'   \item{age.tweetgenie}{Numeric. Age in years estimated by TweetGenie.}
##'   \item{age}{Integer. Age in years assigned by a human annotator.}
##' }
##' @source Anonymised data provided by Eric Sanders (Radboud University),
##'   personal communication, 2 March 2026, in response to a request to include
##'   them in this package.
##' @references Sanders, E., de Gier, M., and van den Bosch, A. (2016). Using
##'   demographics in predicting election results with Twitter. In
##'   \emph{Social Informatics}, Lecture Notes in Computer Science,
##'   pp. 259--268. Springer. \doi{10.1007/978-3-319-47874-6_18}
##' @keywords datasets
"ELECTION"
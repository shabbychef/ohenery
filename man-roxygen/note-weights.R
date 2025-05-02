#' @note
#' To avoid incorrect inference when only the top
#' performers are recorded, and all others are
#' effectively tied, one should use weighting.
#' Set the weights to zero for participants who
#' are tied non-winners, and one for the rest
#' So for example, if you observe the Gold, Silver,
#' and Bronze medal winners of an Olympic event
#' that had a starting field of 12 participants,
#' set weights to 1 for the medal winners, and 0
#' for the others. Note that the weights do not
#' attach to the participants, they attach to
#' the place they took.

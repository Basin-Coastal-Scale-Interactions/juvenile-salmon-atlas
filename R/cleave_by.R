
#-----------------------------------------------------------------------------
#' tidy version of the base `split()` function
#'
#' @param df data.frame to split
#' @param ... unquoted column names from the 'df' data.frame to specify 
#'            the splitting groups
#'
#' @return list of data.frames, one for each group
#-----------------------------------------------------------------------------
cleave_by <- function(df, ...) {
  stopifnot(inherits(df, "data.frame"))
  
  # use tidyeval to get the names ready for dplyr
  grouping <- quos(...)
  
  # Calculate a single number to represent each group
  group_index <- df %>%
    group_by(!!!grouping) %>%
    group_indices()
  
  # do the split by this single group_index variable and return it
  split(df, group_index)
}

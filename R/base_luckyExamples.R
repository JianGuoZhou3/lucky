


luckyEsxamples <- function(){

  ## Parallel for plyr functions
  if(F){
    pddply <- function(.data,
                       .variables,
                       .fun = NULL,
                       ...,
                       njobs = parallel::detectCores() - 1,
                       .progress = "none",
                       .inform = FALSE,
                       .drop = TRUE,
                       .paropts = NULL){

      # Sanity check
      if (njobs < 1) {
        njobs <- 1
      }

      # Non parallel option
      if (njobs == 1) {
        return(plyr::ddply(.data, .variables, .fun = .fun, ...,
                           .progress = .progress, .inform = .inform,
                           .drop = .drop, .parallel = FALSE))
      }

      # Make sure the foreach package is available
      if (!requireNamespace("foreach", quietly = TRUE)) {
        stop("The 'foreach' package is required by pddply() when 'njobs > 1'")
      }

      # Set up the cluster
      cl <- parallel::makeCluster(njobs)
      doParallel::registerDoParallel(cl)

      # Process the results
      o <- tryCatchWE(plyr::ddply(.data, .variables, .fun = .fun, ...,
                                  .progress = .progress, .inform = .inform,
                                  .drop = .drop, .parallel = TRUE, .paropts = .paropts))

      # Shut down the cluster
      parallel::stopCluster(cl)

      # If we have an error, then stop
      if (is(o$value, "error")) {
        stop(o$value)
      }

      # Extract and remove bogus warnings from the plyr package
      if (!is.null(o$warning)) {

        # Seemingly unique strings from the bogus warning
        bogus1 <- c("<anonymous>: ... may be used in an incorrect context:")
        bogus2 <- c(".fun(piece, ...)")

        # Determine whether these strings are present in any warnings
        vWarn <- function(x) {
          return(!(grepl(bogus1, x, fixed = TRUE) & grepl(bogus2, x, fixed = TRUE)))
        }

        # Identify the valid warnings
        validWarnings <- unlist(lapply(o$warning, vWarn))

        # I don't think this portion of the code will ever be really used--
        # because as far as I can tell, both ddply() (and parLapply(), for that
        # matter) don't capture and return warnings from the worker nodes.
        # But if there any other warnings, issue them
        if (any(validWarnings)) {

          nothing <- lapply(o$warning[validWarnings], warning)

        }

      }

      # Return the result
      return(o$value)

    } # pddply
  }


  if(F){
    ## multip
    installed <- as.character(installed.packages()[,"Package"])
    if(!"multidplyr" %in% installed){
      devtools::install_github("hadley/multidplyr")
    }
  }



}
























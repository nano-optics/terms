##' @title Extract timings from log files
##' @description Parses a log file to extract timing information from subroutines (verbosity-dependent)
##' @describeIn utils extract timings from log files
##' @param log filename
##' @return returns a tibble of timings
##' @import tibble
##' @export
extract_time <- function(log = 'log'){
  command <- paste0("grep '(CPU & real in s)' ", log)
  con <- pipe(command)
  timestrings <- readLines(con)
  close(con)
  calculation <- strsplit(timestrings, ">")
  routine <- sapply(calculation, "[", 1)
  detail <- stringr::str_extract(string = sapply(calculation, "[", 2),
                                 pattern = "(?<=\\[).*(?=\\])")
  
  cpu <- sapply(strsplit(timestrings, " "), 
                function(v) as.numeric(v[length(v)-2]))
  total <- sapply(strsplit(timestrings, " "), 
                  function(v) as.numeric(v[length(v)]))
  
  tibble(routine = routine, detail = detail, cpu = cpu, total=total)
}

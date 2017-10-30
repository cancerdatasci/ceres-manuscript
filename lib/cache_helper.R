
cache_helper <- function(object_name, load_function, force=F) {
    cache_filename <- file.path("./cache",
                                str_c(ProjectTemplate:::cache.name(object_name),
                                      ".RData"))

    if (!exists(object_name) | force) {
        message("\t munging ", object_name, "...\n")
        assign(object_name, load_function(), .GlobalEnv)
        message("\t caching ", object_name, "...\n")
        cache(object_name)
    }  else {
        message("\t", object_name, " already loaded into memory. Not re-running.")
    }

}

clear_cache <- function() {
    are_you_sure <- readline(prompt="Are you sure? ")
    if (are_you_sure %in% c("Yes", "yes", "Y", "y")) {
        res <- file.remove(list.files(path="./cache/", pattern=".RData", full.names=T))
    }
}

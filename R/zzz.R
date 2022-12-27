.onLoad <- function(libname, pkgname) {
    # nocov start
    shiny::addResourcePath(
        prefix = "logo",
        directoryPath = system.file(
            "logo",
            package = "CRISPRball"
        )
    )
    # nocov end
}

.onUnload <- function(libname, pkgname) {
    # nocov start
    shiny::removeResourcePath("logo")
    # nocov end
}

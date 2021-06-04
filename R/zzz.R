# Misc functions

# warning about invalid PPS objects in version <0.12.0
.onAttach <- function(libname, pkgname) {
#   packageStartupMessage(
#     paste0('\nNOTE: The PAMpalSettings object has been changed in ',
#            'version 0.12.0. You will get a "Not a validObject()" warning',
#            ' for PAMpalSettings and AcousticStudy objects created prior to ',
#            '0.12.0, this can be fixed with the "updatePamObject() function:\n\n',
#            'x <- updatePamObject(x)')
#   )
}
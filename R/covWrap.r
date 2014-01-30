covWrap <- function(s1, s2) {
    out <- .Call("covWrap", s1, s2)
    return(out)
}

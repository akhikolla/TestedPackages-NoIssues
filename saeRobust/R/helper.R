addAttr <- function(to, what, name) {
    attr(to, name) <- what
    to
}

getter <- function(.x, .wrapper = identity) {
    memoise::memoise(function() {
        .wrapper(.x)
    })
}

storage <- module({
    reformat <- function(storage) {
        .iter <- function(ind) {
            lapply(storage, function(.) attr(.[[ind]], "history")) %>%
            { lapply(seq_along(.), function(i) cbind(.[[i]], i)) } %>% # add iteration
                do.call(what = rbind)
        }

        lapply(seq_along(storage[[1]]), .iter)
    }
})

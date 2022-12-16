data_list <- list(
  linch = jsonlite::fromJSON("raw_data/linch-zhang.json"),
  finn = jsonlite::fromJSON("raw_data/finn-moorhouse.json"),
  gavin = jsonlite::fromJSON("raw_data/gavin-leech-complete.json"),
  jamie = jsonlite::fromJSON("raw_data/jaime-sevilla.json"),
  misha = jsonlite::fromJSON("raw_data/misha-yagudin.json"),
  ozzie = jsonlite::fromJSON("raw_data/ozzie-gooen.json")
)

lengths <- sapply(data_list, nrow)
new <- do.call(rbind, lapply(data_list, \(x) data.frame(
  source = x$source,
  target = x$target,
  distance = x$distance
)))
new$rater <- ""
i <- 1
for (j in seq_along(lengths)) {
  new$rater[i:(lengths[j] + i - 1)] <- names(lengths[j])
  i <- lengths[j] + i
}

rownames(new) <- NULL
new$source <- as.factor(new$source)
new$target <- as.factor(new$target)
dat.research <- new
save(dat.research, file = "data/dat.research.rda")

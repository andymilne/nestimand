## demonstration data, exactly as in partially_nested_variables.md Section 1
set.seed(1); n <- 40
cells <- data.frame(
  chord_type = c(rep(c("dim", "min", "maj"), each = 3), "aug"),
  inversion  = c(rep(c("0", "1", "2"), times = 3), "none"),
  mu         = c(4.16, 4.00, 3.73,
                 4.29, 4.20, 3.71,
                 4.98, 4.62, 4.39,
                 3.92))
dat <- cells[rep(seq_len(nrow(cells)), each = n), ]
dat$participant <- factor(rep(seq_len(n), times = nrow(cells)))
tr <- runif(n, 0, 10)
dat$training <- tr[as.integer(dat$participant)]
dat$response <- rnorm(nrow(dat), dat$mu + 0.12 * dat$training, 1.2)
dat$chord_type <- factor(dat$chord_type, levels = c("aug", "dim", "min", "maj"))
dat$inversion  <- factor(dat$inversion,  levels = c("none", "0", "1", "2"))

# g1_glmm_v05.R -- primary G1 test (protocol v0.5, frozen).
# Binomial GLMM: miscover ~ model * rank(D) + log(n_cells) + (1 | target gene) + (1 | split)
# primary score S1_k50, alpha = 0.1, models mean (reference), linear_pca, linear_go.
# Usage: Rscript g1_glmm_v05.R <ROOT>
suppressPackageStartupMessages({library(lme4); library(arrow); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else "."
an <- file.path(root, "results", "g1_v0.5", "analysis")
d0 <- as.data.frame(read_parquet(file.path(an, "glmm_input.parquet")))
d0 <- d0[d0$model %in% c("mean", "linear_pca", "linear_go"), ]
d0$model <- factor(d0$model, levels = c("mean", "linear_pca", "linear_go"))
d0$target <- factor(d0$target); d0$split <- factor(d0$split)

fit_one <- function(audit, dist, with_ncells) {
  d <- d0[d0$audit == audit, ]
  d$rk <- d[[paste0("rank_", dist)]]
  f <- if (with_ncells) miscover ~ model * rk + log_ncells_c + (1 | target) + (1 | split)
       else miscover ~ model * rk + (1 | target) + (1 | split)
  msgs <- character(0)
  fit <- withCallingHandlers(
    glmer(f, data = d, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
    warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") })
  b <- fixef(fit); V <- as.matrix(vcov(fit))
  lc <- function(w) { est <- sum(w * b); se <- sqrt(as.numeric(t(w) %*% V %*% w)); list(est = est, se = se, p = 2 * pnorm(-abs(est / se))) }
  e <- function(nm) { w <- setNames(rep(0, length(b)), names(b)); w[nm] <- 1; w }
  slopes <- list(mean = lc(e("rk")),
                 linear_pca = lc(e("rk") + e("modellinear_pca:rk")),
                 linear_go = lc(e("rk") + e("modellinear_go:rk")))
  diffs <- list(`linear_pca-mean` = lc(e("modellinear_pca:rk")), `linear_go-mean` = lc(e("modellinear_go:rk")))
  vc <- as.data.frame(VarCorr(fit))
  list(n = nrow(d), n_genes = nlevels(droplevels(d$target)), n_splits = nlevels(droplevels(d$split)),
       slopes = slopes, diffs = diffs, fixef = as.list(b),
       re_sd = setNames(as.list(vc$sdcor), vc$grp), singular = isSingular(fit),
       messages = c(msgs, fit@optinfo$conv$lme4$messages))
}

fits <- list()
for (audit in c("A", "B")) {
  for (dist in c("D_GO", "D_STRING", "D_GO_filt")) {
    key <- paste(audit, dist, "ncells", sep = "_"); t0 <- Sys.time()
    fits[[key]] <- fit_one(audit, dist, TRUE)
    cat(key, round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
  }
  key <- paste(audit, "D_GO_no_ncells", sep = "_")
  fits[[key]] <- fit_one(audit, "D_GO", FALSE); cat(key, "done\n")
}
out <- list(protocol_version = "0.5", lme4 = as.character(packageVersion("lme4")), R = R.version.string, fits = fits)
write(toJSON(out, auto_unbox = TRUE, digits = NA, null = "null"), file.path(an, "glmm_results.json"))
cat("GLMM complete\n")

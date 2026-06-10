require(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript estimate_bulk_min_cov.R <input_file> <output name>")
}

input_file <- args[1]
output_name <- args[2]
min_pct <- 90
# alpha <- 0.05
alpha <- 0.1 # we want to be more conservative in our estimate of the minimum coverage threshold, to avoid picking a threshold that is too high and results in too few positions being retained. By setting alpha to 0.1, we are saying that we want to pick the highest point where we have at least 90% of the plateau, which means that we are willing to accept a lower threshold to study more positions. We would need a panel-of-normals and strict gnomad overlap to remove any remaining germline variants, but this is a good starting point for our analysis. In practice, we get very few, due to the extensive masks from Lawson et al. 
min_threshold <- 6

cov_test <- read.table(file = input_file,header = FALSE,sep = "\t") %>% as_tibble()
colnames(cov_test) <- c("batch","min_cov","n_positions")
print(head(cov_test))
pct <- cov_test %>% group_by(min_cov) %>% summarize(n_pos = sum(n_positions) * 100) %>% mutate(pct_of_all_pos = 100*n_pos/max(n_pos))
pct_nonZero <- pct %>% dplyr::filter(min_cov > 0)
print(pct_nonZero)
est_pct_logistic <- nls(pct_of_all_pos ~ SSlogis(min_cov,Asym,xmid,scal),pct_nonZero)
pct_nonZero$pred <- predict(object = est_pct_logistic,data.frame(min_cov = pct_nonZero$min_cov))

# pick the plateau point where the at least 1-alpha percentage of the plateau are covered. The plateau is the maximum percentage of positions that can have coverage >= min_cov, as estimated by the logistic curve. We want to pick the highest point where we have at least 1-alpha percentage of the plateau, because we want to be conservative in our estimate of the minimum coverage threshold, to avoid picking a threshold that is too high and results in too few positions being retained.
plateau <- min(100,coefficients(est_pct_logistic)["Asym"]) # the plateau of the logistic curve, i.e. the maximum percentage of positions that can have coverage >= min_cov, as estimated by the logistic curve. if this is < min_pct, then we will never reach min_pct, and we should just pick the min_cov that corresponds to the plateau.
# we want to pick the highest point where we have at least 100*(1-alpha)% of the plateau
max_pct <- (1 - alpha) * plateau
chosen_threshold <- pct_nonZero %>% dplyr::filter(pred > max_pct) %>% slice_max(min_cov)
if (nrow(chosen_threshold) == 0) {
  print(paste("No threshold found that meets the criteria. Choosing the minimum coverage threshold such that we can get",min_pct,"percent of bases"))
  # print(paste("No threshold found that meets the criteria. Reverting to default threshold of",min_threshold))
  chosen_threshold <- pct_nonZero %>% dplyr::filter(min_cov == min_threshold)
}

# save model
out_table_file <- paste0(output_name,".logistic_fit.rds")
saveRDS(object=est_pct_logistic,file=out_table_file)
# write table
out_table_file <- paste0(output_name,".logistic_pred.txt")
write.table(pct_nonZero,out_table_file,sep = "\t",col.names = TRUE,row.names = FALSE,quote=FALSE)
# write chosen threshold
out_threshold_file <- paste0(output_name,".threshold.txt")
write.table(chosen_threshold,out_threshold_file,sep = "\t",col.names = TRUE,row.names = FALSE,quote=FALSE)

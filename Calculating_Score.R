tidy_input <- function(file) {
  surv_expr <- fread(file, data.table = F)
  surv_df <- surv_expr %>% dplyr::select(c("sample", "status", "time", sel))
  # surv_df$time <- surv_df$time / 30
  surv_expr <- surv_df
  names(surv_expr)[2:3] <- c("status", "time")
  surv_expr %<>% column_to_rownames("sample")
  sur <- surv_expr %>% dplyr::select(1:2)
  expr <- surv_expr %>% dplyr::select(-c(1:2))
  expr <- zscore(expr) %>% as.data.frame()
  expr <- expr %>% mutate_all(~replace(., is.nan(.), 0))
  res <- cbind(sur, expr)
  return(res)
}

TCGA <- tidy_input("TCGA_surv_expr.txt")
mod.TCGA <- RSF.mod(TCGA)
save(mod.TCGA, file = "model.RData")

pred.TCGA <- RSF.pred(mod.TCGA, TCGA)
df <- pred.TCGA %>% rownames_to_column("sample")
cindex(df)
write.csv(df, "score_surv_TCGA.csv", row.names = F)
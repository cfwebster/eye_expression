#CWebster - quantum catch modeling
# only using two time points as examples (one from seasonal, one from daily)
library(pavo)
fall_spec #from dataset downloaded from external source of irradiance values, from 300-700nm
colnames(fall_spec)[1] <- "wl"

#converting to 1 nm increments, but need to return as not rspec objects
fall_illum <- as.rspec(fall_spec, interp = TRUE, lim = c(300, 700))
#setting up opsin sensitivity peaks, based on literature where SW = 355, RHO = 500, MW/LW = 555
absorb <- sensmodel(c(355, 500, 555), beta = TRUE)
wl <- absorb$wl
illum_fall <- fall_illum[, 2]
names(illum_fall) <- fall_illum$wl

#removing wl column for absorbance
absorb_only <- absorb[, -1]
#multiplying opsin sensitivity by the irradiance
qc_fall_raw <- colSums(absorb_only * illum_fall)

#opsin expression, getting data from DESeq2 results
norm_counts <- counts(eye_dds, normalized = TRUE)
opsin_genes <- c("ENST00000249389.OPN1SW.14", "ENST00000599405.OPN1MW3.4347_evm.model.CM061280.346", "ENST00000296271.RHO.160")
opsin_counts <- norm_counts[rownames(norm_counts) %in% opsin_genes, ]
fall_samples <- colnames(eye_dds)[eye_dds$season == "fall"]
fall_expr <- rowMeans(opsin_counts[, fall_samples])
fall_expr_scaled <- fall_expr / sum(fall_expr)
names(fall_expr_scaled) <- c("OPN1SW", "RHO", "OPN1MW/LW")

#to get relative opsin quantum catch values
qc_fall_combined <- qc_fall_raw * fall_expr_scaled
qc_fall_rel <- qc_fall_combined / sum(qc_fall_combined)
qc_fall_individ_rel <- sapply(fall_samples, function(ind) {
  expr <- opsin_counts[, ind]
  expr_scaled <- expr / sum(expr)
  qc <- qc_fall_raw * expr_scaled
  qc / sum(qc)
})
qc_fall_individ_rel <- t(qc_fall_individ_rel)

#to get total quantum catch for each individual, combining all opsins
qc_fall_individuals <- sapply(fall_samples, function(ind) {
  expr <- opsin_counts[, ind]
  sum(qc_fall_raw * expr)
})
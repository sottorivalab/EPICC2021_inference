get_sample_annotation__ = function(samples, attribute, silent=FALSE) {

  results =
    do.call(rbind, lapply(samples, function(sample) {
      tryCatch({

        annotation = 
          tryCatch({
            annotation_from_barcode(sample, TRUE)},
            error=function(e) {return(NULL)
          })

        # invalid annotation
        if (is.null(annotation)) {
          return(data.frame(value = NA, reason = "Invalid barcode."))
        }

        # sequenza fits (WGS)
        if (exists("sequenza_purity_ploidy_fits_")) {
          if (annotation$sample_barcode %in% rownames(sequenza_purity_ploidy_fits_)) {
            value = sequenza_purity_ploidy_fits_[annotation$sample_barcode, attribute]
            return(data.frame(value = value, reason = NA))
          }
        }

        # genotyping data (LP), only purity
        if (exists("lp_genotyping_data") & attribute == "purity" ) {
          if (annotation$sample_barcode %in% rownames(lp_genotyping_data)) {
            value = lp_genotyping_data[annotation$sample_barcode, attribute]
            return(data.frame(value = value, reason = NA))
          }
        }

        # lowpass fits (LP)
        if (annotation$sample_barcode %in% rownames(lowpass_ploidy_fits)) {
          value = lowpass_ploidy_fits[annotation$sample_barcode, attribute]
          return(data.frame(value = value, reason = NA))
        }

        # throw error, caught below.
        stop("")

      }, error=function(e) {
        return(data.frame(value = NA, reason = "Missing data."))
      })
    }))

  results = cbind(sample=samples, results)

  wh_missing = is.na(results$value)

  if (any(wh_missing) & !silent) {
    cat(sprintf("Failed to get %s estimates for the following samples:\n\n", attribute))
    print(results[wh_missing,c("sample","reason")])
  }

  return(results$value)
}

get_purity = function(samples, silent=FALSE) {
  get_sample_annotation__(samples, "purity", silent)
}

get_ploidy = function(samples, silent=FALSE) {
  get_sample_annotation__(samples, "ploidy", silent)
}

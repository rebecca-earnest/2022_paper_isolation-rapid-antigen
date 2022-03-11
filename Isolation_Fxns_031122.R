
# extracts ggplot legend
get_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend) 
}

# plots survival curves for each covariate
plot_cov_KM <- function(cov_val) {
  
  fit_cov <- survfit(surv_obj ~ cov_val)
  legend_vals <- levels(cov_val)
  num_levels <- length(legend_vals)
  plot(fit_cov, xlab = "Isolation Duration", ylab = "Probability of survival (testing pos)",
       lwd = 2, lty = 1:num_levels)
  
  if (num_levels == 2) {
    legend("topright", bty = "n", lty = 1:num_levels, c(legend_vals[1], legend_vals[2]), lwd = 2)
    
  } else {
    legend("topright", bty = "n", lty = 1:num_levels, c(legend_vals[1], legend_vals[2], legend_vals[3]), lwd = 2)
  }
}

# creates survival analysis results table
build_results_table <- function(dat, dat_all, col_names) {
  
  for (each_col in col_names) {
    
    dat_col <- dat[, each_col]
    if (each_col != "p_vals") {
      dat_col <- round(dat_col, digits = 2)
    } else {
      dat_col <- round(dat_col, digits = 3)
      
    }
    dat_all[which(dat_all$ref_col != "ref"), each_col] <- dat_col
    dat_all[is.na(dat_all)] <- "-"
  }
  dat_all
}

# from AFTtools R package; generates QQ plots from survreg object
AFTplot <- function(model) {
  
  assertive::assert_is_all_of(model, "survreg")
  
  
  
  par_set <- par(mar = c(4, 4, 0.5, 0.5),
                 oma = c(0, 0, 2, 0),
                 mfrow = c(2, 2))
  
  
  
  quantPredict <- cbind(model.frame(model),
                        Q = predict(model, type = "quantile", p = c(1:49 / 50)))
  
  
  
  fixEf <- names(model$xlevels) %>%
    {.[which(!grepl("frailty.", .))]}
  
  if (is.null(fixEf)) {
    
    stop("Only one fixed effect level found in model. Do not know how to produce a Q-Q-plot.")
    
  }
  
  ii <- 0
  for (i in seq_along(fixEf)) {
    
    vals <- lapply(model$xlevels[[i]], function(x) {
      
      subset(quantPredict, quantPredict[, fixEf[i]] == x) %>%
        {.[, grep("Q.", names(.))]}
      
    })
    
    cbs <- combn(length(model$xlevels[[i]]), 2)
    for (j in 1:ncol(cbs)) {
      
      qqplot(x = as.matrix(vals[[cbs[1, j]]]),
             y = as.matrix(vals[[cbs[2, j]]]),
             xlab = paste(fixEf[i], "==", model$xlevels[[i]][cbs[1, j]]),
             ylab = paste(fixEf[i], "==", model$xlevels[[i]][cbs[2, j]]))
      
      XX <- quantile(as.matrix(vals[[cbs[1, j]]]), probs = seq(0, 1, 0.1))
      YY <- quantile(as.matrix(vals[[cbs[2, j]]]), probs = seq(0, 1, 0.1))
      bb <- lm(YY ~ 0 + XX)$coefficients
      
      stopifnot(length(bb) == 1)
      
      abline(a = 0, b = bb, lty = 2)
      
      if (assertive::is_whole_number(ii / 4)) {
        
        title(main = "Q-Q-Plot of predicted survival times",
              outer = TRUE)
        
      }
      ii <- ii + 1
      
    }
    
  }
  
  par(par_set)
  
}
  


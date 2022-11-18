

count.mean.sd.fun <- function(data) {
    # Calculate the count, median, mean, min, max ,sd, cv for the given data
    #
    # Args:
    # 	data: Numeric data table
    #
    # Returns:
    # 	Calculated values will be appended at the end of data table

    # doc.table <- flex.table.fun(data)

    n <- c() # count
    med <- c() # median
    min <- c() # min
    max <- c() # max
    mean <- c() # mean
    sd <- c() # standard deviation
    cv <- c() # coefficient variation

    for (i in 1:(ncol(data) - 1)) {
      n[i] <- length(data[data[,1+i]!="NA", 2])
       # n[i] <- sum(!is.na(data[, 1 + i]))
        med[i] <- median(as.numeric(data[, 1 + i]), na.rm = TRUE)
        min[i] <- min(as.numeric(data[, 1 + i]), na.rm = TRUE)
        max[i] <- max(as.numeric(data[, 1 + i]), na.rm = TRUE)
        # mean[i] <- round(sum(data[,1+i], na.rm=TRUE)/length(data[data[,1+i]!="NA",2]),2)
        mean[i] <- round(mean(as.numeric(data[, 1 + i]), na.rm = TRUE), 2)
        # sd[i]   <- round(sd(data[,1+i], na.rm=TRUE),2)
        sd[i] <- round(sd(as.numeric(data[, 1 + i]), na.rm = TRUE), 2)
        cv[i] <- round(sd(as.numeric(data[, 1 + i]), na.rm = TRUE) / mean(as.numeric(data[, 1 + i]), na.rm = TRUE), 2)
    }

    all.list <- list(n, med, min, max, mean, sd, cv)
    all.names <- c("N", "Median", "Min", "Max", "Mean", "SD", "CV")

    for (j in 1:length(all.list)) {
        all.list[[j]][all.list[[j]] == Inf] <- "-"
        all.list[[j]][all.list[[j]] == -Inf] <- "-"
        all.list[[j]][is.na(all.list[[j]])] <- "-"
        all.list[[j]][is.nan(all.list[[j]])] <- "-"
        all.list[[j]][all.list[[j]] == NaN] <- "-" #updated 22.10.04.
    }
    return(list(list_ = all.list, names_ = all.names))
}

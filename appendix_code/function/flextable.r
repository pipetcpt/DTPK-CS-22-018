
"%ni%" <- Negate("%in%")

flex.table.fun <-
    function(data,
             head.font = 9,
             body.font = 9,
             calc_ = FALSE) {
        # Convert data frame into FlexTable format
        #
        # Args:
        #	data: Merged data table
        #	head.font: column header font size (default is 9)
        #	body.font: body content font size  (default is 9)
        #
        # Returns:
        #	The flextable format of dataset
        
        if ("SUBJID" %in% names(data)) {
            data$SUBJID <- as.character(data$SUBJID)
        }
        big_border = fp_border(color = "black", width = 2)
        border_v = fp_border(color = "black")
        border_h = fp_border(color = "black")
        data[is.na(data)] <- "NA"
        data[data == ""] <- "NA"
        data[data == " "] <- "NA"
        if (calc_ == TRUE) {
            print("Calculating the count, mean, sd, max, min, etc..")
            count.mean.sd.list <- count.mean.sd.fun(data)
            for (i in 1:length(count.mean.sd.list$list_)) {
                data <-
                    rbind(data, c(
                        count.mean.sd.list$names_[i],
                        unlist(count.mean.sd.list$list_[i])
                    ))
            }
            data[is.na(data)] <- "NA"
            MyFTable <- flextable(data) %>%
                bold(part = "header") %>%
                fontsize(part = "header", size = head.font) %>%
                fontsize(part = "body", size = body.font) %>%
                align(align = "center", part = "all") %>%
                border_outer(part = "all", border = big_border) %>%
                border_inner_h(part = "all", border = border_h) %>%
                border_inner_v(part = "all", border = border_v) %>%
                bg(
                    i = (nrow(data) - 6):nrow(data),
                    bg = "#DDDDDD",
                    part = "body"
                ) %>%
                autofit()
            
        } else{
            if ("Age" %in% names(data) && c('Height', 'Weight') %ni% names(data)) {
                MyFTable <- flextable(data) %>%
                    colformat_num(, digits = 0, col_keys = 'Age') %>%
                    #colformat_num(, digits = 1, col_keys = c("Height", "Weight")) %>%
                    bold(part = "header") %>%
                    fontsize(part = "header", size = head.font) %>%
                    fontsize(part = "body", size = body.font) %>%
                    align(align = "center", part = "all") %>%
                    border_outer(part = "all", border = big_border) %>%
                    border_inner_h(part = "all", border = border_h) %>%
                    border_inner_v(part = "all", border = border_v) %>%
                    autofit()
            } else if ('Height' %in% names(data) && 'Weight' %in% names(data)) {
                MyFTable <- flextable(data) %>%
                    colformat_num(, digits = 0, col_keys = 'Age') %>%
                    colformat_num(,digits = 1,col_keys = c("Height", "Weight")) %>%
                    bold(part = "header") %>%
                    fontsize(part = "header", size = head.font) %>%
                    fontsize(part = "body", size = body.font) %>%
                    align(align = "center", part = "all") %>%
                    border_outer(part = "all", border = big_border) %>%
                    border_inner_h(part = "all", border = border_h) %>%
                    border_inner_v(part = "all", border = border_v) %>%
                    autofit()
            } else {
                MyFTable <- flextable(data) %>%
                    bold(part = "header") %>%
                    fontsize(part = "header", size = head.font) %>%
                    fontsize(part = "body", size = body.font) %>%
                    align(align = "center", part = "all") %>%
                    border_outer(part = "all", border = big_border) %>%
                    border_inner_h(part = "all", border = border_h) %>%
                    border_inner_v(part = "all", border = border_v) %>%
                    autofit()
            }
        }
        return (MyFTable)
    }
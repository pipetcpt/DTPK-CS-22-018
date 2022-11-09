

create.table <-
    function(csvfile,
             code_id = NULL,
             fullname = NULL,
             filter_var = NULL,
             value_var,
             period_ = NULL,
             type_ = FALSE) {
        # Create the pivot table for infile dataset
        #
        # Args:
        # 	infile: the numeric dataset to be used to create the pivot table
        # 	code_id: specific code id that is required to filter out from the infile dataset
        # 	fullname: the full name of code_id
        # 	filter_var: the filter column name
        # 	value_var: the variable needed to be put as value
        # 	period_: whether period 1 or period 2
        #
        # Returns:
        # 	Pivot table for the dataset

        # Match with VISIT CODE
        csvfile <- visit.match.fun(infile = csvfile, visitData = VST)

        if (!is.null(code_id) && !is.null(filter_var)) {
            if (!is.null(period_)) {
                if (period_ == "P1") {
                    data <- csvfile %>%
                        mutate(Folder2 = substr(Folder, 1, 2)) %>%
                        filter((
                            Folder2 == as.character(period_) | Folder2 == "SC"
                        ) &
                            (csvfile[filter_var] == as.character(code_id))) %>%
                        select(
                            #changed subject -> SID/ #add RID
                            # updated 22.02.03
                            SID,RID,
                            InstanceRepeatNumber,
                            FolderSeq,
                            value_var
                        ) %>%  # updated 22.02.03
                        arrange(SID,RID, InstanceRepeatNumber)
                } else if (period_ == "P2") {
                    data <- csvfile %>%
                        mutate(Folder2 = substr(Folder, 1, 2)) %>%
                        filter(Folder2 == as.character(period_) &
                            csvfile[filter_var] == as.character(code_id)) %>%
                        select(  # updated 22.02.03
                            SID,RID,
                            InstanceRepeatNumber,
                            FolderSeq,
                            value_var
                        ) %>%
                        arrange(SID,RID, InstanceRepeatNumber)  # updated 22.02.03
                }
            } else {
                data <- csvfile %>%
                    filter(csvfile[filter_var] == as.character(code_id)) %>%
                    select(SID,RID, VISIT, value_var) %>%
                    arrange(SID,RID, VISIT) # updated 22.02.03
            }
        } else {
            data <- csvfile %>%
                select(SID,RID, InstanceRepeatNumber, FolderSeq, value_var) %>%
                arrange(SID,RID, InstanceRepeatNumber) # updated 22.02.03
        }

        if (type_ == "numeric") {
            data[value_var][[1]] <-
                as.numeric(as.character(data[value_var][[1]]))
        }

        # Remove duplicated row by Subject and FolderSeq
        # data <- remove.dup.fun(data)

        # Create Pivot table by value_var
        # data.pivot <-
        #     dcast(data,
        #           SUBJID ~ factor(VISIT, levels = str_sort(unique(VISIT), numeric = T)),
        #           value.var = as.character(value_var))
        distinct.data <- data %>% distinct()
        data.pivot <- spread(distinct.data, key = "VISIT", value = value_var)

        # Change Column order by Screening
        data.pivot <- data.pivot %>%
            select(SID, RID, starts_with("Scr"), starts_with("투여"), everything()) # updated 22.10,04

        # Map to full name
        if (!is.null(fullname)) {
            colnames(data.pivot)[-1] # <-
                #paste0(as.character(fullname), "_", colnames(data.pivot)[-1])
        }

        data.table <- data.pivot
        # data post-processing (Can be consistently added)
        data.table[is.na(data.table)] <- "NA"
        data.table[data.table == ""] <- "NA"

        data.table[data.table == "1 - 3"] <- "1-3"
        data.table[data.table == "0 - 1"] <- "0-1"
        # data.table[data.table=="10?? 19??"] <- "10-19"
        data.table[data.table == "negative"] <- "Negative"
        data.table[data.table == "Neg"] <- "Negative"

        return(data.table)
    }

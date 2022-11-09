
is.defined <- function(obj_) {
    # Check if the object exists or not
    # Return: TRUE if exists, FALSE if not

    obj_ <- deparse(substitute(obj_))
    # env <- parent.frame()
    exists(obj_)
    # exists(obj_, env)
}



visit.match.fun <- function(infile, visitData) {
    for (i in 1:nrow(infile)) {
        for (j in 1:nrow(visitData)) {
            if (infile$VISIT[i] == as.numeric(visitData$VISIT_CODE[j])) {
                infile$VISIT[i] <- as.character(visitData$VISIT_LABEL[j])
            }
        }
    }
    return(infile)
}


# remove.dup.fun <- function(data){
#     # Remove the duplicated row by Subject and FolderSeq

#     dup.idx <- c()
#     for(i in 2:nrow(data)){
#         if(data[i,]$Subject==data[i-1,]$Subject && data[i,]$FolderSeq==data[i-1,]$FolderSeq){
#             dup.idx <- c(dup.idx, i)
#         }
#     }
#     if(!length(dup.idx)==0){
#         data <- data[-(dup.idx-1), ]
#     }else{
#         return(data)
#     }
# }

# loop code for change Class annotation from 1,0 to pos,neg after pCRE-------
library(dplyr)
# library(pgirmess)
library(stringr)

dir <- './Data/ML'
folder <- list.dirs(dir,
                    full.names = FALSE,
                    recursive = FALSE)
middle <- 'neg'
ends <- '.txt'



for (j in 1:length(folder)) {
  folder_name <- "Test"
  for ( i in 1:40){
    number <- i
    df <- read.delim(paste0(dir, '/',folder_name,'/',
                            paste( "data_Test", number, sep = '_'),ends))
    df$Class <- str_replace_all(df$Class, '1', 'pos') # class 1 replace with postive
    df$Class <- str_replace_all(df$Class, '0', 'neg')# class 0 replace with negative
    
    # output result
    write.table(df, paste0(dir, '/',folder_name,'/',
                       paste(middle,"data_Test", number, sep = '_'),ends),
    row.names = F,
    quote = F,
    sep = '\t')

    }
}

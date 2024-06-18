# loop code for change Class annotation from 1,0 to pos,neg after pCRE-------
library(dplyr)
# library(pgirmess)
library(stringr)
library(argparse)

dir <- './Data/Kmer/Train'
folder <- list.dirs(dir,
                    full.names = FALSE,
                    recursive = FALSE)
middle <- 'neg'
ends <- '.txt.fa.pcre_df_p0.01.txt'

parser <- ArgumentParser()
parser$add_argument(
  "-TN_filename",
  help="Your TN data file name (NOT PATH)",
  required=TRUE
)
args <- parser$parse_args()
tn_filename <- args$TN_filename

for (j in 1:length(folder)) {
  folder_name <- folder[j]
  for ( i in 1:40){
    number <- i
    df <- read.delim(paste0(dir, '/',folder_name,'/',
                            paste( paste(tn_filename, "_Train", sep = ""), number, sep = '_'), ends))
    df$Class <- str_replace_all(df$Class, '1', 'pos') # class 1 replace with postive
    df$Class <- str_replace_all(df$Class, '0', 'neg')# class 0 replace with negative
    
    # output result
    write.table(df, paste0(dir, '/',folder_name,'/',
                       paste(middle,paste(tn_filename, "_Train", sep = ""), number, sep = '_'), ends),
    row.names = F,
    quote = F,
    sep = '\t')
    }


}



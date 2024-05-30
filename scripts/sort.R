library(data.table)

# take the output of bedtools coverage with argument -d. For each element in the query, calculate its mean coverage (column 9).

file <- "intact_cenh3.coverage"

# read the file
dt <- fread(file, header = FALSE)

# make new column of bed name (Name=???;)
dt[, element := sub(".*Name=([^;]+);.*", "\\1", V9)]

# summarise the mean coverage for each element
dt <- dt[, .(mean_coverage = mean(V11)), by = element]

fwrite(dt, "intact_cenh3.coverage.mean", sep = "\t", quote = FALSE, row.names = FALSE)

#!/usr/bin/env Rscript

# This script performs Kendall tests on supplied input



### Environment set-up ---------------------------------------------------

# Start with clean environment
rm(list = ls())

# Define current time
date_str <- format(Sys.time(), "%Y-%m-%d")
print(paste("Start", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# Define the current user
user <- Sys.info()[["user"]]



### List of required packages -----------------------------------------------

# Load required packages
library("magrittr")
library("optparse")

# User specified arguments
my_data_file <- "H:/data/Synechocystis_6frame/2022-01-27_Copy_numbers/Scy004_copy_numbers_norm.txt"
id_cols <- "Proteins"



### Data import ----------------------------------------------------------

my_data <- data.table::fread(
    input = my_data_file, sep = "\t", quote = "", header = TRUE,
    stringsAsFactors = FALSE, colClasses = "character", data.table = FALSE)

my_data[, colnames(my_data) != id_cols] %<>%
    lapply(., as.double) %>%
    as.data.frame(.)

my_data[is.na(my_data) | my_data == 0] <- NA

if (any(duplicated(my_data[[id_cols]]))) {
    stop("Cannot have duplicate IDs!")
}

my_data_format <- my_data
rownames(my_data_format) <- my_data_format[[id_cols]]
my_data_format[[id_cols]] <- NULL
my_data_format %<>%
    t(.) %>%
    as.data.frame(.)

#my_data_format %<>%
#    log10(.)



### Regression analysis --------------------------------------------------

all_combi_list <- expand.grid(
    colnames(my_data_format), colnames(my_data_format)) %>%
    set_colnames(c("response", "terms")) %>%
    split(x = ., f = rep(x = 1:ceiling(nrow(.)/10000), each = 10000))

for (i in 1:length(all_combi_list)) {
    
    all_combi <- all_combi_list[[i]]
    
    my_regression <- apply(X = all_combi, MARGIN = 1, FUN = function(x) {
        
        model <- try(
            expr = lm(
                formula = as.formula(
                    paste(x[["response"]], x[["terms"]], sep = "~")),
                data = my_data_format,
                na.action = "na.omit"),
            silent = TRUE)
        res <- summary(model, correlation = TRUE)
        
        if (
            x[[1]] == x[[2]] ||
            class(model) == "try-error" ||
            class(res) == "try-error" ||
            nrow(res$coefficients) < 2) {
            data.table::data.table(
                Estimate = NA,
                `Std. Error` = NA,
                `t value` = NA,
                `Pr(>|t|)` = NA,
                sigma = NA,
                df = NA,
                fstatistic = NA,
                r.squared = NA,
                adj.r.squared = NA,
                correlation = NA)
        } else {
            data.table::data.table(
                Estimate = res$coefficients[2, 1],
                `Std. Error` = res$coefficients[2, 2],
                `t value` = res$coefficients[2, 3],
                `Pr(>|t|)` = res$coefficients[2, 4],
                sigma = res$sigma,
                df = paste0(res$df, collapse = ";"),
                fstatistic = res$fstatistic[[1]],
                r.squared = res$r.squared,
                adj.r.squared = res$adj.r.squared,
                correlation = res$correlation[2, 1])
        }
        
    }) %>%
        plyr::ldply(., dplyr::bind_rows, .id = NULL) %>%
        dplyr::bind_cols(all_combi, .)
    
    add_header <- FALSE
    if (!file.exists(paste0(my_data_file, ".Regression_lm.txt"))) {
        add_header <- TRUE
    }
    
    data.table::fwrite(
        x = my_regression,
        file = paste0(my_data_file, ".Regression_lm.tmp", i),
        quote = FALSE, sep = "\t",
        row.names = FALSE, col.names = add_header,
        append = TRUE)
    
}

lm_files <- list.files(
    path = dirname(my_data_file),
    pattern = ".Regression_lm.tmp", full.names = TRUE)

my_regression <- lapply(X = lm_files, function(x) {
    data.table::fread(
        input = x, sep = "\t", quote = "", header = TRUE,
        stringsAsFactors = FALSE, colClasses = "character",
        data.table = TRUE)
}) %>%
    plyr::ldply(., data.table::data.table, .id = NULL) %>%
    dplyr::mutate(
        ., BH = p.adjust(p = `Pr(>|t|)`, method = "BH"),
        Bonferroni = p.adjust(p = `Pr(>|t|)`, method = "bonferroni"))



### END ------------------------------------------------------------------

data.table::fwrite(
    x = my_regression,
    file = paste0(my_data_file, ".Regression_lm.txt"),
    quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE,
    append = FALSE)

save.image(paste0(my_data_file, ".Regression_lm_session.RData"))

# Define end time
print(paste("Complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))



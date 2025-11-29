# Load KEGGREST
if (!requireNamespace("KEGGREST", quietly = TRUE)) install.packages("KEGGREST")
library(KEGGREST)

# 📌 Load KO IDs from CSV

file_path <- "down_trimmed.csv"  
df <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)

# Function to get KEGG pathway names
get_pathway_names <- function(KEGG_ko) {
  ko_split <- unlist(strsplit(KEGG_ko, ","))  
  pathways <- c()
  
  for (ko in ko_split) {
    ko <- trimws(ko)  
    
    # Fetch pathway links
    ko_pathways <- tryCatch({
      keggLink("pathway", paste0("ko:", ko))
    }, error = function(e) {
      message(paste("Error fetching pathways for:", ko))
      return(NULL)
    })
    
    # Fetch pathway names if found
    if (!is.null(ko_pathways) && length(ko_pathways) > 0) {
      pathway_ids <- unique(ko_pathways)
      
      pathway_names <- tryCatch({
        Sys.sleep(2)  # Increase delay to avoid rate limiting
        sapply(pathway_ids, function(pw) {
          tryCatch(keggList(pw), error = function(e) NA)
        })
      }, error = function(e) {
        message(paste("Error fetching pathway names for:", ko))
        return(NULL)
      })
      
      pathways <- c(pathways, pathway_names)
    }
  }
  
  if (length(pathways) == 0 || all(is.na(pathways))) {
    return("No pathway found")
  }
  
  return(paste(unique(na.omit(pathways)), collapse = "; "))
}

# Apply function to all KO entries
df$Pathway_Name <- sapply(df$KEGG_ko, get_pathway_names)

# Save results
output_file <- "KO_to_Pathway_Output_down.csv"
write.csv(df, output_file, row.names = FALSE)

# Print first few rows
print(head(df))

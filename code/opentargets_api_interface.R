# Install relevant library for HTTP requests
library(tidyjson)
library(httr)
library(tidyverse)

# Set gene_id variable
# dz_id <- "EFO_0000249" # alzheimer's
dz_id <- "MONDO_0004975" #late-onset AD
# dz_id <- "EFO_1001870" #late-onset AD

# API sandbox link:
# https://api.platform.opentargets.org/api/v4/graphql/browser?query=query%20targetInfo%20%7B%0A%20%20disease%28efoId%3A%20%22EFO_0000249%22%29%20%7B%0A%20%20%20%20id%0A%20%20%20%20name%0A%20%20%20%20associatedTargets%28%20page%3A%7Bsize%3A7213%2Cindex%3A0%7D%2C%20enableIndirect%3A%20true%29%20%7B%0A%20%20%20%20%20%20count%0A%20%20%20%20%20%20rows%20%7B%0A%20%20%20%20%20%20%20%20target%20%7B%0A%20%20%20%20%20%20%20%20%20%20id%0A%20%20%20%20%20%20%20%20%20%20approvedSymbol%0A%20%20%20%20%20%20%20%20%20%20approvedName%0A%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20score%0A%20%20%20%20%20%20%20%20datatypeScores%20%7B%0A%20%20%20%20%20%20%20%20%20%20id%0A%20%20%20%20%20%20%20%20%20%20score%0A%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%7D%0A%20%20%20%20%7D%0A%20%20%7D%0A%7D%0A

# Build query string
query_string = "
query disease($efoId: String!) {
  disease(efoId: $efoId) {
    id
    name
    associatedTargets( page:{size:8000,index:0}, enableIndirect: true) {
      count
      rows {
        target {
          id
          approvedSymbol
          approvedName
        }
        score
        datatypeScores {
          id
          score
        }
      }
    }
  }
}
"

# Set base URL of GraphQL API endpoint
base_url <- "https://api.platform.opentargets.org/api/v4/graphql"

# Set variables object of arguments to be passed to endpoint
variables <- list("efoId" = dz_id)

# Construct POST request body object with query string and variables
post_body <- list(query = query_string, variables = variables)

# Perform POST request
r <- POST(url=base_url, body=post_body, encode='json')

# gather data to tibble
open.targets <- content(r)$data[1] %>% 
  enter_object(associatedTargets) %>% 
  enter_object(rows) %>% gather_array() %>% spread_all() %>% 
  rename(tg.index=array.index, tg.score=score) %>% 
  enter_object(datatypeScores) %>% gather_array() %>% spread_all() %>%
  as_tibble() %>% 
  pivot_wider(id_cols = c(target.id, target.approvedSymbol, tg.score), names_from = id, values_from = score) 

rm(post_body, r, base_url, variables, query_string, dz_id)


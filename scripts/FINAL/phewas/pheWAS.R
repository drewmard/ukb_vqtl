require(httr)
require(jsonlite)
require(purrr)

url <- "http://genetics-api.opentargets.io/graphql"
GQL <- function(query, 
                ..., 
                .token = NULL,
                .variables = NULL, 
                .operationName = NULL, 
                .url = url){
  pbody <- list(query = query, variables = .variables, operationName = .operationName)
  if(is.null(.token)){
    res <- POST(.url, body = pbody, encode="json", ...)
  } else {
    auth_header <- paste("bearer", .token)
    res <- POST(.url, body = pbody, encode="json", add_headers(Authorization=auth_header), ...)
  }
  res <- content(res, as = "parsed", encoding = "UTF-8")
  if(!is.null(res$errors)){
    warning(toJSON(res$errors))
  }
  res$data
}

pheWAS <- function(snp) {
  s <- paste0(
    '{
  pheWAS(variantId: "',snp,'") {
    associations {
      study {
        traitReported
        traitCategory
      }
      beta
      pval
      nTotal
      nCases
    }
  }
}'
  )
  res <- flatten(GQL(s))$associations;
}

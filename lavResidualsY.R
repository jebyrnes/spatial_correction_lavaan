lavResidualsY <- function(object,
                          ynames = lavNames(object, "ov.nox"),
                          xnames = lavNames(object, "ov.x")){
  pred <- lavPredictY(object,
                      ynames = ynames,
                      xnames = xnames) |> 
    as.data.frame()
  
  d <- inspect(object, "data") |> 
    as.data.frame()
  
  r <- lapply(names(pred), function(x){
    pred[[x]] - d[[x]]
  })
  
  res <- do.call(cbind, r) |> 
    as.data.frame()
  names(res) <- names(pred)
  
  res
}

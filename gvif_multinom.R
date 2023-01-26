gvif.multinom <- function(model){
  (terms <- attr(model$terms, "term.labels"))
  if(any(grepl(":" , terms))){stop("Interactions are not supported.")}
  (y.levs <- model$lev)  
  (V.all <- vcov(model))
  (V.noIntercepts <- V.all[!grepl("\\(Intercept\\)$", rownames(V.all), perl = T), 
                           !grepl("\\(Intercept\\)$", colnames(V.all), perl = T)])
  (R <- cov2cor(V.noIntercepts))
  (gvif <- numeric(length = length(terms)))
  (DF <- integer(length = length(terms)))
  (SE.multiplier <- numeric(length = length(terms)))
  (names(gvif) <- names(DF) <- names(SE.multiplier) <- terms)
  
  #Now we start looping through the predictors:
  for(i in terms){
    if(is.null(model$xlevels[[i]])){(x.lev.regex <- "$")}else{
      (x.lev.regex <- paste0("(", paste(model$xlevels[[i]][2:length(model$xlevels[[i]])], collapse = "|"), ")"))}
    (RegexToMatch <- paste0("^(", paste(y.levs[2:length(y.levs)], collapse = "|") ,"):", i, x.lev.regex ))
    (indices <- grep(RegexToMatch, rownames(R), perl = T))
    #Below is the actual calculation:
    (gvif[i] <- det(R[indices, indices]) * det(R[-indices, -indices]) / det(R))
    (DF[i] <- length(indices))
    (SE.multiplier[i] <- gvif[i]^(1/(2*DF[i])))
  }
  #And then we output the results, ordered by degree of SE inflation:
  (result <- cbind(GVIF = gvif, DF = DF, `GVIF^(1/(2df))` = SE.multiplier))
  return(result[order(result[,"GVIF^(1/(2df))"], decreasing = T),])}
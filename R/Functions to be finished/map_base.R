map.base<-
  function(
    KUDest,
    var = TRUE,
    id = NULL)
    {

    require(rworldmap)
    require(rworldxtra)

    if(is.null(KUDest$Spatial.Objects)){
      stop("Can't find saved spatial objects.\nMake sure you've run the HRSummary() function with 'storepoly = TRUE'")
    }

    if(full){
      plotstack<-unlist(KUDest$Spatial.Objects)[grep("*_full", names(unlist(KUDest$Spatial.Objects)))]
      names(plotstack)<-unlist(lapply(strsplit(names(plotstack), "[.]"), `[[`, 1))
    }




  }

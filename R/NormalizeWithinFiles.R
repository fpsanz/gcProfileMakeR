#' NormalizeWithinFiles
#'
#'Reads and pre-processes xls and csv files, filling in blank fields
#' and grouping compounds with similar retention time and CAS numbers in common.
#'
#' @param path Complete directory path where the files are located.
#' @param type Either 'xls' or 'csv'. Input file type. Default: xls
#' @param thr Establishes the threshold of retention time to be considered in each file
#'  for each peak in order to unify compounds. Default 0.2
#' @param savefiles A boolean. If you want to save the resulting files on disk. Default: FALSE
#' @param filterLowQual Numeric. Filters out secondary hits with low 
#' quality comparing to the first hit in a particular retention time. 
#' Percentage of quality similarity of secondary hits with first hit for 
#' a given retention time. Values in range 0-1. filterLowQual = 0, means do not filter. 
#' Default: 0.99, means remove hits with lower quality than 99% 
#' of the first hit (with the highest quality for that retention time).
#' @param cas2rm A Vector or NULL. Character vector containing the CAS to be removed
#' from the analysis. Default: NULL means do not filter out.
#' @param minQuality. Numeric or NULL. Filter out hits with quality lower than minQuality value.
#' @return A normdata class object
#'
#' 
#' @examples
#'
#' ## 
#'
#' path <- paste(system.file(package = "gcProfileMakeR"), "/extdata", sep="")
#' out1 <- NormalizeWithinFiles(path = path)
#' 
#' ## or with custom parameters
#' path <- paste(system.file(package = "gcProfileMakeR"), "/extdata", sep="")
#' out1 <- NormalizeWithinFiles(path = path, type = 'xls', thr = 0.2, savefile = FALSE,
#'                           filterLowQual = 0.99)
#'
#'
#' @author  Fernando Pérez-Sanz \code{fernando.perez8@um.es}, Victoria Ruiz \code{victoria.ruiz@@upct.es}
#' @export

NormalizeWithinFiles <- function(path, type = "xls", thr = 0.2, savefiles = FALSE,
                                 filterLowQual=0.99, cas2rm=NULL, minQuality = NULL ){

    # Check filterLowQual
  if (!is.numeric(filterLowQual) | filterLowQual<0 | filterLowQual >1){
    stop("FilterLowQual is not between 0-1 or is not numeric")
  }
  # Check type
  if (type=="xls"){
    filelist <- list.files(path, pattern = "(.xls)$", full.names = TRUE)
  } else if (type == "csv"){
    filelist <- list.files(path, pattern = "(.tsv)$|(.csv)$")
  } else {stop("Wrong type. Valid types are xls or csv")}
  # Check thr
  if(!is.numeric(thr) | thr<0){
    stop("Threshold is not numeric or is negative number")
  }
  #check savefiles
  if(!is.logical(savefiles)){
    stop("Value must be TRUE or FALSE")
  }
  #check savefiles
  if(!is.null(cas2rm) & !is.vector(cas2rm)){
      stop("cas2rm should be a character vector with CAS id.")
  }
  #check minQuality
  if(!is.null(minQuality)){
    if( !is.numeric(minQuality) | minQuality<0 | minQuality>99){
     stop("minQuality is not a number or is out of range 0-99.")
    }
  }
  outdf <- data.frame()
  cas <- vector()
  hit <- vector()
  cat("\nPlease wait, this may take a while\n")
  #cat("Processing files ")
  pgr <- "a"
  for (file in filelist){
    if (type =="xls"){
      kk <- readxl::read_excel(file,  sheet = "LibRes", skip = 8, col_names = TRUE)
      df<-kk[,c(1,2,4,9,10,12)]
    } else if (type == "csv"){
      #kk <- read.csv2( file, header = TRUE, sep = ";", stringsAsFactors = FALSE)
      kk <- readr::read_delim(file,  delim = "\t", escape_double = FALSE, trim_ws = TRUE, col_names = TRUE)
      df <- kk
    }
     
    # progress points
    if(pgr=="a"){
        cat("\r","Processing files    "); pgr="b"} 
    else if(pgr=="b"){
        cat("\r","Processing files .  "); pgr="c"    
    } else if(pgr=="c"){
        cat("\r","Processing files .. "); pgr="d"
    } else if(pgr == "d"){ 
        cat("\r","Processing files ..."); pgr="a" 
            }
    ## 
    
    df <- as.data.frame(df) # si cargamos excel en lugar de csv el objeto creado es un tibble (dataframe extraño) con esto se pasa a dataframe normal
    colnames(df)<-c("Compound","RT","Area","Hitname","Quality","CAS")
    hit <- append(hit, df$Hitname)
    cas <- append(cas, df$CAS)
    # bucle para completar los huecos
    for ( i in seq( 1, dim( df )[1] ) ){
      if( is.na( df[i, 1] ) ){ df[i,1:3]<- df[i-1, 1:3] }
    }

    df <- df[,-4]
    RowsNoCAS <- which(is.na(df$CAS)) #buscar lineas sin CAS (las hay)
    if (length(RowsNoCAS!=0) ) {df <- df[-RowsNoCAS, ] } # eliminarlas

    df2 <- df
    
    # En este punto es donde se eliminan los CAS de la lista para limpiar (añadido para versión V221)
    if(is.vector(cas2rm) & length(cas2rm)>=1){
        rmcases <- data.frame(cas=unique(cas2rm))
        df2 <- df2 %>% dplyr::filter(!CAS %in% rmcases$cas)
    }
    #...............................#
    
    # Eliminar hits con valores de calidad por debajo de un threshold definido por minQuality
    if( ! is.null(minQuality)){
      df2 <- df2 %>% dplyr::filter(Quality >= minQuality)
    }
    
    df2$revisado <-0
    dfaux <- list()
    #Generar lista de dataframes por RT
    for (i in seq(1,nrow(df2) ) ){
      if (df2$revisado[i] != 1){
        dfaux[[as.character(i) ]] <- df2[df2$RT == df2$RT[i] , ]
        df2$revisado[df2$RT == df2$RT[i] ] <- 1
      }
    }
    
    # unificar dataframes con dif de RT <=0.2 y CAS común
    for (i in seq(1, length(dfaux) - 1)) {
      if (plyr::empty(dfaux[[i]]) == TRUE) {
        #print("empty")
      }
      else {
        for (j in seq((i + 1), length(dfaux))) {
          if (plyr::empty(dfaux[[j]]) == TRUE) {
            #print("empty")
          }
          else{
            #print(paste("j", j))
            kk <- which(dfaux[[j]]$CAS %in% dfaux[[i]]$CAS)
            if ((abs(dfaux[[j]]$RT[length(dfaux[[j]]$RT)] -  dfaux[[i]]$RT[length(dfaux[[i]]$RT)]) <= thr) &
                length(kk) != 0) {
              dfaux[[i]] <- rbind(dfaux[[i]], dfaux[[j]])
              dfaux[[j]] <- data.frame()
            }
          }
        }
      }
    }

    # # Eliminar compuestos "colegas" con quality < xx% de la mayor quality de ese RT
    # # y eliminar líneas idénticas
    if (filterLowQual != 0 ){
      for (i in seq(1,length(dfaux) ) ){
            if(!plyr::empty(dfaux[[i]])){
                maxQ = max(dfaux[[i]]$Quality)
                remove = which( dfaux[[i]]$Quality/maxQ < filterLowQual)
                if (length(remove)!=0){
                    dfaux[[i]] <- dfaux[[i]][-remove, ]
                }
                dfaux[[i]] <- unique(dfaux[[i]])
            }
      }
    }

    # Recalcular RT de cada grupo con la media y cada area con la suma,
    # si hay CAS iguales nos quedamos con el de mayores valores
    for (i in seq( 1, length(dfaux) ) ){
      if ( plyr::empty( dfaux[[i]] ) == TRUE ) {
        #print("empty")
      }
      else{
        dfaux[[i]]$RT <- mean( unique( dfaux[[i]]$RT ) )
        dfaux[[i]]$Area <- sum( unique( dfaux[[i]]$Area ) )
        dfaux[[i]] <- aggregate( .~CAS ,data=dfaux[[i]],  FUN = max )
      }
    }

    dfToExport <- plyr::ldply(dfaux, data.frame)
    dfToExport <- dfToExport[c("Compound","RT","Area","Quality","CAS")]

    if (type =="xls"){
      filecol <- sub(".xls","",file)
    } else if (type == "csv"){
      filecol <- sub(".csv","",file)
    }

    dfToExport$filecol <- filecol
    outdf <- rbind(outdf, dfToExport)

    # Guardar fichero en disco
    if (savefiles){
      fileout <- paste( filecol, "_step1.txt", sep="" )
      write.table(dfToExport, fileout, sep=";", row.names = FALSE, quote = FALSE)
    }


  }
  out <- new("normdata")
  hitfinal <- data.frame(cas=c(cas),hit=c(hit), stringsAsFactors = FALSE )
  hitfinal <- hitfinal[!duplicated(hitfinal[,c('cas')]),]
  out@data <- outdf
  out@hitname <- hitfinal
  return(out)
}


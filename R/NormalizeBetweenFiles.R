#' NormalizeBetweenFiles
#'
#' The function groups together compounds with similar retention time among all files,
#'  determining for each group which CAS will be the most representative for that group.
#'
#' @param data Object class normdata from NormalizeWithinFiles() function
#' @param thr Threshold of time to be considered between all files in order to unify compounds. Default 2
#' @param savefiles A boolean. If you want to save the resulting files on disk. Default: FALSE
#' @param filterLowQual Numeric. Filters out hits with quality below the established parameter.
#'  Values in range 0-100. Default: 0, means do not filter.
#'
#' @return A normdata class object
#'
#' @examples
#'
#' ## This example runs the function with default parameters. 
#' ## Where out1 is a normdata class object output from the NormalizeWithinFiles function.
#' path <- paste(system.file(package = "gcProfileMakeR"), "/extdata", sep="")
#' out1 <- NormalizeWithinFiles(path = path)
#' out2 <- NormalizeBetweenFiles(data = out1)
#' 
#' ## full param configuration
#' out2 <- NormalizeBetweenFiles(data = out1, thr = 2, savefiles = FALSE, filterLowQual = 0)
#'
#' @seealso \code{\link{NormalizeWithinFiles}}
#'
#' @author  Fernando Pérez-Sanz \code{fernando.perez8@um.es}, Victoria Ruiz \code{victoria.ruiz@upct.es}
#' @export

NormalizeBetweenFiles <- function(data, savefiles=FALSE, thr = 2, filterLowQual = 0){
    . <- NULL
    vals <- NULL
  # Check filterLowQual
  if (!is.numeric(filterLowQual) | filterLowQual<0 | filterLowQual >100){
    stop("FilterLowQual is not between 0-100 or is not numeric")
  }
  # Check thr
  if(!is.numeric(thr) | thr<0){
    stop("Threshold is not numeric or is negative number")
  }
  #check savefiles
  if(!is.logical(savefiles)){
    stop("Value must be TRUE or FALSE")
  }


  cat("\nPlease wait, this may take a while\n")
  #cat("Processing data ")

  #Cargar datos en dataframe
  df <- data@data

  #Ordenar dataframe por RT
  df <- df[with(df, order(RT)),]

  #Añadir columna revisado/no revisado
  df$revisado <- rep(0, nrow(df))

  #Eliminar registros (compuestos) con calidad <70, no debe haber xq en paso 1 se eliminaron colegas
  #con Quality < que el 95% del de mayor qual de su RT
  if (filterLowQual>0){
    df = df[df$Quality>=filterLowQual, ]
  }
  # convertir data.frame en lista tribble y despues en lista convencional
  ldf2 <- df %>% dplyr::group_by(df$RT, df$filecol) %>%
    dplyr::do(vals = data.frame(.)) %>%
    select(vals) %>% lapply(function(x) {(x)})

  ld <- as.list(ldf2)
  resultados <- list()
  j=1
  pgr="a"
    # Crear tablas intermedias de registros con algún CAS común y sus colegas de RT
  while(any(unlist(lapply(ld$vals, function(x) any(x$revisado == 0))) == TRUE )){
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
    
    dfaux <- data.frame()
    #print("nuevo dfaux")
    for (i in seq(1, length(ld$vals))){
      #print(i)
      # si dfaux esta vacio y en i esta revisado pasar al siguiente
      if(plyr::empty(dfaux) & any(ld$vals[[i]]$revisado %in% 1)){
        next
      }
      #si dfaux esta vacio y ld$vals en i esta sin revisar
      else if( plyr::empty(dfaux) & any(ld$vals[[i]]$revisado %in% 0) ){
        dfaux <- rbind(dfaux, ld$vals[[i]])
        ld$vals[[i]]$revisado <- 1
        #print("por arriba")
      } #si no, if dfaux no esta vacio y  ld$vals en i esta sin revisar
      #..y tienen dfaux tiene algun cas comun con ld$vals y la diferencia de Rt es menor que thr
      else if( !plyr::empty(dfaux)  & any(ld$vals[[i]]$revisado %in% 0) &
               any(ld$vals[[i]]$CAS %in% dfaux$CAS) & (ld$vals[[i]]$RT[1] - tail(dfaux$RT,1) <= thr) ) {
        dfaux <- rbind(dfaux, ld$vals[[i]])
        ld$vals[[i]]$revisado <- 1
      }
      #si no, si la diferencia en RT es mayor de thr sale del for
      else if( !plyr::empty(dfaux) & ( ( ld$vals[[i]]$RT[1] - tail(dfaux$RT,1) ) > thr) ){
        break()
      }
    }
    dfaux <- dfaux[,c("Compound","RT","Area","Quality","CAS","filecol","revisado")]
    resultados[[j]]<-dfaux
    ld$vals = ld$vals[ -which(unlist(lapply(ld$vals, function(x) all(x$revisado==1)))  )]
    j=j+1
  }


  ## pretablas generadas, asignar CAS comun

  resultados2 <- resultados
  for (i in seq(1, length(resultados2))) {
    kkmean <- aggregate(Quality ~ CAS, data = resultados2[[i]], mean)
    kklength <- aggregate(Quality ~ CAS, data = resultados2[[i]], length)
    kkmedian <- aggregate(Quality ~ CAS, data = resultados2[[i]], median)
    maxQual <- which(kkmean$Quality == max(kkmean$Quality)) # que elementos tienen la max Qual media
    if (length(maxQual) > 1) { # si hay mas de una qual maxima
      maxN <- which(kklength$Quality[maxQual] == max(kklength$Quality[maxQual])) # calcula cuantos valores tenía cada CAS
      if (length(maxN) == 1) {
        Norder <- maxQual[maxN]
      }
      else{
        maxQualMedian <- which(kkmedian$Quality[maxQual] == max(kkmedian$Quality[maxQual]))
        Norder <- maxQual[maxQualMedian[1]]
      }
    }
    else {
      Norder <- maxQual
    }
    CASasignado <- kkmean$CAS[Norder]

    resultados2[[i]]$CAS <- CASasignado

    #V3 al agrupar por area y quality me quedo con un solo registro por fichero con los valores más altos de las dos variable
    tmp <- aggregate(cbind(Quality,Area)~., data=resultados2[[i]][,-c(1,2)], max) #V3

    resultados2[[i]] <- data.frame()
    resultados2[[i]] <- tmp # V2
    resultados2[[i]] <- resultados2[[i]][,c(1,2,5,3,4)] # V2
  }
  dfout<-data.frame()
  #convertir lista en data.frame
  dfout <- plyr::ldply(resultados2 )

  dfout <- dfout[,c("Area","Quality","CAS","filecol")]
  colnames(dfout) <- c("area", "quality", "cas", "filecol")
  # Regenerar archivos desde las pretablas

  filelist <- unique(sort(df$filecol))

  if (savefiles){
    for (file in filelist){
      dfaux <- data.frame()
      for (i in seq(1, length(resultados2) ) ){
        wt <- which( resultados2[[i]]$filecol %in% file)
        dfaux <- rbind(dfaux, resultados2[[i]][wt,] )
      }

      write.table(dfaux, file = paste(file,"_step2.txt", sep = ""),
                  row.names = FALSE, col.names = TRUE, sep = "\t", quote =FALSE )}
  }
  out <- new("normdata")
  out@data <- dfout
  out@hitname <- data@hitname
  return(out)
}


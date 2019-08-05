#' getGroups
#'
#'A function to finally determine emission profiles and non-constitutive compounds emitted.
#'
#' @param data Object class normdata from NormalizeBetweenFiles() function
#' @param savefiles A boolean indicating whether to save the results in the working directory. Default FALSE
#' @param verbose A boolean indicating whether to show the profiles in console. Default FALSE
#' @param qcutoff A numerical integer value indicating the average quality cut-off point to determine which compounds should be
#'  listed as "Profile" and "Non Constitutive by Quality". Default value = 85
#' @param ncFreqCutoff A numerical float value between 0-1 indicating minimal frequency that establishes compounds
#' considered as non constitutive by frequency. This value always must be lower than pFreqCutoff.
#' Default value = 0.3
#' @param pFreqCutoff A numerical float value between 0-1 indicating minimal frequency that establishes compounds
#' considered as profile. This value always must be higher than ncFreqCutoff. Default value = 0.9
#' 
#'
#' @return A gcprofile class object
#'
#'
#' @examples
#' ## This example run the function with default parameters. 
#' ## Where out2 is a normdata class object output from the NormalizeBetweenFiles function.
#' path <- paste(system.file(package = "gcProfileMakeR"), "/extdata", sep="")
#' cat(path)
#' out1 <- NormalizeWithinFiles(path = path)
#' out2 <- NormalizeBetweenFiles(data = out1)
#' out <- getGroups(data = out2)
#'
#' ## full param configuration
#' out <- getGroups(data = out2, savefiles = FALSE, verbose = TRUE, 
#'                    qcutoff = c(85), ncFreqCutoff = c(0.3), pFreqCutoff = c(0.9) )
#'
#' @seealso \code{\link{NormalizeBetweenFiles}}
#'
#' @author  Fernando Pérez-Sanz \code{fernando.perez8@@um.es}, Victoria Ruiz \code{victoria.ruiz@@upct.es}


# gsub("^(0*)", "", "0000000064-76-2") eliminar los ceros delante de los CAS

getGroups <- function(data, savefiles = FALSE, verbose = TRUE, qcutoff = 85,
                       ncFreqCutoff = 0.3, pFreqCutoff = 0.9 ){


  # Check qcutoff
  if(!is.numeric(qcutoff) | qcutoff<0 | qcutoff>100){
    stop("qcutoff is not numeric or is out of range (0-100)")
  }
  #check savefiles
  if(!is.logical(savefiles)){
    stop("Value must be TRUE or FALSE")
  }
  #check ncFreqCutoff
  if(!is.numeric(ncFreqCutoff) | length(ncFreqCutoff) != 1 | ncFreqCutoff >= pFreqCutoff
     | TRUE %in% (ncFreqCutoff > 1) | TRUE %in% ncFreqCutoff < 0 ){
    stop("Value out of range. Value must be in the range 0 - 1 and must be lower than pFreqCutoff")
  }
  #check pFreqCutoff
  if(!is.numeric(pFreqCutoff) | length(pFreqCutoff) != 1 | pFreqCutoff <= ncFreqCutoff
    | TRUE %in% (pFreqCutoff > 1) | TRUE %in% pFreqCutoff < 0 ){
        stop("Values out of range. Values must be in the range 0 - 1 and pFreqCutoff must be higher than ncFreqCutoff")
    }
  #check verbose
  if(!is.logical(verbose)){
    stop("Value must be TRUE or FALSE")
  }

  # definir clase de tipo gcprofile
  #' @export
  out = new("gcprofile")
  out@data <- data@data
  # Esta parte corresponde al paso 3

  #ordenar fichero por CAS
  us3 <- data@data
  us3<-us3[with(us3, order(cas)),]

  #juntar bajo mismo CAS todas las áreas, todas las qualities y todos los file names
  us4A<-aggregate(area~cas, us3, function(x){ paste(x, sep = " ", collapse = ",")})
  us4Q<-aggregate(quality~cas, us3, function(x){ paste(x, sep = " ", collapse = ",")})
  us4F<-aggregate(filecol~cas, us3, function(x){ paste(x, sep = " ", collapse = ",")})

  ######### OK #### alternativa para no tener que escribir en disco y despues leer #####

  k1 <- as.data.frame(stringr::str_split_fixed(us4A$area, ",", n = Inf))
  k2 <- data.frame(cas = as.character(us4A$cas))
  k2$cas <- as.character(k2$cas)
  k1 <- as.data.frame(apply(k1, MARGIN = 2, function(x) as.numeric(as.character(x))))
  us5A <- cbind(k2,k1)

  k1 <- as.data.frame(stringr::str_split_fixed(us4Q$quality, ",", n = Inf))
  k2 <- data.frame(cas = as.character(us4Q$cas))
  k2$cas <- as.character(k2$cas)
  k1 <- as.data.frame(apply(k1, MARGIN = 2, function(x) as.numeric(as.character(x))))
  us5Q <- cbind(k2,k1)

  k1 <- as.data.frame(stringr::str_split_fixed(us4F$filecol, ",", n = Inf))
  k2 <- data.frame(cas = as.character(us4F$cas))
  k2$cas <- as.character(k2$cas)
  k1 <- as.data.frame(apply(k1, MARGIN = 2, function(x) as.character(x)))
  us5F <- cbind(k2,k1)

  tableaux <- rbind( us5A, us5Q, us5F )
  tableaux$type <- rep(c("areas","quals","Files"), each=dim(us5A)[1] )
  tableaux <- tableaux[ with(tableaux, order(cas)), ]

  filelist <- unique(data@data$filecol[order(data@data$filecol)] )
  if (dim(us5A)[2]-1 < length(filelist) ){
    colrest <- length(filelist) - ( dim(us5A)[2]-1 )
    df3aux <- matrix(data = 0, nrow = dim(us5A)[1], ncol = colrest)
    us5A <- cbind(us5A, df3aux)
    us5F <- cbind(us5F, df3aux)
    us5Q <- cbind(us5Q, df3aux)
  }


  #Sustituir Na por 0
  us5A[is.na(us5A)]<-0
  us5Q[is.na(us5Q)]<-0

  out@areas <- us5A
  out@quals <- us5Q
  out@AuxTable <- tableaux 

  #Escribir ficheros finales y borrar temporales
  if (savefiles){
    #write.table(us5A, file="areas.csv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep=",")
    #write.table(us5Q, file="qual.csv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep=",")
    write.table(tableaux, file = "AuxTable.csv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep=";")
  }
  rownames(us5Q) <- us5Q$cas
  us5Q <- us5Q[,-c(colnames(us5Q)=='cas')]
  #return(us5Q)

  # Esta parte corresponde al paso 4
  ###########################################################################################################
  qual <- us5Q
  qcff <- qcutoff
  minf <- ncFreqCutoff
  maxf <- pFreqCutoff

  numQuals <- apply(qual, MARGIN = 1, function(x)sum(x!=0) )
  maxQuals <- max(numQuals)
  qualmeans <- apply(qual, MARGIN = 1, function(x) mean(x[x!=0]))

  #Calcular el perfil: CAS con calidad media >qcutoff (80) y que aparecen >90% de las veces
  profileNumber <- which((numQuals/maxQuals) >= maxf & qualmeans >= qcff)
  profile <- rownames(qual)[profileNumber]

  #Calcular no constitutivo por frecuencia: CAS que tienen calidad media >80  y frecuencia <0.9
  newQual <- qual[-c(profileNumber),]
  newqualmeans <- apply(newQual, MARGIN = 1, function(x) mean(x[x!=0]))
  nonConstitutiveByFreqNumber <- which(newqualmeans >= qcff)
  nonConstitutiveByFreq <- rownames(newQual)[nonConstitutiveByFreqNumber]
  
  #Calcular no constitutivo por calidad: CAS que tienen frecuencia >0.3 y calidad <80
  toerase <- c(profile,nonConstitutiveByFreq)
  numQuals <- numQuals[-c(which(names(numQuals) %in% toerase) ) ]
  nonConstitutiveByQualNumber <- which((numQuals/maxQuals) >= minf)
  nonConstitutiveByQual <-names(nonConstitutiveByQualNumber)
  
  #Asignar CAS a  perfiles
  hitname <- data@hitname
  profileFinal <- hitname[which(hitname$cas %in% profile),]
  nonConstitutiveQualFinal <- hitname[which(hitname$cas %in% nonConstitutiveByQual),]
  nonConstitutiveFreqFinal <- hitname[which(hitname$cas %in% nonConstitutiveByFreq),]
  
  #Generar objeto de salida
  out@profile <- profileFinal
  out@nonConstitutiveFreq <- nonConstitutiveFreqFinal
  out@nonConstitutiveQual <- nonConstitutiveQualFinal

  if(savefiles){
    write.table(profileFinal, file="profile.csv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep=";")
    write.table(nonConstitutiveFreqFinal, file="nonConstitutiveFreq.csv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep=";")
    write.table(nonConstitutiveQualFinal, file="nonConstitutiveQual.csv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep=";")
  }


  #Imprimir resultados por pantalla
  if (verbose){
    cat("Profile:\n")
    cat(paste("Compounds with quality >",qcff,"and frequency >",maxf,"\n"))
    cat(profile, sep = " ; ", fill = TRUE)
    cat("\n")
    cat("Non Constitutive by quality:\n")
    cat(paste("Compounds with quality <",qcff,"and frequency >",minf,"\n"))
    cat(nonConstitutiveByQual, sep =" ; ", fill = TRUE)
    cat("\n")
    cat("Non Constitutive by frequency:\n")
    cat(paste("Compounds with quality >",qcff,"and frequency <",maxf,"\n"))
    cat(nonConstitutiveByFreq, sep = " ; ", fill = TRUE)
  }
  if(savefiles){
  cat("\n\nFiles 'profile.csv', 'nonConstitutiveFreq.csv', 'nonConstitutiveQual.csv'\nhave been saved in the working directory.\n")
  cat("\nFiles 'areas.csv', 'qual.csv', 'AuxTable.csv'\nhave been saved in the working directory.\n")
    }
  return(out)
}

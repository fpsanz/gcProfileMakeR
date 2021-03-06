#' plotGroup
#'
#' The function plots both line chart with average qualities and bars chart with average areas of
#'  compounds grouped as “profile”, “non constititive by frequency” or “non constitutive by quality”
#'
#' @param data Output object from getGroups() function
#' @param compoundType Set of compounds to be plotted. Valid values are: "profile", "nonconstitutivebyfreq"
#' and "nonconstitutivebyqual". Its abbreviations "p", "ncf", "ncq" are also valid. Default value is "profile"
#' @param maxBars A numerical integer value indicating the maximum number of bars to be displayed. By default displays all bars
#' 
#' @return ggplot object
#'
#' @examples
#'
#' ## Where out is a gcprofile class object output from the getGroups function.
#' path <- paste(system.file(package = "gcProfileMakeR"), "/extdata", sep="")
#' out1 <- NormalizeWithinFiles(path = path)
#' out2 <- NormalizeBetweenFiles(data = out1)
#' out <- getGroups(data = out2)
#' p <- plotGroup(data = out, compoundType = "p")
#' p
#' ## or directly
#' plotGroup(data = out, compoundType = "p")
#' 
#'
#' @author  Fernando Pérez-Sanz \code{fernando.perez8@@um.es}, Victoria Ruiz \code{victoria.ruiz@@upct.es}
#' @export

plotGroup <- function(data, compoundType = "profile", maxBars = NA ){
    require(ggplot2)
    require(egg)
    require(tidyverse)
    require(plotly)
    quals <- NULL
    CAS <- NULL
    todas_areas <- data@areas
    todas_qual <- data@quals
    titulo <- ""
    if(compoundType=="profile" | compoundType == "p"){
        outareas <- data@profile
        titulo <- "Profile compounds"
    }
    else if(compoundType=="nonconstitutivebyfreq" | compoundType == "ncf"){
        outareas <- data@nonConstitutiveFreq   
        titulo <- "Non constitutive compounds by frequency"
    }
    else if(compoundType=="nonconstitutivebyqual" | compoundType == "ncq"){
        outareas <- data@nonConstitutiveQual
        titulo <- "Non constitutive compounds by quality"
    }
    else{
        stop("Not valid option. Please see help")
        
    }
    if(is.na(maxBars)){
        maxBars <- dim(data@areas)[1]
    }
    #seleccionar datos
    areas <- todas_areas[which(todas_areas[,1] %in% outareas$cas), ]
    medias <- round(apply(areas[,2:ncol(areas)], MARGIN =1 , function(x) mean(x[x!=0]) ),1)
    sd <- round(apply(areas[,2:ncol(areas)], MARGIN =1 , function(x) sd(x[x!=0]) ),1)
    n <- apply(areas[,2:ncol(areas)], MARGIN =1 , function(x) length(x[x!=0]) )
    se <- sd/sqrt(n)
    
    # quals del perfil
    qual <- todas_qual[which(todas_qual[,1] %in% outareas$cas), ]
    mediasq <- round(apply(qual[,2:ncol(qual)], MARGIN =1 , function(x) mean(x[x!=0]) ),1)
    sdq <- round(apply(qual[,2:ncol(qual)], MARGIN =1 , function(x) sd(x[x!=0]) ),1)
    nq <- apply(qual[,2:ncol(qual)], MARGIN =1 , function(x) length(x[x!=0]) )
    seq <- sdq/sqrt(nq)
    
    # preparar datos
    mydf <- data.frame( cbind(areas = medias,quals = mediasq, stderr = se, stderrq = seq),
                        stringsAsFactors = FALSE)
    mydf <- cbind(CAS=areas$cas, mydf)
    mydf$CAS <- as.character(mydf$CAS)
    mydf <- mydf[with(mydf, order(-areas)),]
    # Si los datos son > 30 y largeData = FALSE me quedo con los 50 primeros
    if(dim(mydf)[1] > maxBars){
        mydf <- mydf[1:maxBars,]
    }
    
    tmp <- left_join(mydf, outareas, by = c("CAS"="cas"))
    tmp <- tmp %>% mutate(lowerErr = areas-stderr,
                          upperErr = areas+stderr,
                          lowerErrQ = quals-stderrq,
                          upperErrQ = quals+stderrq )
    tmp$hit <- paste0("* ",tmp$hit)
    
    ## definir grafico
    dodge <- position_dodge(width = 0.9)

    p1 <- ggplot(tmp, aes(x = CAS, y = areas ) ) +
        geom_bar(stat = "identity", aes(fill=areas), show.legend = F) +
        labs(y = "Average area", x = "Compounds") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        scale_x_discrete(limits = tmp$CAS) +
        theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), units = "points"))+
        geom_errorbar(aes(ymin=lowerErr, ymax=upperErr), position = dodge, width = 0.25, na.rm=TRUE, colour="red")
    
    p2 <- ggplot(tmp, aes(x = CAS, y = quals, group = 1)) +
        geom_line(col = "red") +
        #geom_point(colour="red") +
        labs(y = "Average quality") +
        theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), units = "points"),
              axis.title.y = element_text(vjust = 0.75),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) +
        ggtitle(label = titulo) +
        theme(plot.title = element_text(hjust = 0.5, size = 15)) +
        scale_x_discrete(limits = tmp$CAS)+
        geom_errorbar(aes(ymin=lowerErrQ, ymax=upperErrQ), position = dodge, width = 0.25, na.rm=TRUE)
    
    #ggarrange(p2,p1, heights = c(1,5))
    p11 <- ggplotly(p1)    
    p22 <- ggplotly(p2)
    
    for( i in c( 0:( length( tmp$CAS )-1 ) ) ){
        p11$x$data[[( length( tmp$CAS)-i ) ]]$text <- gsub( "\\$\\$", "\n*", tmp$hit[(i+1)] ) 
    } 
    subplot(p22, p11, nrows = 2, heights = c(0.25,0.75), which_layout = 1) 
 
}


#' Class gcprofile
#'
#' This is a class representation grouping different outputs from \code{getGroups} function
#'
#' @slot data Data frame from \code{NormalizeBetweenFiles} function
#' @slot profile Data frame with compounds classified as "Profile".
#'  List of CAS numbers associated to a list of proposed names.
#' @slot nonConstitutiveFreq Data frame with compounds classified as "Non Constitutive by Frequency".
#'  List of CAS numbers associated to a list of proposed names.
#' @slot nonConstitutiveQual Data frame with compounds classified as "Non Constitutive by Quality".
#'  List of CAS numbers associated to a list of proposed names.
#' @slot areas Data frame with CAS and areas of each CAS
#' @slot quals Data frame with CAS and qualities of each CAS
#' @slot AuxTable Data frame with CAS numbers associated to their areas and qualities in each file.
#'  Data presented is the result of the consensus reached by NormalizeWithinFiles and NormalizeBetweenFiles.
#'
#' @seealso Examples in \code{\link{getGroups}}
#'
#' @author Fernando Pérez-Sanz \code{fernando.perez8@um.es}, Victoria Ruiz \code{victoria.ruiz@upct.es1}
#'
#' @import readxl stringr readxl dplyr ggplot2 egg tidyr
#' @importFrom plyr empty
#' @importFrom methods new
#' @importFrom stats aggregate median
#' @importFrom utils read.csv2 tail write.table
#' 
#' @export

setClass("gcprofile", slots= c(data = "data.frame", profile = "data.frame", nonConstitutiveFreq = "data.frame",
                               nonConstitutiveQual = "data.frame", areas = "data.frame", quals = "data.frame",
                               AuxTable = "data.frame"))


setClass("normdata", slot = c(data = "data.frame", hitname = "data.frame"))

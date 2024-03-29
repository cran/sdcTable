% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\docType{class}
\name{safeObj-class}
\alias{safeObj-class}
\title{S4 class describing a safeObj-object}
\description{
Objects of class \code{safeObj} are the final result after protection a
tabular structure. After a successful run of \code{\link{protectTable}}
an object of this class is generated and returned. Objects of class
\code{safeObj} contain a final, complete data set (slot \code{finalData})
that has a column showing the anonymization state of each cell and the
complete information on the dimensional variables that have defined the table
that has been protected (slot \code{dimInfo}). Also, the number of
non-duplicated table cells (slot \code{nrNonDuplicatedCells}) is returned
along with the number of primary (slot \code{nrPrimSupps}) and secondary
(slot \code{nrSecondSupps}) suppressions. Furthermore, the number of cells
that can be published (slot \code{nrPublishableCells}) and the algorithm that
has  been used to protect the data (slot \code{suppMethod}) is returned.
}
\details{
\describe{
\item{slot \code{finalData}:}{a data.frame (or NULL) featuring columns for each variable defining the table (with their original codes), the cell counts and values of any numerical variables and the anonymization status for each cell with
\itemize{
\item \code{s, z}: cell can be published
\item \code{u}: cell is a primary sensitive cell
\item \code{x}: cell was selected as a secondary suppression}
}
\item{slot \code{dimInfo}:}{an object of class \code{\link{dimInfo-class}} holding all information on variables defining the table }
\item{slot \code{nrNonDuplicatedCells}:}{numeric vector of length 1 (or NULL) showing the number of non-duplicated table cells. This value is different from 0 if any dimensional variable features duplicated codes. These codes have been re-added to the final dataset.}
\item{slot \code{nrPrimSupps}:}{numeric vector of length 1 (or NULL) showing the number of primary suppressed cells}
\item{slot \code{nrSecondSupps}:}{numeric vector of length 1 (or NULL) showing the number of secondary suppressions}
\item{slot \code{nrPublishableCells}:}{numeric vector of length 1 (or NULL) showing the number of cells that may be published}
\item{slot \code{suppMethod}:}{character vector of length 1 holding information on the protection method}
}
}
\note{
objects of class \code{safeObj} are returned after the function \code{\link{protectTable}} has finished.
}
\author{
Bernhard Meindl \email{bernhard.meindl@statistik.gv.at}
}

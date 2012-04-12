\name{LumiBatch-class}
\docType{class}
\alias{class:LumiBatch}
\alias{LumiBatch-class}
\alias{LumiBatch}

\alias{se.exprs,LumiBatch-method}
\alias{se.exprs<-,LumiBatch,matrix-method}

\alias{beadNum}
\alias{beadNum<-}
\alias{beadNum,LumiBatch-method}
\alias{beadNum<-,LumiBatch,matrix-method}

\alias{detection}
\alias{detection<-}
\alias{detection,LumiBatch-method}
\alias{detection<-,LumiBatch,matrix-method}

\alias{initialize,LumiBatch-method}
\alias{summary,LumiBatch-method}
\alias{show,LumiBatch-method}
\alias{combine,LumiBatch,LumiBatch-method}
\alias{[,LumiBatch-method}
\alias{getHistory,LumiBatch-method}
\alias{getHistory}

\title{Class LumiBatch: contain and describe Illumina microarray data}
\description{This is a class representation for Illumina microarray data. It extends \code{\link[Biobase]{ExpressionSet}}.}
\section{Extends}{
   Directly extends class \code{\link[Biobase]{ExpressionSet}}.
}

\section{Creating Objects}{

  \code{new('LumiBatch',
    exprs = [matrix],
	se.exprs = [matrix],
	beadNum = [matrix],
	detection = [matrix],
    phenoData = [AnnotatedDataFrame],
    history = [data.frame],
    ...)
  }

  \code{LumiBatch} instances are usually created through
  \code{new("LumiBatch", ...)}. The arguments to \code{new}
  should include \code{exprs} and \code{se.exprs}, others can be 
  missing, in which case they are assigned default values.

  Objects can be created using the function \code{\link{lumiR}}.
}
\section{Slots}{
   Slot specific to \code{LumiBatch}:
   \describe{
   		\item{\code{history}:}{ a data.frame recording the operation history of the LumiBatch object.}
   		\item{\code{controlData}:}{ a data.frame with first two columns as "controlType" and "ProbeID". The rest columns are the control probe expression amplitudes for individual samples.}
   		\item{\code{QC}:}{ a the quality control information of the LumiBatch object, returned by lumiQ function.}
   }

   Slots inherited from \code{ExpressionSet}:
   \describe{
      \item{\code{assayData}}{ contains equal dimensional matrices:
 		\code{\bf exprs} (contains gene expression level, which is the mean of its bead replicates.),
 		\code{\bf se.exprs} (contains gene expression standard error, which is the standard error of its bead replicates.), 
 		\code{\bf beadNum} (records the number of beads for the probe.),
 		\code{\bf detection} (records the detection p-value of the probe. The number is from [0,1]. By default, < 0.01 indicates good detection.). 
		For more details of \code{assayData}, please see \code{\link[Biobase]{ExpressionSet}} } 
      \item{\code{phenoData}:}{See \code{\link[Biobase]{eSet}}}
      \item{\code{experimentData}:}{See \code{\link[Biobase]{eSet}}}
      \item{\code{annotation}:}{See \code{\link[Biobase]{eSet}}}
  	}
}

\section{Methods}{

  \bold{Class-specific methods:}
  \describe{
     \item{\code{se.exprs(LumiBatch)}, \code{se.exprs(LumiBatch,matrix)<-}:}{Access and
       set elements named \code{se.exprs} in the \code{AssayData-class} slot.}
     \item{\code{beadNum(LumiBatch)}, \code{beadNum(LumiBatch,matrix)<-}:}{Access and
       set elements named \code{beadNum} in the \code{AssayData-class} slot.}
     \item{\code{detection(LumiBatch)}, \code{detection(LumiBatch,matrix)<-}:}{Access and
       set elements named \code{detection} in the \code{AssayData-class} slot.}
     \item{\code{getHistory(LumiBatch)}:}{Access the operation history of \code{LumiBatch} object.}
  }

   \bold{Derived from \code{\link[Biobase]{ExpressionSet}}} (For the directly inherited methods, 
	 please see \code{\link[Biobase]{ExpressionSet}} and \code{\link[Biobase]{eSet}}):
   \describe{
     \item{\code{combine(LumiBatch,missing)}:}{Combine two LumiBatch objects, including \code{history} slot. See \code{\link[Biobase]{eSet}}}
	 \item{\code{exprs(LumiBatch)}, \code{exprs(LumiBatch,matrix)<-}:}{Access and
	    set elements named \code{exprs} in the \code{AssayData-class} slot.}
	 \item{\code{object[(i,j)}:}{Conduct subsetting of the data in a LumiBatch object }
   }

   \bold{Standard generic methods} (For the directly inherited methods, 
	 please see \code{\link[Biobase]{ExpressionSet}} and \code{\link[Biobase]{eSet}}):
  \describe{
    \item{\code{initialize(LumiBatch)}:}{Object instantiation, used
      by \code{new}; not to be called directly by the user.}
    \item{\code{validObject(LumiBatch)}:}{Validity-checking method, ensuring that \code{exprs} and \code{se.exprs} 
	  is a member of \code{assayData}. Other validity check is the same as \code{checkValidity(ExpressionSet)}.}
    \item{\code{show(LumiBatch)}}{A summary of the \code{LumiBatch} object.}
   }
}

\author{ Pan Du, Simon Lin }
\seealso{ \code{\link{lumiR}}, \code{\link{lumiT}}, \code{\link{lumiN}},
	\code{\link{boxplot-methods}}, \code{\link{pairs-methods}}, \code{\link{MAplot-methods}}
}
\examples{
## load example data
data(example.lumi)

## show the summary of the data
# summary(example.lumi)
example.lumi	

## get express matrix
temp <- exprs(example.lumi)

## get a subset
temp <- example.lumi[,1]	## retrieve the first sample

## get the probe id
featureNames(example.lumi)[1:3]

## combine LumiBatch objects
temp <- combine(example.lumi[,1], example.lumi[,3])
temp

}

\keyword{classes}

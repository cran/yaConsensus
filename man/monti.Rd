\name{monti}
\alias{monti}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute and display the Monti's statistics for class discovery.
%%  ~~function to do ... ~~
}
\description{The function computes Monti's statistics, and/or displays the correponding plots.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
monti(obj, gMax = nclass.Sturges(obj$labels), just_compute = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{An object of 'yaConsensus' class
%%     ~~Describe \code{obj} here~~
}
  \item{gMax}{an integer value indicating the maximum number cluster to be explored
%%     ~~Describe \code{gMax} here~~
}
  \item{just_compute}{A logical value indicating if Monti's statistics have to be just computed (TRUE) or, in addition, displayed (FALSE, default)
%%     ~~Describe \code{just_compute} here~~
}
}
\details{
If the 'fname' slot of the input object is instantiated, the input object is updated with the Monti's statistics and saved.%%  ~~ If necessary, more details than the description above ~~
}
\value{
The same input object of 'yaConsensus' class with a named list in the new 'monti' slot
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{monti}{
A named list with the following slots:}
  \item{x}{a sequence of 500 knots from 0 to 1}
  \item{y}{a real matrix of 500 rows and gMax - 1 colums. Each columns stores the values of the empirical distribution function corresponing to the number of clusters from 2 to gMax. }
  \item{area}{a list of real values, each of them corresponing to the area under the empirical distribution function (as stored in y)}
%% ...
}
\references{
Monti et al. (2003) - Consensus Clustering: A Resampling-Based Method for Class Discovery and Visualization of Gene Expression Microarray Data - Machine Learning 52(1-2):91-118 <DOI: 10.1023/A:1023949509487>
}
\author{
Stefano M. Pagnotta
}
\note{
In case the 'monti' slot is instantiated, the function provides the graphical result.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{yaConsensus}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Generate data and annotation
n <- 50; m <- 3000
ddata <- matrix(rnorm(n * m), ncol = m)  
ddata[1:20, ] <- ddata[1:20, ] + 0.2
ddata[21:35, ] <- ddata[21:35, ] + 0.4
row.names(ddata) <- c(paste0("A", 1:20), paste0("B", 1:15), paste0("C", 1:15))
ddist <- dist(ddata)

annotation <- data.frame(row.names = rownames(ddata), clust = substr(rownames(ddata), 1, 1))
annotation.colorCode <- c("red", "blue", "green")
names(annotation.colorCode) <- c("A", "B", "C")

####### run in sequential mode
####### sampling the samples ....
aConsensus <- yaConsensus(ddist)
ans <- plot(aConsensus, G = 3, 
            annotation = annotation, 
            annotation.colorCode = annotation.colorCode)
monti(aConsensus)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

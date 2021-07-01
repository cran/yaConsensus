\name{yaConsensus}
\alias{yaConsensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{yaConsensus computes a number of hierarchical clusterings by sampling either samples or features.
}
\description{
This function mainly generates a list of "hclust" objects for downstream analysis.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
yaConsensus(ddata, runs = 1000, epsilon = 0.65, is_by_sample = TRUE, 
            distMethod = "euclidean", hcMethod = "ward.D2", prefix = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ddata}{either a data matrix (samples in rows, and features in columns), or a "dist" object.
%%     ~~Describe \code{ddata} here~~
}
  \item{runs}{an integer value for the number of samplings.
%%     ~~Describe \code{runs} here~~
}
  \item{epsilon}{a real value indicating the sampling rate.
%%     ~~Describe \code{epsilon} here~~
}
  \item{is_by_sample}{a logical value indicating if the sampling is by samples (TRUE) or features (FALSE).
%%     ~~Describe \code{is_by_sample} here~~
}
  \item{distMethod}{a character indicating the kind of distance for the inner clustering. It can be any of the methods from the \code{\link{dist}} function.  
%%     ~~Describe \code{distMethod} here~~
}
  \item{hcMethod}{a character indicating the linkage mathod of the inner clustering. It can be any of the methods from the \code{\link{hclust}} function.
%%     ~~Describe \code{hcMethod} here~~
}
  \item{prefix}{string specifying a prefix to store the results in a .RData file.  
%%     ~~Describe \code{prefix} here~~
}
}
\details{
This function can run sequentially or in parallel. In this case, it is necessary to register a cluster of CPUs according to doParallel protocol.  

To get the consensus clustering, the output of the function has to be processed with the plot() function. The consensus dissimilarity follows from the algorithm of Monti et al. (2003). The consensus clustering is from a hierarchical procedure (hclust) with "complete" linkage (outer hc method). 

%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
A named list with the following slots: 
\item{distMethod}{matches the input}
\item{hcMethod}{matches the input}
\item{lables}{a string list with the names of the samples}
\item{bySample}{matches 'is_by_sample' input parameter}
\item{epsilon }{matches the input}
\item{subsetDimension}{actual dimension of the subsets}
\item{runs}{matches the input}
\item{hclust}{a list of 'hclust' objects}
\item{elapsed_time}{time (in seconds) required}
\item{ncores}{the number of cores used}
%% ...
}
\references{
Monti et al. (2003) - Consensus Clustering: A Resampling-Based Method for Class Discovery and Visualization of Gene Expression Microarray Data - Machine Learning 52(1-2):91-118 <DOI: 10.1023/A:1023949509487>
}
\author{
Stefano M. Pagnotta
}
\note{
%%  ~~further notes~~

The plot function in the example provides an invisible result with detail ans statistics of the experiment.  
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{dist}}, \code{\link{hclust}}, \code{\link{doParallel}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Generate data and annotation
n <- 50; m <- 3000
ddata <- matrix(rnorm(n * m), ncol = m)  
ddata[1:20, ] <- ddata[1:20, ] + 0.2
row.names(ddata) <- c(paste0("A", 1:20), paste0("B", 1:30))
ddist <- dist(ddata)

annotation <- data.frame(row.names = rownames(ddata), clust = substr(rownames(ddata), 1, 1))
annotation.colorCode <- c("red", "blue")
names(annotation.colorCode) <- c("A", "B")

####### run in sequential mode
####### sampling the samples ....
(aConsensus <- yaConsensus(ddist))
plot(aConsensus, G = 2)

ans <- plot(aConsensus, G = 2, 
            annotation = annotation, 
            annotation.colorCode = annotation.colorCode)
summary(ans)
summary(ans, given = "clust")

####### sampling the features ....
(aConsensus <- yaConsensus(ddata, runs= 20, epsilon = 0.2, is_by_sample = FALSE))
ans <- plot(aConsensus, G = 2, 
            annotation = annotation, 
            annotation.colorCode = annotation.colorCode,
            matching_clustering = "clust")

summary(ans, given = "clust")


####### run in parallel mode
## uncomment to run

# require(doParallel)
# cpu_cluster <- makeCluster(3)
# registerDoParallel(cpu_cluster)

# (aConsensus <- yaConsensus(ddist))
# plot(aConsensus, G = 2)

#stopCluster(cpu_cluster)





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

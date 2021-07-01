\name{consensus.diss}
\alias{consensus.diss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes the consensus dissimilarity matrix.
%%  ~~function to do ... ~~
}
\description{
Computes the consensus dissimilarity according to the algorithm of Monti et al. (2003).%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
consensus.diss(cclusters, similarity = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cclusters}{a matrix of integers where the column are the samples, and the rows are different clusterings of the samples.
%%     ~~Describe \code{cclusters} here~~
}
  \item{similarity}{a logical value signaling if the similarity matrix is required. 
%%     ~~Describe \code{similarity} here~~
}
}
\details{
In any row of the ccluster matrix, the value 0 means that the corresponding sample is not assigned to any cluster. In this case, the dissimilarity is computed consistently.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
An object of the 'dist' class.
}
\references{
Monti et al. (2003) - Consensus Clustering: A Resampling-Based Method for Class Discovery and Visualization of Gene Expression Microarray Data - Machine Learning 52(1-2):91-118 <DOI: 10.1023/A:1023949509487>
}
\author{
Stefano M. Pagnotta
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{dist}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
clusters <- rep(1:3, c(3, 9, 18))
clusterings <- matrix(NA, ncol = 30, nrow = 50)
for(k in 1:50) clusterings[k,] <- sample(clusters)
ddist <- consensus.diss(clusterings)
class(ddist)
attr(ddist, "method")
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

\name{plot.yaConsensus}
\alias{plot.yaConsensus}
%%\alias{plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot the consensus dissimilarity matrix and provide the consensus clustering.
%%  ~~function to do ... ~~
}
\description{
This function processes the ouput of yaConsensus and acts as a wrapper to the pheatmap function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{plot}{yaConsensus}(x, G = 2, 
                            annotation = NULL, annotation.colorCode = NULL, 
                            matching_clustering = NULL, consensus_colors = NULL, 
                            reduce_to = -1, main = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object coming from yaConsensus().
%%     ~~Describe \code{obj} here~~
}
  \item{G}{an integer value indicating the number of clusters required for the consensus clustering. Default is 2.
%%     ~~Describe \code{G} here~~
}
  \item{annotation}{a data frame where the variables are annotations (as labels) of samples. The row-names have to match the names of the samples.
%%     ~~Describe \code{G} here~~
}
  \item{annotation.colorCode}{a string named list of color names. The names have to be values stored in the annotation data-frame.
%%     ~~Describe \code{G} here~~
}
  \item{matching_clustering}{a string value matching one of the annotation valiables in the annotation data-frame. The function tries to match at best the color coding of the selected variable in the data-frame.
%%     ~~Describe \code{G} here~~
}
  \item{consensus_colors}{a list of color provided to annotate the consensus clustering. If provided, the matching_clustering is overrided.
%%     ~~Describe \code{G} here~~
}
  \item{reduce_to}{an integer value; default = -1 means all samples are included in the graphical representation; a value greater than 10 forces the graphical output to show 'reduce_to' number of brahch of the consensusu dendrogram. See also detail.
%%     ~~Describe \code{G} here~~
}
  \item{main}{a main title for the plot; default is NULL (no title).
}
  \item{\dots}{parameters compatible with pheatmap function. 
%%     ~~Describe \code{\dots} here~~
}
}
\details{In the slot 'statistics', the function returns the same statistics of yaConsensus().

If the consensus analysis concerns single-cell, the'reduce_to' parameter has to be set to 500 (suggestion), given that the number of cells is often thousands. In this case, the original annotation is compressed to the number of items displayed.  If the consensus color annotation has to match an existing a-priori coloring code, the match is computed concerning the unreduced annotation. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
A named list with the following slots: 
\item{annotation}{a data frame. It is the same given in input, with 'consensus' and 'consensus.col' more variables.}
\item{ann_colors}{a named list of colors associated with each variable in the annotation data-frame.}
\item{hclust}{an object of hclust clust. It's the result of the hclust() applied to the consensus dissimilarity with the complete linkage.}
\item{statistics}{see Note}
}
\references{
Risso and Pagnotta (2021) - Per-sample standardization and asymmetric winsorization lead to accurate clustering of RNA-seq expression profiles - Bioinformatics, btab091, <DOI: 10.1093/bioinformatics/btab091>
%% ~put references to the literature/web site here ~
}
\author{
Stefano M. Pagnotta
}
\note{In case an 'annotation' is provided, with summary(), additional statistics are provided. 

The variable specified as given in summary(..., given = "some clustering") is assumed as a theoretical clustering, while the consensus clustering is the empirical one. The entropy accuracy, precision, and average (Risso and Pagnotta, 2021) are computed with the former assumptions. summary() itself returns all the statistics. "some clustering" has to be one of the column of the annotation data.frame. 
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[pheatmap:pheatmap]{pheatmap::pheatmap()}}, \code{\link{yaConsensus}}
}
\examples{
# see the examples in yaConsensus help.
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

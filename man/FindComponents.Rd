\name{FindComponents}
\alias{FindComponents}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find all the components in the disconnected graph.
}

\description{
To generate the graph with each cluster is connected.
}

\usage{
FindComponents(adj, label)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adj}{
    `[matrix]`: The adjacency matrix of the data.
}
  \item{label}{
    `[numeric]`: The label of the data.
}

}
\details{
The function is used to generate the connected graph by implementing the HCV algorithm to inspect the connectedness of each cluster, if the cluster is disconnected, then we split the cluster into several sub-clusters simplicity.
}
\value{
  \item{label }{The label of all the components}
  \item{neighbors }{`[list]`: The data in same clusters are placed in the same list components}
%% ...
}

\author{
DongDong-Zoez <lbry5230100@gmail.com> University of Taiwan NSYSU.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  `HCV`
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

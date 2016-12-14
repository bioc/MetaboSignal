\name{MS_interactionType}
\alias{MS_interactionType}
\title{Get interaction subtype}
\description{
This function generates a 3 column-matrix indicating the interaction subtype (e.g.
activation, binding) between node pairs involved in the signaling pathways of interest.
Each row of the matrix represents an edge between 2 nodes. The first column contains
the in-nodes,the second column contains the out-nodes, and the third column contains
the interaction subtype.
}
\usage{
MS_interactionType(paths, expand_genes = FALSE)
}
\arguments{
  \item{paths}{character vector containing the KEGG IDs for the
  signaling pathways of interest (organism-specific). For example, the KEGG ID
  for the pathway "insulin signaling pathway" in the rat is "rno04910".
  }
  \item{expand_genes}{logical scalar indicating whether the gene nodes will
  represent orthology IDs (FALSE) or organim-specific gene IDs (TRUE). By default,
  expand_genes = FALSE.
  }
}
\value{
A three-column matrix where each row represents an edge between two nodes.
}

\references{
Zhang, J.D. & Wiemann, S. (2009). KEGGgraph: a graph approach to KEGG PATHWAY in
R and Bioconductor. Bioinformatics, 25, 1470-1471.

http://www.kegg.jp/kegg/docs/keggapi.html
}
\examples{\donttest{
# Get interaction-type matrix of nodes involved in "rno04910" and "rno04151"
interaction_type  <- MS_interactionType (paths = c("rno04910", "rno04151"))
}
}




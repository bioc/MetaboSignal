\name{MS_GetShortestpaths}
\alias{MS_GetShortestpaths}
\title{Calculate shortest paths}
\description{
This function calculates the shortest path(s) between any two reachable nodes of
a network-table.
}
\usage{
MS_GetShortestpaths(network_table, source_node, target_node, mode = "out",
                    type = "first", network = FALSE, interaction_type = NULL)
}
\arguments{
  \item{network_table}{two-column matrix where each row represents and edge
  between two nodes. See function "MetaboSignal_matrix ( )".
  }
  \item{source_node}{character vector containing the KEGG ID of the node from which
  the shortest paths will be calculated.
  }
  \item{target_node}{character vector containing the KEGG ID of node to which the
  shortest path will be calculated.
  }
  \item{mode}{character constant indicating whether a directed or an undirected
  network will be considered. "all" indicates that all the edges of the network
  will be considered as undirected. "out" indicates that all the edges of the
  network will be considered as directed. "SP" indicates that all network will
  be considered as directed except the edges linked to target metabolite, which
  will be considered as undirected. The difference between the "out" and "SP"
  options, is that the latter aids reaching target metabolites that are substrate
  of irreversible reactions. By default, mode = "out".
  }
  \item{type}{indicates whether all shortest paths or a single shortest path will
  be considered when there are several shortest paths between the source_node and
  the target_node.  If type = "all", all shortest paths will be considered. If
  type = "first" a single path will be considered. If type = "bw" the path with
  the highest betweenness score will be considered. The betweenness score is
  calculated as the average betweenness of the gene nodes of the path.
  Using type = "bw" increases the time required to compute this function. By
  default, type = "first".
  }
  \item{network}{logical constant indicating the output of the function. If network
  = TRUE, the shortest path(s) will be represented as a two-column matrix where
  each row represents and egde between two nodes. If network = FALSE, the shortest
  path(s) will be represented as a vector (single path), or as a matrix (multiple paths),
  where each row represents a shortest path. By default, network = FALSE.
  }
  \item{interaction_type}{optional three-column matrix generated using the function
  "MS_interactionType( )". This argument is ignored when network = FALSE.
  }
}
\value{
By default, a vector or a matrix where each row contains a shortest path from the
source_node to the target_node. If network = TRUE, the function returns matrix where
each row represents an edge between two nodes. KEGG IDs can be transformed into
common names using the function "MS_ChangeNames( )".
}
\references{
G. Csardi and T. Nepusz (2015). igraph package, The Comprehensive R Archive
Network, v1.0.1.
}
\examples{
data(MetaboSignal_table)

# Shortest path from HK ("K00844") to a-D-Glucose ("cpd:C00267")

path1 <- MS_GetShortestpaths(MetaboSignal_table, "K00844", "cpd:C00267", mode = "SP")
path2 <- MS_GetShortestpaths(MetaboSignal_table, "K00844", "cpd:C00267", mode = "out")
\donttest{
# Shortest paths from G6PC ("K01084") to pyruvate ("cpd:C00022")

path3 <- MS_GetShortestpaths(MetaboSignal_table, "K01084", "cpd:C00022", type = "all")
path4 <- MS_GetShortestpaths(MetaboSignal_table, "K01084", "cpd:C00022", type = "bw")

# Shortest paths from G6PC ("K01084") to pyruvate ("cpd:C00022") as network

path3 <- MS_GetShortestpaths(MetaboSignal_table, "K01084", "cpd:C00022", type = "all",
                             network = TRUE)
}
}


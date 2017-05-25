\name{MS_keggNetwork}
\alias{MS_keggNetwork}
\title{Build MetaboSignal network-table}
\description{
This function generates a directed network-table (i.e. three-column matrix), where
each row represents an edge connecting two nodes (from source to target). The third
column indicates the type of interaction.Nodes represent four different molecular
entities: metabolic-genes (i.e. genes encoding enzymes that catalyze metabolic
reactions), signaling-genes (e.g. kinases), reactions and compounds (metabolites,
drugs or glycans).

The network-table generated with this function can be customized based on several
criteria. For instance, undesired nodes can be removed or replaced using the
functions "MS_removeNode( )" or "MS_replaceNode( )" respectively. Also, the
network can be filtered according to different topological parameters (e.g. node
betweenness) using the function "MS_filterNetwork( )".
}
\usage{
MS_keggNetwork(metabo_paths, signaling_paths, expand_genes = FALSE)
}
\arguments{
  \item{metabo_paths}{character vector containing the KEGG IDs of the metabolic
  pathways of interest (organism-specific). For example, the KEGG ID for the
  pathway "glycolysis/gluconeogenesis" in the rat is"rno00010". See the function
  "MS_keggFinder( )".
  }
  \item{signaling_paths}{character vector containing the KEGG IDs for the
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
\note{
Reaction directionality reported in KEGG has been cross-validated with published
literature (Duarte et al., 2007).
For examples of metabolic and signaling KEGG pathway check the dataset rno_pathways.
}
\references{
Carlson, M. org.Hs.eg.db: Genome wide annotation for Human.R package version 3.2.3.

Davidovic, L., Navratil, V., Bonaccorso, C.M., Catania, M.V., Bardoni, B. & Dumas,
M.E. (2011). A metabolomic and systems biology perspective on the brain of the
fragile X syndrome mouse model. Genome Research, 21, 2190-2202.

Duarte, N.C., Becker, S.A., Jamshidi, N., Thiele, I., Mo, M.L., Vo, T.D., Srivas,
R. & Palsson, O.B. (2007). Global reconstruction of the human metabolic network
based on genomic and bibliomic data. Proceedings of the National Academy of Sciences,
104, 1777-1782.

Mark, A., Thompson, R. and Wu, C. (2014) mygene: Access MyGene.Info_ services.R
package version 1.6.0.

Posma, J.M., Robinette, S.L., Holmes, E. & Nicholson, J.K. (2014). MetaboNetworks,
an interactive Matlab-based toolbox for creating, customizing and exploring
sub-networks from KEGG. Bioinformatics, 30, 893-895.

Zhang, J.D. & Wiemann, S. (2009). KEGGgraph: a graph approach to KEGG PATHWAY in
R and Bioconductor. Bioinformatics, 25, 1470-1471.

http://www.kegg.jp/kegg/docs/keggapi.html
}
\examples{\donttest{
# MetaboSignal network-table with organism-specific gene nodes

MS_netIsoforms <- MS_keggNetwork(metabo_paths = c("rno00010", "rno00562"),
                                   signaling_paths = c("rno04910", "rno04151"),
                                   expand_genes = TRUE)

# MetaboSignal network-table with orthology gene nodes

MS_net <- MS_keggNetwork(metabo_paths = c("rno00010", "rno00562"),
                         signaling_paths = c("rno04910", "rno04151"))

}
}




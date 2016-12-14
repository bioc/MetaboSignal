#################### get_interactiontype ####################
get_interactiontype = function(path) {
    file = paste("http://rest.kegg.jp/get/", path, "/kgml", sep = "")
    interactions = parseKGML2DataFrame(file, reactions = FALSE)

    if(nrow(interactions) == 0) {
        return(NULL)
    } else {
        return(interactions)
    }
}

#################### collapse_interactions ####################
collapse_interactions = function(interaction, all_interactions) {
    ind = which(all_interactions[, 1] == interaction[1] &
                    all_interactions[, 2] == interaction[2])
    subtypes = sort(unique(all_interactions[ind, 3]))
    subtype = paste(subtypes, collapse = ";")
    subtype = gsub(" ", "_", subtype)
    new_line = c(interaction[1], interaction[2], subtype)

}

#################### MS_interactionType ####################
MS_interactionType = function(paths, expand_genes = FALSE) {

    ## Check that all paths belong to the same organism
    input_paths = unique(paths)
    organism_code = substr(input_paths, 1, 3)
    organism_code = unique(organism_code)

    if (length(organism_code) > 1) {
        stop("All paths have to belong to the same organism: check path IDs")
    }

    ## Get all paths of the organism of interest
    lines = try (keggList("pathway", organism_code), silent = TRUE)

    if (grepl("Error", lines)[1]) { # example: metabo_paths = 'X'
        stop("Incorrect path IDs")
    }
    all_paths = substr(names(lines), 6, 13)

    ## Paths included
    paths_org = intersect(input_paths, all_paths)

    if (length(paths_org) == 0) {
        stop ("Path IDs are incorrect")
    }

    paths_orgDF = do.call(rbind, lapply(paths_org, get_interactiontype))

    if(is.null(paths_orgDF)) {
        stop ("Impossible to build an interaction type table with the current
              path IDs")
    }
    interactionM = as.matrix(paths_orgDF, ncol = 3)
    rownames(interactionM) = NULL

    if (expand_genes == FALSE) {
        to_print = ("Transforming gene IDs into orthology IDs")
        message(to_print)

        file_ko = paste("http://rest.kegg.jp/link/ko/", organism_code, sep = "")
        response_ko = getURL(file_ko)
        koTable = convertTable(response_ko)
        koTable[, 2] = substr(koTable[, 2], 4, 9)

        for (r in 1:nrow(interactionM)) {
            for (c in 1:ncol(interactionM[, 1:2])) {
                index = which(koTable[, 1] == interactionM[r, c])
                if (length(index) > 0) {
                    interactionM[r, c] = koTable[index, 2]
                }
            }
        }
    }

    interactionM = unique(interactionM)

    ## Collapse interactions ##
    IM_lines = split(interactionM, row(interactionM))
    interactionM_collapsed = unique(do.call(rbind, lapply(IM_lines, collapse_interactions,
                                                          interactionM)))
    interactionM_collapsed = as.matrix(interactionM_collapsed, ncol = 3)
    colnames(interactionM_collapsed) = c("node1", "node2", "interaction_type")
    rownames(interactionM_collapsed) = NULL

    return(interactionM_collapsed)
}

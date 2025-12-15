
runGSEA <-
  function(geneList,OrgDb,keyType,

           ontologies   = c("BP","MF","CC"),
           minGSSize    = 10,
           maxGSSize    = 500,
           pvalueCutoff = 0.05,
           seed         = TRUE,
           by           = "fgsea",
           nPermSimple  = 10000,
           eps          = 1e-50,
           verbose      = FALSE,
           ...
  ) {


    stopifnot(all(ontologies %in% c("BP","MF","CC")))

    # return a list of results, one per ontology
    sapply(ontologies,
           function(ont) {
             clusterProfiler::gseGO(
               geneList = geneList,
               OrgDb = OrgDb,
               keyType = keyType,

               ont=ont,
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               seed         = TRUE,
               by           = "fgsea",
               nPermSimple  = 10000,
               eps          = 1e-50,
               verbose      = FALSE)
           },simplify=FALSE)

  }

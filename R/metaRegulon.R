#' Inferring the gene regulatory network and estimating the transcription factor activity using metaRegulon
#'
#' Inferring the gene regulatory network and estimating the transcription factors activity using metaRegulon.
#'
#' @param obj An \code{SingleCellExperiment} object.
#' @param use_species The species used for data analysis. Default: human.
#' @param project_name The name of present project. Default: Regulon_analysis.
#' @param prior_net The prior-network used for regulon analysis. To use a custom
#' network, please use standard three-column data.frame with regulators, targets,
#' and weights. Default: NULL, pan-tissue network.
#' @param tf_gene The TF gene symbols used for regulon analysis. To use a custom
#' TF list, please ensure to match the symbol in expression matrix and prior-network.
#' Default: NULL, human TF symbol.
#' @param grn_net_type The method used for extracting TF-target network from gene
#' regulatory networks. The "auto" representing "aracne", while the "top5" and
#' "top10" representing the top 5% and 10% targets of each regulator are retained
#' for regulon network integration. Default: "auto".
#' @param force_repeat Whether to overwrite the files for the intermediate results.
#' If you want to continue with the results calculated from last time, please set
#' this option to FALSE. Default: FALSE.
#' @param out_path The output files for the project.
#' @param ncores number of cores used for present project.
#' @param ... Other parameters passed to scATFR.

#'
#' @return An \code{SingleCellExperiment} object.
#'
#' @import scATFR
#'
#' @export
#'
#' @examples
#' data(test_obj)
#' test_obj <- metaRegulon(obj = test_obj, use_species = "mouse", ncores = 4, force_repeat = TRUE)
#'
#'
metaRegulon <- function(obj,
                        use_species = c("human","mouse","other"),
                        project_name = "Regulon_analysis",
                        prior_net = NULL,
                        tf_gene = NULL,
                        grn_net_type = c("auto","top5","top10"),
                        force_repeat = FALSE,
                        out_path = "./",
                        ncores=1, ...){
  use_species <- match.arg(use_species)
  grn_net_type <- match.arg(grn_net_type)
  out_path <- file.path(out_path,project_name)
  if(!is(obj,"SingleCellExperiment")) stop("The 'obj' must be a SingleCellExperiment object!")
  #--- 1. load data
  message("Step 1. Loading data...")

  if(use_species == "human"){
    tf_gene <- TF_Symbols$human
    if(is.null(prior_net)){
      data("human_pantissue_net")
      prior_net <- human_pantissue_net
    }else{
      prior_net <- as.data.frame(prior_net)
      colnames(prior_net)[1:3] <- c("tf","target","weight")
    }
    prior_net <- prior_net[as.character(prior_net$tf) %in% tf_gene,]
  }else if(use_species == "mouse"){
    tf_gene <- TF_Symbols$mouse
    if(is.null(prior_net)){
      data("mouse_pantissue_net")
      prior_net <- mouse_pantissue_net
    }else{
      prior_net <- as.data.frame(prior_net)
      colnames(prior_net)[1:3] <- c("tf","target","weight")
    }
    prior_net <- prior_net[as.character(prior_net$tf) %in% tf_gene,]
  }else if(use_species == "other"){
    if(is.null(prior_net) || is.null(tf_gene)){
      stop("The 'prior_net' and 'tf_gene' should not be NULL when 'use_species' is set to other")
    }
    prior_net <- prior_net[as.character(prior_net$tf) %in% tf_gene,]
    if(nrow(prior_net)==0) stop("No gene symbol detected in 'prior_net', please check your 'tf_gene'!")
  }

  if(length(intersect(tf_gene,row.names(obj)))==0){
    stop("No TF gene symbol detected in your data, please check the row names of your expression matrix!")
  }
  dir.create(path = out_path, recursive = TRUE,showWarnings = FALSE)
  saveRDS(obj,file.path(out_path,"01_object.rds"))

  #--- inferring GRNs
  message("Step 2. inferring GRNs ...")
  if(("02_inferred_GRNs.rds" %in% list.files(out_path)) && !force_repeat){
    message("File '02_inferred_GRNs.rds' exist in your output dir, jumping this step...")
    inferred_grns <- readRDS(file.path(out_path,"02_inferred_GRNs.rds"))
  }else{
    obj <- scATFR::inferGRNs(x = obj, method="puic", ncores=ncores, ...)
    inferred_grns <- as.matrix(scATFR::GRN(obj))
    saveRDS(inferred_grns,file.path(out_path,"02_inferred_GRNs.rds"))
    saveRDS(obj,file.path(out_path,"01_object.rds"))
  }

  #--- integrating GRNs
  message("Step 3. integrating regulons...")
  if(("03_integrated_GRNs.rds" %in% list.files(out_path)) && !force_repeat){
    message("File '03_integrated_GRNs.rds' exist in your output dir, jumping this step...")
    integrated_grns <- readRDS(file.path(out_path,"03_integrated_GRNs.rds"))
  }else{
    if(grn_net_type == "auto"){
      grn_net <- PUIC::matToNet(weightMat = inferred_grns, methods = "aracne",cutoff = 0)
    }else if(grn_net_type == "top5"){
      grn_net <- GENIE3::getLinkList(weightMatrix = inferred_grns)
      top_num <- round(0.05*ncol(grn_net))
      grn_net <- grn_net %>% dplyr::group_by(regulatoryGene) %>% dplyr::top_n(top_num,weight)
    }else if(grn_net_type == "top10"){
      grn_net <- GENIE3::getLinkList(weightMatrix = inferred_grns)
      top_num <- round(0.1*ncol(grn_net))
      grn_net <- grn_net %>% dplyr::group_by(regulatoryGene) %>% dplyr::top_n(top_num,weight)
    }
    colnames(grn_net)[1:3] <- c("tf","target","weight")
    grn_net <- grn_net[grn_net$tf %in% tf_gene,]
    integrated_grns <- scATFR::integraNets(grnNets = list(prior_net, grn_net), join_type = "full",
                                           min_target_num = 10L,norm_weight = TRUE, return_sig=TRUE,
                                           maxInter = 2)
    saveRDS(integrated_grns,file.path(out_path,"03_integrated_GRNs.rds"))
  }

  #--- filtering regulons
  message("Step 4. filtering positive regulons...")
  if(("04_filtered_regulons.rds" %in% list.files(out_path)) && !force_repeat){
    message("File '04_filtered_regulons.rds' exist in your output dir, jumping this step...")
    filtered_regulons <- readRDS(file.path(out_path,"04_filtered_regulons.rds"))
  }else{
    integrated_regulons <- scATFR::dfToList(df = integrated_grns)
    obj <- scATFR::filterRegulons(x = obj, gene_list = integrated_regulons, maxSize=1000L, ncores=ncores)
    filtered_regulons <- as.list(scATFR::Regulon(x = obj))
    saveRDS(filtered_regulons,file.path(out_path,"04_filtered_regulons.rds"))
    saveRDS(obj,file.path(out_path,"01_object.rds"))
  }

  #--- activity evaluation
  message("Step 5. estimating TF activity...")
  obj <- scATFR::regulonActivity(x = obj, gene_list = filtered_regulons)
  act_mat <- assay(altExp(obj))
  saveRDS(act_mat, file.path(out_path,"05_TF_activity_matrix.rds"))
  saveRDS(obj, file.path(out_path,"01_object.rds"))
  message("Done!")
  obj
}



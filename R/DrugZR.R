#' DrugZR
#' Identify drug-gene interactions in paired sample genomic perturbation screens
#'
#' Version 1.0.0
#'
#' Last Modified Oct 14, 2024
#'
#'
#' inputs: SummarizeExperiment (se) / DrugZ countTable
#'
#' input parameters:
#'
#' pseduocount:Constant value added to all reads (default 5) - prevents log(0) problems
#'
#' min_observations: An integer value to act as a threshold for the minimum number observations to be included in the analysis (default=1)
#'
#' seTrue: whether the input is SummarizedExperiment (default TRUE)
#'
#'
#'
#' @export
#' @import dplyr tidyr


DrugZR<-function(se,control,treatment,min_observations=1,pseudocount=5,seTrue=TRUE){


  if(seTrue == TRUE){
    print("converting se to DrugZR input Tab")
    reads<-.DrugZCT(se)
    head(reads)
  } else {
    reads<-se
  }

  #set rownames
  rownames(reads)<-reads$sgRNA
  reads<-reads[order(reads$sgRNA),]

  #Make sure the length of control & treatment is the same
  if(length(control) != length(treatment)){
    stop("Must have the same number of control and drug samples to run the paired approach")
  } else {
    replicate<-1:length(control)
  }

  print(paste0("Controls : ",paste(control,collapse = ", ")))
  print(paste0("Treatments : ",paste(treatment,collapse = ", ")))

  #Normalize Samples Analyzed to read/10 Million
  normalized_counts<-.normalize_readcounts(reads,control,treatment)


  #List of foldchange per replicate
  print("Caculating guide-level Zscores for replicates")
  foldc<-list()
  foldchangez<-data.frame(sgRNA=reads$sgRNA,GENE=reads$Gene)

  for(i in seq_along(replicate)){

    fc<-.calculate_fold_change_eb(reads, normalized_counts, control, treatment, pseudocount, replicate=i,half_window_size=500)
    #fc<-fc[order(fc$GENE),]

      foldc[[paste0("rep",i)]]<-fc

      foldchangez[paste0("zscore_fc_",i)]<-fc[,grep("zscore_fc",colnames(fc))]
  }

  bar_lfc<-data.frame(sgRNA=reads$sgRNA,Gene=reads$Gene)
   for(l in seq_along(foldc)){
     col<-grep("eb_std|zscore_fc|fc",colnames(foldc[[l]]))
     bar_lfc<-cbind(bar_lfc,foldc[[l]][,col])
   }

  print("Calculating gene-level Zscores")
  columns<-colnames(foldchangez)[grep("zscore_fc_",colnames(foldchangez))]
  GeneLevel<-.calculate_drugz_score(fold_change=foldchangez, min_observations, columns=columns)

  return(list(bar_lfc,GeneLevel))

  }



#Normalize ReadCounts to reads/10^7

.DrugZCT<-function(se){
  countmat<-as.matrix(assays(se)$counts)
  se$TREATMENT_NAME <- gsub(" ", "", se$TREATMENT_NAME)
  colnames(countmat)<-se$TREATMENT_NAME
  se
  sgRNA<-rownames(countmat)
  Gene<- str_remove(sgRNA,'_[:alnum:]+$')

  dzct<-cbind(sgRNA,Gene,countmat)

  return(dzct)

}

.normalize_readcounts <- function(reads, treatment, control) {
  # Assuming norm_value is a global variable or predefined
  norm_value <- 1e7

  # Select the columns for treatment and control groups
  reads_to_normalize <- reads[, c(control, treatment)]
  rownames(reads_to_normalize)<-reads$sgRNA

  # Normalize raw read counts using norm_value
  # Sum across columns, replicate the sum to match dimensions, and perform element-wise division
  column_sums <- colSums(reads_to_normalize)
  normalized_counts <- (norm_value * reads_to_normalize) / matrix(rep(column_sums, each = nrow(reads_to_normalize)), nrow = nrow(reads_to_normalize))

  return(normalized_counts)
}

#calculate fold change per replicate
#Fold Change & empirical bayes
.calculate_fold_change_eb <- function(reads, normalized_counts, control, treatment, pseudocount,replicate,half_window_size=500) {

  # Create a dataframe with guide ids as row names
  fold_change <- data.frame(GENE = reads[["Gene"]], row.names = rownames(reads))

  # Generate foldchange, estimated_variance, and foldchange zscore column ids for each replicate
  fc_replicate_id <- sprintf("fc_%s", replicate)
  fc_zscore_id <- sprintf("zscore_fc_%s", replicate)
  empirical_bayes_id <- sprintf("eb_std_%s", replicate)

  # Get the control and treatment sample ids for each replicate
  control_sample <- control[replicate]
  treatment_sample <- treatment[replicate]

  # Add the control sample column to the fold change dataframe and sort by this column
  fold_change[control_sample] <- reads[,control_sample]
  fold_change[treatment_sample] <- reads[,treatment_sample]
  fold_change <- fold_change[order(fold_change[[control_sample]], decreasing = TRUE),]

  # Extract the number of rows (number of guides) of the reads dataframe
  no_of_guides <- nrow(reads)

  # Fill in the estimated_variance and foldchange_zscore columns with 0s for now
  fold_change[empirical_bayes_id] <- rep(0, no_of_guides)
  fold_change[fc_zscore_id] <- rep(0, no_of_guides)

  # Calculate the log2 ratio of treatment normalized read counts to control - foldchange
  normalized_counts<-normalized_counts[rownames(fold_change),]
  fold_change[fc_replicate_id] <- log2((normalized_counts[treatment_sample] + pseudocount) / (normalized_counts[control_sample] + pseudocount))


  # Calculate the standard deviation of foldchange based on a 2 * half_window_size range
  std_dev <- sd(na.omit(as.matrix(fold_change[1:(half_window_size * 2), fc_replicate_id])))
  fold_change[1:half_window_size, empirical_bayes_id] <- std_dev


  # Iterate in a range(half_window_size, no_of_guides - half_window_size, 25)
  for (i in seq(half_window_size, no_of_guides - half_window_size + 25, by = 25)) {

    # In every bin calculate stdev
    std_dev <- sd(na.omit(as.matrix(fold_change[(i - half_window_size + 1):(i + half_window_size), fc_replicate_id])))


    if(std_dev >= fold_change[i - 1, empirical_bayes_id]){
      fold_change[i:min(i + 25, no_of_guides), empirical_bayes_id]<-std_dev
    } else {
      fold_change[i:min(i + 25, no_of_guides), empirical_bayes_id]<-fold_change[i - 1, empirical_bayes_id]
    }
  }



  # Get the variation estimate for the final bin and set the remaining values in the empirical bayes column
  final_std_dev <- fold_change[no_of_guides - half_window_size, empirical_bayes_id]
  fold_change[(no_of_guides - half_window_size + 1):no_of_guides, empirical_bayes_id] <- final_std_dev

  # Calculate the z_score for each guide (fc/eb_std)
  fold_change[, fc_zscore_id] <- fold_change[, fc_replicate_id] / fold_change[, empirical_bayes_id]
  fold_change<-fold_change[order(rownames(fold_change)),]

  return(fold_change)
}



.calculate_drugz_score <- function(fold_change, min_observations, columns) {
  # Aggregate z-scores per gene, replacing NAs with zeros

  sumZs<-c()
  numObvs<-c()
  for(i in seq_along(columns)){
    tmp<-fold_change %>%
      group_by(GENE) %>%
      summarize(sumZ = sum(.data[[columns[i]]]),
                numObs=sum(!is.na(.data[[columns[i]]]) & .data[[columns[i]]] != 0))

    sumZs<-cbind(sumZs,tmp$sumZ)
    numObvs<-cbind(numObvs,tmp$numObs)

    rownames(sumZs)<-tmp$GENE
    rownames(numObvs)<-tmp$GENE
  }

  per_gene_scores<-data.frame(sumZ=rowSums(sumZs),numObs=rowSums(numObvs))

  # Filter genes based on the minimum number of observations
  per_gene_results <- filter(per_gene_scores, numObs >= min_observations)

  # Calculate normalized gene z-score
  per_gene_results <- per_gene_results %>%
    mutate(normZ = scale(sumZ / sqrt(numObs), center = TRUE, scale = TRUE))

  colnames(per_gene_results)[grep("normZ",colnames(per_gene_results))]<-"normZ"

  # Calculate p-values and FDR for synthetic interactions
  per_gene_results <- per_gene_results %>%
    arrange(normZ) %>%
    mutate(pval_synth = pnorm(-normZ, lower.tail = FALSE),
           rank_synth = row_number(),
           fdr_synth = p.adjust(pval_synth, method = "BH"))

  per_gene_results<-per_gene_results[order(per_gene_results$normZ),]
  pnorm((per_gene_results$normZ * -1), lower.tail = FALSE)


  # Calculate p-values and FDR for suppressor interactions
  per_gene_results <- per_gene_results %>%
    arrange(desc(normZ)) %>%
    mutate(pval_supp = pnorm(normZ, lower.tail = FALSE),
           rank_supp = row_number(),
           fdr_supp = p.adjust(pval_supp, method = "BH"))

  # Return the results sorted by normalized z-scores (ascending)
  arrange(per_gene_results, normZ)
}



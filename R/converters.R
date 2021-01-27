#' Import Skyline files
#'
#' @param input name of MSstats input report from Skyline, which includes feature-level data.
#' @param annotation name of 'annotation.txt' data which includes Condition, BioReplicate, Run. If annotation is already complete in Skyline, use annotation=NULL (default). It will use the annotation information from input.
#' @param removeiRT TRUE (default) will remove the proteins or peptides which are labeld 'iRT' in 'StandardType' column. FALSE will keep them.
#' @param filter_with_Qvalue TRUE(default) will filter out the intensities that have greater than qvalue_cutoff in DetectionQValue column. Those intensities will be replaced with zero and will be considered as censored missing values for imputation purpose.
#' @param qvalue_cutoff Cutoff for DetectionQValue. default is 0.01.
#' @param ... additional parameters to `data.table::fread`.
#'
#' @return data.frame with the required format of MSstats.
#'
#' @author Meena Choi, Olga Vitek
#'
#' @export
#'
SkylinetoMSstatsGlycansFormat = function(
    input, removeiRT = TRUE, filter_with_Qvalue = TRUE,
    qvalue_cutoff = 0.01, removeFewMeasurements = TRUE,
    removeProtein_with1Feature = FALSE, use_log_file = TRUE, append = FALSE,
    verbose = TRUE, log_file_path = NULL, session_info_path = NULL, ...
) {
    MSstatsConvert::MSstatsLogsSettings(use_log_file, append, verbose,
                                        log_file_path)
    MSstatsConvert::MSstatsSaveSessionInfo(session_info_path, append = TRUE)

    input = MSstatsConvert::MSstatsImport(list(input = input),
                                          "MSstatsGlycans", "Skyline", ...)
    input = MSstatsConvert::MSstatsClean(input)
    input$IsotopeLabelType = "L"
    annotation = MSstatsConvert::MSstatsMakeAnnotation(input, NULL,
                                                       "Run" = "FileName")

    irt_filter = list(col_name = "StandardType",
                      filter_symbols = "iRT",
                      filter = removeiRT,
                      behavior = "remove",
                      fill_value = NULL,
                      drop_column = FALSE)

    truncated_filter = list(col_name = "Truncated",
                            filter_symbols = "TRUE",
                            behavior = "fill",
                            fill_value = NA_real_,
                            filter = TRUE,
                            drop_column = TRUE)

    qval_filter = list(score_column = "DetectionQValue",
                       score_threshold = qvalue_cutoff,
                       direction = "smaller",
                       behavior = "fill",
                       fill_value = 0,
                       handle_na = "keep",
                       filter = filter_with_Qvalue,
                       drop_column = TRUE)

    feature_columns = c("PeptideSequence", "PrecursorCharge",
                        "FragmentIon", "ProductCharge")
    input = MSstatsConvert::MSstatsPreprocess(
        input,
        annotation,
        feature_columns,
        remove_shared_peptides = TRUE,
        remove_single_feature_proteins = removeProtein_with1Feature,
        score_filtering = list(qval = qval_filter),
        exact_filtering = list(irt = irt_filter,
                               truncated = truncated_filter),
        aggregate_isotopic = FALSE,
        feature_cleaning = list(
            remove_features_with_few_measurements = removeFewMeasurements,
            summarize_multiple_psms = sum))
    input = MSstatsConvert::MSstatsBalancedDesign(input, feature_columns)

    msg_final = paste("** Finished preprocessing. The dataset is ready",
                      "to be processed by the dataProcess function.")
    getOption("MSstatsLog")("INFO", msg_final)
    getOption("MSstatsMsg")("INFO", msg_final)
    getOption("MSstatsLog")("INFO", "\n")
    input
}

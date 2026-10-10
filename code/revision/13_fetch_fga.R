# Fraction genome altered (FGA) for the TCGA PanCanAtlas 2018 studies, used
# as the scalar copy number baseline requested by Reviewer 2.
#
# Source: cBioPortal public API, sample attribute FRACTION_GENOME_ALTERED of
# the 32 *_tcga_pan_can_atlas_2018 studies. cBioPortal computes it from each
# study's segment file as the total length of segments with
# |segment mean| >= 0.2 divided by the total length of the measured segments
# (cbioportal-core, FractionGenomeAlteredCalculator), rounded to four decimals.
# One sample per patient: the primary tumour (sample type 01, or 03 for LAML)
# or, when there is none, the metastasis (06). The second case concerns almost
# only SKCM, whose TCGA samples are mostly metastatic: without it FGA would be
# missing for most melanomas, i.e. not at random.
#
# Run from the repository root (needs internet access; the API is occasionally
# unavailable, so each request is retried):
#   Rscript code/revision/13_fetch_fga.R
# Output: data/processed/fga_pancanatlas.csv and fga_pancanatlas_provenance.txt

suppressMessages({library(httr); library(jsonlite)})
api <- "https://www.cbioportal.org/api"
out <- "data/processed"
studies <- paste0(c("acc", "blca", "brca", "cesc", "chol", "coadread", "dlbc", "esca", "gbm",
                    "hnsc", "kich", "kirc", "kirp", "laml", "lgg", "lihc", "luad", "lusc",
                    "meso", "ov", "paad", "pcpg", "prad", "sarc", "skcm", "stad", "tgct",
                    "thca", "thym", "ucec", "ucs", "uvm"), "_tcga_pan_can_atlas_2018")

get_json <- function(url, tries = 8) {
  for (i in seq_len(tries)) {
    res <- tryCatch(GET(url, accept_json(), timeout(90)), error = function(e) NULL)
    if (!is.null(res) && status_code(res) == 200) {
      return(fromJSON(content(res, "text", encoding = "UTF-8")))
    }
    Sys.sleep(5 * i)
  }
  stop("Request failed after ", tries, " attempts: ", url)
}

fga <- do.call(rbind, lapply(studies, function(study) {
  x <- get_json(paste0(api, "/studies/", study,
                       "/clinical-data?clinicalDataType=SAMPLE&attributeId=FRACTION_GENOME_ALTERED"))
  data.frame(patient = x$patientId, sample = x$sampleId, study = study,
             fga = as.numeric(x$value), stringsAsFactors = FALSE)
}))
stopifnot(!anyNA(fga$fga), all(fga$fga >= 0 & fga$fga <= 1))
fga$sample_type <- substr(fga$sample, 14, 15)
primary <- fga[fga$sample_type %in% c("01", "03"), ]
if (anyDuplicated(primary$patient)) stop("More than one primary sample per patient.")
metastasis <- fga[fga$sample_type == "06" & !fga$patient %in% primary$patient, ]
if (anyDuplicated(metastasis$patient)) stop("More than one metastatic sample per patient.")
chosen <- rbind(primary, metastasis)
chosen <- chosen[order(chosen$patient), c("patient", "sample", "sample_type", "study", "fga")]
utils::write.csv(chosen, file.path(out, "fga_pancanatlas.csv"), row.names = FALSE)

meta <- get_json(paste0(api, "/studies/", studies[1], "/clinical-attributes/FRACTION_GENOME_ALTERED"))
writeLines(c(
  paste0("Retrieved: ", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  paste0("API: ", api, " (sample attribute ", meta$clinicalAttributeId, ", '", meta$displayName, "')"),
  paste0("Studies: ", length(studies), " TCGA PanCanAtlas 2018 studies"),
  paste0("Samples returned: ", nrow(fga), " (sample types: ",
         paste(names(table(fga$sample_type)), table(fga$sample_type), sep = " = ", collapse = ", "), ")"),
  paste0("One sample per patient: primary tumour (01; 03 for LAML) for ", nrow(primary),
         " patients; metastasis (06) when no primary was available for ", nrow(metastasis),
         " patients (", paste(names(table(metastasis$study)), table(metastasis$study), sep = " = ",
                              collapse = ", "), ")"),
  "Definition: total length of segments with |segment mean| >= 0.2 divided by the total length of measured segments, from each study's segment file (cbioportal-core FractionGenomeAlteredCalculator; cbioportal v5.4.0 DaoCopyNumberSegment)."),
  file.path(out, "fga_pancanatlas_provenance.txt"))
message("FGA for ", nrow(chosen), " patients written to ", file.path(out, "fga_pancanatlas.csv"))

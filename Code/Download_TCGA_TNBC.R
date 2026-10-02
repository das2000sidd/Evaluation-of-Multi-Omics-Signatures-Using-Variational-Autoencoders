setwd("~/Desktop/Related_to_Haider_Lab/TCGA_TNBC")

library(TCGAbiolinks)
library(dplyr)

# 1) download and prepare expression (example: HTSeq counts)
query.rna <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = "RNA-Seq",
  sample.type = "Primary Tumor"
)



rna.results <- getResults(query.rna)

dim(rna.results)
table(rna.results$sample_type)
table(rna.results$experimental_strategy)


GDCdownload(
  query = query.rna,
  directory = "TCGA_BRCA"
)

rna <- GDCprepare(
  query = query.rna,
  directory = "TCGA_BRCA"
)

list.files(
  "TCGA_BRCA",
  recursive = TRUE,
  full.names = TRUE
)[1:20]


length(
  list.files(
    "TCGA_BRCA",
    recursive = TRUE,
    full.names = TRUE
  )
)

object.size(query.rna)


rna.files <- list.files(
  "TCGA_BRCA/TCGA-BRCA",
  recursive = TRUE,
  full.names = TRUE
)

length(rna.files)

rna.files[grep("RNA|Gene_Expression|STAR", rna.files)][1:20]


rna.files <- rna.files[
  grepl(
    "Transcriptome_Profiling|Gene_Expression_Quantification",
    rna.files
  )
]

length(rna.files)
head(rna.files)

list.dirs(
  "TCGA_BRCA/TCGA-BRCA",
  recursive = TRUE
)[1:10]

length(list.dirs(
  "TCGA_BRCA/TCGA-BRCA",
  recursive = TRUE
))

rna.test <- read.delim(
  rna.files[1],
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  check.names = FALSE
)

dim(rna.test)
colnames(rna.test)
head(rna.test)


rna.results <- getResults(query.rna)

dim(rna.results)
colnames(rna.results)

grep(
  "file|case|sample|submitter|barcode",
  colnames(rna.results),
  ignore.case = TRUE,
  value = TRUE
)

rna.map <- rna.results[, c(
  "file_id",
  "file_name",
  "cases.submitter_id",
  "sample.submitter_id",
  "sample_type"
)]

head(rna.map)


length(unique(rna.map$cases.submitter_id))
length(unique(rna.map$sample.submitter_id))

table(duplicated(rna.map$cases.submitter_id))

table(table(rna.map$cases.submitter_id))

head(basename(rna.files))
head(rna.map$file_name)


sum(basename(rna.files) %in% rna.map$file_name)


wes.patients <- unique(maf$patient_id)
length(wes.patients)



rna.patients <- unique(rna.map$cases.submitter_id)

length(rna.patients)
length(intersect(wes.patients, rna.patients))

dup.patients <- names(
  table(rna.map$cases.submitter_id)[
    table(rna.map$cases.submitter_id) > 1
  ]
)

rna.map[
  rna.map$cases.submitter_id %in% dup.patients,
  c(
    "file_id",
    "file_name",
    "cases.submitter_id",
    "sample.submitter_id",
    "sample_type"
  )
]

table(rna.map$sample_type)


id="53184"
sample.counts <- table(rna.map$sample.submitter_id)

sample.counts[sample.counts > 1]


dup.samples <- names(sample.counts[sample.counts > 1])

rna.map[
  rna.map$sample.submitter_id %in% dup.samples,
  c(
    "file_id",
    "file_name",
    "cases.submitter_id",
    "sample.submitter_id"
  )
]

rna.map[
  rna.map$sample.submitter_id %in% dup.samples,
  c(
    "file_id",
    "file_name",
    "cases.submitter_id",
    "sample.submitter_id"
  )
]

## determine whether the duplicate files are identical
for (s in dup.samples) {
  
  rows <- rna.map$sample.submitter_id == s
  
  cat("\n", s, "\n")
  print(
    rna.map[
      rows,
      c(
        "file_id",
        "file_name",
        "cases.submitter_id",
        "sample.submitter_id"
      )
    ]
  )
}

for (s in dup.samples) {
  
  rows <- which(rna.map$sample.submitter_id == s)
  
  f1 <- rna.files[
    basename(rna.files) == rna.map$file_name[rows[1]]
  ]
  
  f2 <- rna.files[
    basename(rna.files) == rna.map$file_name[rows[2]]
  ]
  
  x1 <- read.delim(
    f1,
    header = TRUE,
    sep = "\t",
    comment.char = "#",
    check.names = FALSE
  )
  
  x2 <- read.delim(
    f2,
    header = TRUE,
    sep = "\t",
    comment.char = "#",
    check.names = FALSE
  )
  
  # Remove STAR summary rows
  x1 <- x1[!grepl("^N_", x1$gene_id), ]
  x2 <- x2[!grepl("^N_", x2$gene_id), ]
  
  same <- identical(
    x1$unstranded,
    x2$unstranded
  )
  
  cat(s, ":", same, "\n")
}

## compare their metadata and sequencing depth for duplicate files
for (s in dup.samples) {
  
  rows <- which(rna.map$sample.submitter_id == s)
  
  cat("\n====================\n")
  cat(s, "\n")
  
  for (i in rows) {
    
    f <- rna.files[
      basename(rna.files) == rna.map$file_name[i]
    ]
    
    x <- read.delim(
      f,
      header = TRUE,
      sep = "\t",
      comment.char = "#",
      check.names = FALSE
    )
    
    # STAR summary rows
    star <- x[grepl("^N_", x$gene_id), ]
    
    cat(
      "\nFile:", rna.map$file_name[i],
      "\nTotal unstranded counts:",
      sum(x$unstranded[!grepl("^N_", x$gene_id)]),
      "\n"
    )
    
    print(
      star[, c("gene_id", "unstranded")]
    )
  }
}

## For duplicate sample.submitter_id
qc <- data.frame()

for (s in dup.samples) {
  
  rows <- which(rna.map$sample.submitter_id == s)
  
  for (i in rows) {
    
    f <- rna.files[
      basename(rna.files) == rna.map$file_name[i]
    ]
    
    x <- read.delim(
      f,
      header = TRUE,
      sep = "\t",
      comment.char = "#",
      check.names = FALSE
    )
    
    star <- x[grepl("^N_", x$gene_id), ]
    gene <- x[!grepl("^N_", x$gene_id), ]
    
    qc <- rbind(
      qc,
      data.frame(
        sample = s,
        file_name = basename(f),
        assigned = sum(gene$unstranded),
        multimapping = star$unstranded[
          star$gene_id == "N_multimapping"
        ],
        noFeature = star$unstranded[
          star$gene_id == "N_noFeature"
        ],
        unmapped = star$unstranded[
          star$gene_id == "N_unmapped"
        ],
        ambiguous = star$unstranded[
          star$gene_id == "N_ambiguous"
        ]
      )
    )
  }
}

qc

qc$assigned_fraction <- qc$assigned / (
  qc$assigned +
    qc$multimapping +
    qc$noFeature +
    qc$unmapped +
    qc$ambiguous
)

best_dup <- qc[
  ave(
    qc$assigned_fraction,
    qc$sample,
    FUN = function(x) x == max(x)
  ),
]

bad_files <- qc$file_name[
  qc$assigned_fraction <
    ave(
      qc$assigned_fraction,
      qc$sample,
      FUN = max
    )
]

rna.map.clean <- rna.map[
  !rna.map$file_name %in% bad_files,
]

nrow(rna.map.clean)
length(unique(rna.map.clean$sample.submitter_id))

patient.counts <- table(rna.map.clean$cases.submitter_id)

multi.patients <- names(
  patient.counts[patient.counts > 1]
)

length(multi.patients)

rna.multi <- rna.map.clean[
  rna.map.clean$cases.submitter_id %in% multi.patients,
  c(
    "cases.submitter_id",
    "sample.submitter_id",
    "file_name"
  )
]

rna.multi[
  order(
    rna.multi$cases.submitter_id,
    rna.multi$sample.submitter_id
  ),
]

# 2) get clinical table (GDC clinical or legacy file)
clin <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")  # inspect columns

grep(
  "days|death|vital|survival",
  colnames(clin),
  ignore.case = TRUE,
  value = TRUE
)

clin[, grep(
  "days|death|vital|survival",
  colnames(clin),
  ignore.case = TRUE
)]

saveRDS(
  clin,
  file = "TCGA_BRCA_clinical.rds",
  compress = FALSE
)


query.wes <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

getResults(query.wes)

getResults(query.wes)[, c(
  "file_name",
  "data_type",
  "experimental_strategy",
  "workflow_type"
)]

GDCdownload(
  query.wes,
  method = "api",
  files.per.chunk = 5,
  directory = "TCGA_BRCA"
)

maf <- GDCprepare(
  query = query.wes,
  directory = "TCGA_BRCA"
)


query.wes <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open",
  experimental.strategy = "WXS"
)

nrow(getResults(query.wes))


wes.results <- getResults(query.wes)

dim(wes.results)

table(wes.results$data_type)
table(wes.results$experimental_strategy)
table(wes.results$data_format)

length(unique(wes.results$project))


maf$patient_id <- substr(
  maf$Tumor_Sample_Barcode,
  1,
  12
)

dim(maf)

length(unique(maf$patient_id))

table(maf$Variant_Classification)

## keep protein coding changes
keep <- maf$Variant_Classification %in% c(
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Splice_Site",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Translation_Start_Site",
  "Nonstop_Mutation"
)

maf2 <- maf[keep, ]

dim(maf2)

length(unique(maf2$patient_id))

table(maf2$Variant_Classification)


gene_counts <- sort(
  table(maf2$Hugo_Symbol),
  decreasing = TRUE
)

length(gene_counts)

head(gene_counts, 30)



## keep 1095 patients with uniqu WES and RNA
# Keep 01A when a patient has multiple primary-tumor samples
rna.final <- do.call(
  rbind,
  lapply(
    split(rna.map.clean, rna.map.clean$cases.submitter_id),
    function(x) {
      if (nrow(x) == 1) {
        return(x)
      }
      
      x[x$sample.submitter_id == 
          paste0(unique(x$cases.submitter_id), "-01A"), ]
    }
  )
)

rownames(rna.final) <- NULL

dim(rna.final)
length(unique(rna.final$cases.submitter_id))

table(table(rna.final$cases.submitter_id))

## overlap RNA and WES samples
rna.patients <- unique(rna.final$cases.submitter_id)
wes.patients <- unique(maf$patient_id)

matched.patients <- intersect(
  rna.patients,
  wes.patients
)

length(rna.patients)
length(wes.patients)
length(matched.patients)

##create the final 1,095-patient mapping
# Keep one sample per patient.
# For patients with multiple samples, retain 01A.

rna.final <- do.call(
  rbind,
  lapply(
    split(rna.map.clean, rna.map.clean$cases.submitter_id),
    function(x) {
      
      if (nrow(x) == 1) {
        return(x)
      }
      
      x[x$sample.submitter_id == 
          paste0(unique(x$cases.submitter_id), "-01A"), ]
    }
  )
)

rownames(rna.final) <- NULL

dim(rna.final)
length(unique(rna.final$cases.submitter_id))


## Build the expression matrix
# Get local file paths corresponding to the final 1095 samples
rna.final$file_path <- rna.files[
  match(
    rna.final$file_name,
    basename(rna.files)
  )
]

# Read first file to establish gene information
x <- read.delim(
  rna.final$file_path[1],
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  check.names = FALSE
)

# Remove STAR summary rows
x <- x[!grepl("^N_", x$gene_id), ]

# Keep gene ID, gene name and raw unstranded counts
gene_id <- x$gene_id
gene_name <- x$gene_name

# Initialize matrix
rna.counts <- matrix(
  NA_integer_,
  nrow = length(gene_id),
  ncol = nrow(rna.final)
)

rownames(rna.counts) <- gene_id
colnames(rna.counts) <- rna.final$cases.submitter_id

# First sample
rna.counts[, 1] <- x$unstranded

rm(x)


for (i in 2:nrow(rna.final)) {
  
  if (i %% 100 == 0) {
    cat("Processing", i, "of", nrow(rna.final), "\n")
  }
  
  x <- read.delim(
    rna.final$file_path[i],
    header = TRUE,
    sep = "\t",
    comment.char = "#",
    check.names = FALSE
  )
  
  x <- x[!grepl("^N_", x$gene_id), ]
  
  # Make sure gene order is identical
  if (!identical(x$gene_id, gene_id)) {
    x <- x[match(gene_id, x$gene_id), ]
  }
  
  rna.counts[, i] <- x$unstranded
  
  rm(x)
}

dim(rna.counts)

sum(is.na(rna.counts))

length(unique(rownames(rna.counts)))

saveRDS(
  rna.counts,
  file = "rna_counts_1095_patients.rds",
  compress = FALSE
)

gene.annotation <- data.frame(
  gene_id = gene_id,
  gene_name = gene_name,
  stringsAsFactors = FALSE
)

saveRDS(
  gene.annotation,
  file = "na_gene_annotation.rds",
  compress = FALSE
)


saveRDS(
  maf,
  file = "TCGA_BRCA_WES_MAF.rds",
  compress = FALSE
)

# Save the protein-altering subset
saveRDS(
  maf2,
  file = "TCGA_BRCA_WES_protein_altering_MAF.rds",
  compress = FALSE
)




## download TCGA clinical supplement data
library(TCGAbiolinks)

query.brca.clin <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)

GDCdownload(query.brca.clin)

clinical.BCR <- GDCprepare(query.brca.clin)

names(clinical.BCR)


brca.patient <- clinical.BCR$clinical_patient_brca

dim(brca.patient)

colnames(brca.patient)

brca.patient.use <- brca.patient[ -c(1:2),]

brca.patient.final.cols <- brca.patient.use[, c(
  "bcr_patient_barcode",
  "age_at_diagnosis",
  "er_status_by_ihc",
  "pr_status_by_ihc",
  "her2_status_by_ihc",
  "ajcc_pathologic_tumor_stage",
  "ajcc_tumor_pathologic_pt",
  "ajcc_nodes_pathologic_pn",
  "ajcc_metastasis_pathologic_pm",
  "vital_status",
  "death_days_to",
  "last_contact_days_to",
  "days_to_patient_progression_free",
  "days_to_tumor_progression"
)]

write.csv(brca.patient.final.cols,file="TCGA_BRCA_subset_of_clinical_supplement_information_all_patients.csv",col.names = T,row.names = F, quote = FALSE)

write.csv(brca.patient.use,file="TCGA_BRCA_all_clinical_supplement_information_patients.csv",col.names = T,row.names = F, quote = FALSE)



## download TCGA BRCA DNA meth data
query.meth <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 27",
  sample.type = "Primary Tumor"
)

GDCdownload(query.meth)

library(sesameData)

data.meth <- GDCprepare(
  query = query.meth
)


#matched.samples <- matchedMetExp(
#  project = "TCGA-BRCA"
#)

#head(matched.samples)
#dim(matched.samples)

coldata.meth <- colData(data.meth)

colnames(coldata.meth)
head(coldata.meth)


library(SummarizedExperiment)

assayNames(data.meth)

head(assay(data.meth)[, 1:5])


meth.mat <- assay(data.meth)

cpg.na.count <- apply(meth.mat, 1, function(x){sum(is.na(x))/length(x)})

no.missing <- cpg.na.count[cpg.na.count == 0]

meth.mat.no.missing <- meth.mat[ names(no.missing), ]

var.per.cpg <- apply(meth.mat.no.missing, 1, var)

summary(var.per.cpg)

cpg.keep <- names(var.per.cpg)[ var.per.cpg > 3.790e-03 ]

meth.samples <- query.meth[[1]][[1]]

#meth.mat.keep <- meth.mat.no.missing[ cpg.keep, ]

stopifnot(colnames(meth.mat.no.missing) == meth.samples$cases)

colnames(meth.mat.no.missing) <- meth.samples$cases.submitter_id


write.csv(meth.mat.no.missing,file="TCGA_BRCA_CpG_from_Human_Meth27K.csv",col.names = T,row.names = T, quote = FALSE)

meth.annotation <- rowData(data.meth)

dim(meth.annotation)
colnames(meth.annotation)

meth.annotation <- as.data.frame(meth.annotation)

head(
  meth.annotation[, c(
    "gene",
    "gene_HGNC",
    "chrm_A",
    "beg_A",
    "probeType"
  )]
)

write.csv(meth.annotation,file="TCGA_BRCA_CpG_from_Human_Meth27K_annotation.csv",col.names = T,row.names = T, quote = FALSE)


query.cnv <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Copy Number Variation",
  data.type = "Copy Number Segment",
  sample.type = "Primary Tumor"
)


query.cnv

data.cnv <- GDCprepare(
  query = query.cnv
)

cnv.dir <- "TCGA_BRCA_CNV"

dir.create(cnv.dir, showWarnings = FALSE)

GDCdownload(
  query.cnv,
  method = "api",
  directory = cnv.dir,
  files.per.chunk = 50
)

data.cnv <- GDCprepare(
  query = query.cnv,
  directory = cnv.dir
)

# Basic coverage
length(unique(data.cnv$Sample))
table(table(data.cnv$Sample))

head(unique(data.cnv$Sample))

head(data.cnv[, c("Chromosome", "Start", "End", "Segment_Mean")])

unique(data.cnv$Chromosome)

# confirm the barcode overlap
cnv.sample <- sapply(
  strsplit(data.cnv$Sample, "-"),
  function(x) paste(x[1:4], collapse = "-")
)

cnv.sample <- unique(cnv.sample)

length(cnv.sample)


# Get gene coordinates
library(biomaRt)

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

#gene.annotation <- getBM(
#  attributes = c(
#    "ensembl_gene_id",
#    "hgnc_symbol",
#    "chromosome_name",
#    "start_position",
#    "end_position"
#  ),
#  filters = "ensembl_gene_id",
#  values = sub("\\..*$", "", colnames(X)),
#  mart = ensembl
#)

"org.Hs.eg.db" %in% rownames(installed.packages())
"TxDb.Hsapiens.UCSC.hg38.knownGene" %in% rownames(installed.packages())


library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(GenomicRanges)

names(query.cnv[[1]][[1]])

head(query.cnv[[1]][[1]])

unique(query.cnv[[1]][[1]]$platform)

cnv.metadata <- query.cnv[[1]][[1]]

cnv.metadata$sample.submitter_id[1:10]


cnv.patient <- unique(query.cnv[[1]][[1]]$cases.submitter_id)

length(cnv.patient)
length(intersect(rownames(expr_mut_count), cnv.patient))

head(intersect(rownames(expr_mut_count), cnv.patient))


# create one CNV profile per patient
cnv.metadata <- query.cnv[[1]][[1]]

# Number of CNV samples per patient
sample.count <- table(cnv.metadata$cases.submitter_id)

summary(sample.count)
table(sample.count)

# identify patients with multiple CNV samples
multi.patient <- names(sample.count[sample.count > 1])

head(
  cnv.metadata[
    cnv.metadata$cases.submitter_id %in% multi.patient,
    c(
      "cases.submitter_id",
      "sample.submitter_id",
      "sample_type",
      "file_name"
    )
  ],
  30
)

# Create a clean CNV sample-selection table
cnv.metadata <- query.cnv[[1]][[1]]

# Keep primary tumour samples
cnv.primary <- cnv.metadata[
  cnv.metadata$sample_type == "Primary Tumor",
]

# Prefer 01A samples
cnv.primary$sample_type_code <- substr(
  cnv.primary$sample.submitter_id,
  14,
  16
)

table(cnv.primary$sample_type_code)


# Create a clean CNV sample-selection table
cnv.metadata <- query.cnv[[1]][[1]]

# Keep primary tumour samples
cnv.primary <- cnv.metadata[
  cnv.metadata$sample_type == "Primary Tumor",
]

# Prefer 01A samples
cnv.primary$sample_type_code <- substr(
  cnv.primary$sample.submitter_id,
  14,
  16
)

table(cnv.primary$sample_type_code)


# Patients with more than one Primary Tumor CNV record
sample.count <- table(cnv.primary$cases.submitter_id)

sum(sample.count > 1)

patients <- unique(cnv.primary$cases.submitter_id)

selected.cnv.metadata <- lapply(
  patients,
  function(patient) {
    
    d <- cnv.primary[
      cnv.primary$cases.submitter_id == patient,
    ]
    
    # Prefer 01A if available
    d.01A <- d[
      d$sample_type_code == "01A",
    ]
    
    if (nrow(d.01A) > 0) {
      d <- d.01A
    }
    
    # If multiple records remain, use the most recently updated
    d$updated_datetime <- as.POSIXct(
      d$updated_datetime,
      format = "%Y-%m-%dT%H:%M:%OS"
    )
    
    d <- d[
      order(
        d$updated_datetime,
        decreasing = TRUE
      ),
    ]
    
    d[1, ]
  }
)

selected.cnv.metadata <- do.call(
  rbind,
  selected.cnv.metadata
)

rownames(selected.cnv.metadata) <- NULL


table(
  selected.cnv.metadata$sample_type_code
)

length(unique(selected.cnv.metadata$cases.submitter_id))

nrow(selected.cnv.metadata)

## check no 01A record for 01B patients
patients_01B <- selected.cnv.metadata$cases.submitter_id[
  selected.cnv.metadata$sample_type_code == "01B"
]

cnv.primary[
  cnv.primary$cases.submitter_id %in% patients_01B,
  c(
    "cases.submitter_id",
    "sample.submitter_id",
    "sample_type_code",
    "file_id",
    "file_name"
  )
]

# check the linkage between the selected metadata and data.cnv
selected.file.ids <- selected.cnv.metadata$file_id

sum(
  selected.file.ids %in% unique(data.cnv$GDC_Aliquot)
)

length(selected.file.ids)

length(unique(data.cnv$GDC_Aliquot))


head(
  selected.cnv.metadata[, c(
    "id",
    "file_id",
    "submitter_id",
    "cases.submitter_id",
    "sample.submitter_id"
  )]
)


head(data.cnv[, c(
  "GDC_Aliquot",
  "Sample"
)])

head(
  unique(data.cnv$Sample)
)

head(
  selected.cnv.metadata$sample.submitter_id
)

data.cnv$sample.submitter_id <- substr(
  data.cnv$Sample,
  1,
  16
)

head(
  data.cnv[, c("Sample", "sample.submitter_id")]
)

sum(
  selected.cnv.metadata$sample.submitter_id %in%
    unique(data.cnv$sample.submitter_id)
)

length(
  unique(selected.cnv.metadata$sample.submitter_id)
)


sample.profile.count <- table(
  data.cnv$sample.submitter_id
)

table(sample.profile.count > 1)


sample.aliquot.count <- tapply(
  data.cnv$Sample,
  data.cnv$sample.submitter_id,
  function(x) length(unique(x))
)

table(sample.aliquot.count)


ambiguous.samples <- names(
  sample.aliquot.count[
    sample.aliquot.count > 1
  ]
)

ambiguous.samples

selected.ambiguous <- selected.cnv.metadata[
  selected.cnv.metadata$sample.submitter_id %in%
    ambiguous.samples,
  c(
    "cases.submitter_id",
    "sample.submitter_id",
    "sample_type_code",
    "file_id",
    "file_name"
  )
]

selected.ambiguous

data.cnv[
  data.cnv$sample.submitter_id %in%
    ambiguous.samples,
  c(
    "Sample",
    "sample.submitter_id"
  )
] |>
  unique()


cnv.primary[
  cnv.primary$sample.submitter_id %in%
    c(
      "TCGA-A7-A26E-01A",
      "TCGA-A7-A26J-01A"
    ),
  c(
    "cases.submitter_id",
    "sample.submitter_id",
    "sample_type_code",
    "file_id",
    "file_name",
    "created_datetime",
    "updated_datetime",
    "id"
  )
]

cnv.primary$updated_datetime <- as.POSIXct(
  cnv.primary$updated_datetime,
  format = "%Y-%m-%dT%H:%M:%OS"
)

cnv.primary[
  cnv.primary$sample.submitter_id %in%
    c(
      "TCGA-A7-A26E-01A",
      "TCGA-A7-A26J-01A"
    ),
  c(
    "cases.submitter_id",
    "sample.submitter_id",
    "file_id",
    "file_name",
    "updated_datetime"
  )
]


# Check whether GDC_Aliquot corresponds to metadata id
sum(
  selected.cnv.metadata$id %in%
    unique(data.cnv$GDC_Aliquot)
)

head(
  selected.cnv.metadata$id
)

head(
  data.cnv$GDC_Aliquot
)

selected.cnv.metadata[
  selected.cnv.metadata$cases.submitter_id %in%
    c("TCGA-A7-A26E", "TCGA-A7-A26J"),
  c(
    "cases.submitter_id",
    "sample.submitter_id",
    "id",
    "file_id",
    "file_name"
  )
]

unique(
  data.cnv[
    data.cnv$sample.submitter_id %in%
      c(
        "TCGA-A7-A26E-01A",
        "TCGA-A7-A26J-01A"
      ),
    c(
      "GDC_Aliquot",
      "Sample",
      "sample.submitter_id"
    )
  ]
)


# map the selected metadata to the full CNV aliquot
cnv.primary[
  cnv.primary$sample.submitter_id %in%
    c(
      "TCGA-A7-A26E-01A",
      "TCGA-A7-A26J-01A"
    ),
  c(
    "cases.submitter_id",
    "sample.submitter_id",
    "file_id",
    "file_name",
    "id"
  )
]

# confirm that all 1,096 selected files exist in download directory
selected.file.names <- selected.cnv.metadata$file_name

downloaded.files <- list.files(
  cnv.dir,
  recursive = TRUE,
  full.names = TRUE
)

downloaded.file.names <- basename(downloaded.files)

sum(
  selected.file.names %in% downloaded.file.names
)

length(selected.file.names)

sum(duplicated(selected.file.names))


# get the selected file paths
selected.file.paths <- downloaded.files[
  basename(downloaded.files) %in%
    selected.cnv.metadata$file_name
]

length(selected.file.paths)

# inspect one selected CNV file
test.cnv <- read.table(
  selected.file.paths[1],
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

head(test.cnv)

c(
  org.Hs.eg.db =
    "org.Hs.eg.db" %in% rownames(installed.packages()),
  
  TxDb.hg38 =
    "TxDb.Hsapiens.UCSC.hg38.knownGene" %in%
    rownames(installed.packages())
)

rownames(installed.packages())[
  grepl(
    "EnsDb|TxDb.Hsapiens|org.Hs.eg",
    rownames(installed.packages())
  )
]

library(EnsDb.Hsapiens.v86)
library(ensembldb)
library(GenomicRanges)

gene.annotation <- genes(
  EnsDb.Hsapiens.v86,
  columns = c(
    "gene_id",
    "gene_name",
    "seq_name",
    "start",
    "end"
  )
)

head(
  as.data.frame(gene.annotation)
)

gene.annotation.df <- as.data.frame(
  gene.annotation
)

gene.annotation.df <- gene.annotation.df[
  ,
  c(
    "gene_id",
    "gene_name",
    "seqnames",
    "start",
    "end"
  )
]

gene.annotation.df$midpoint <- floor(
  (
    gene.annotation.df$start +
      gene.annotation.df$end
  ) / 2
)

gene.annotation.df <- gene.annotation.df[
  gene.annotation.df$seqnames %in%
    c(as.character(1:22), "X", "Y"),
]

colnames(gene.annotation.df)[
  colnames(gene.annotation.df) == "seqnames"
] <- "Chromosome"

# Add the base Ensembl ID
gene.annotation.df$gene_id_base <- sub(
  "\\..*$",
  "",
  gene.annotation.df$gene_id
)


#build the CNV gene matrix
cnv.base.dir <- file.path(
  cnv.dir,
  "TCGA-BRCA",
  "Copy_Number_Variation",
  "Copy_Number_Segment"
)

selected.file.paths <- file.path(
  cnv.base.dir,
  selected.cnv.metadata$file_id,
  selected.cnv.metadata$file_name
)

sum(file.exists(selected.file.paths))


all(file.exists(selected.file.paths))

# Verify the UUID → filename pairing
file.exists(
  file.path(
    cnv.base.dir,
    selected.cnv.metadata$file_id,
    selected.cnv.metadata$file_name
  )
) |> table()



# gene-level conversion
# Define the midpoint mapping function
library(GenomicRanges)



map.cnv.to.genes <- function(
    cnv.file,
    gene.annotation
) {
  
  cnv <- read.table(
    cnv.file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  cnv.gr <- GRanges(
    seqnames = cnv$Chromosome,
    ranges = IRanges(
      start = cnv$Start,
      end = cnv$End
    )
  )
  
  gene.gr <- GRanges(
    seqnames = gene.annotation$Chromosome,
    ranges = IRanges(
      start = gene.annotation$midpoint,
      end = gene.annotation$midpoint
    )
  )
  
  hits <- findOverlaps(
    gene.gr,
    cnv.gr,
    select = "first"
  )
  
  result <- rep(
    NA_real_,
    length(gene.gr)
  )
  
  hit.idx <- !is.na(hits)
  
  result[hit.idx] <- cnv$Segment_Mean[
    hits[hit.idx]
  ]
  
  result
}


test.cnv.genes <- map.cnv.to.genes(
  selected.file.paths[1],
  gene.annotation.df
)


sum(is.na(test.cnv.genes))

summary(test.cnv.genes)

sum(!is.na(test.cnv.genes)) /
  length(test.cnv.genes)

# Check whether any gene midpoint has >1 overlapping CNV segment
cnv.matrix <- matrix(
  NA_real_,
  nrow = nrow(gene.annotation.df),
  ncol = length(selected.file.paths)
)

rownames(cnv.matrix) <- gene.annotation.df$gene_id

colnames(cnv.matrix) <-
  selected.cnv.metadata$cases.submitter_id

for (i in seq_along(selected.file.paths)) {
  
  if (i %% 50 == 0 || i == 1) {
    cat(
      "Processing",
      i,
      "of",
      length(selected.file.paths),
      "\n"
    )
  }
  
  cnv.matrix[, i] <- map.cnv.to.genes(
    selected.file.paths[i],
    gene.annotation.df
  )
}

cnv.data <- as.data.frame(
  t(cnv.matrix),
  check.names = FALSE
)

## there are 58613 columns for CNV
rna_Seq = read.csv(file="TCGA_BRCA_VSD_transformed_expression_matrix_above_75th_quantile_variance_overlapping_patients.csv", header = T, stringsAsFactors = F)

rownames(rna_Seq) <- rna_Seq$X

rna_Seq <- rna_Seq[, -c(1)]

rna_Seq_genes <- rownames(rna_Seq)

rna_Seq_genes_base <- sub("\\..*$", "", rna_Seq_genes)

overlapping.genes <- intersect(colnames(cnv.data), rna_Seq_genes_base)

cnv.data.rna.genes <- cnv.data[,c(overlapping.genes)]

cnv.data.missing.ness <- apply(cnv.data.rna.genes, 2, function(x){sum(is.na(x))/length(x)})

summary(cnv.data.missing.ness)


cnv.keep <- cnv.data.missing.ness <= 0.05

cnv.data.rna.genes.filtered <- cnv.data.rna.genes[
  ,
  cnv.keep
]


cnv.sample.missingness <- apply(
  cnv.data.rna.genes,
  1,
  function(x) sum(is.na(x)) / length(x)
)

summary(cnv.sample.missingness)

c(
  `>1%`  = sum(cnv.sample.missingness > 0.01),
  `>5%`  = sum(cnv.sample.missingness > 0.05),
  `>10%` = sum(cnv.sample.missingness > 0.10),
  `>20%` = sum(cnv.sample.missingness > 0.20),
  `>50%` = sum(cnv.sample.missingness > 0.50)
)

cnv.sample.keep <- cnv.sample.missingness <= 0.10

cnv.data.rna.genes.qc <- cnv.data.rna.genes[
  cnv.sample.keep,
  ,
  drop = FALSE
]

cnv.variance <- apply(
  cnv.data.rna.genes.qc,
  2,
  var,
  na.rm = TRUE
)

summary(cnv.variance)

cnv.variance <- cnv.variance[
  !is.na(cnv.variance)
]

# remove the zero/NA-variance feature
cnv.data.rna.genes.qc <- cnv.data.rna.genes.qc[
  ,
  !is.na(cnv.variance) & cnv.variance > 0,
  drop = FALSE
]


write.csv(cnv.data.rna.genes.qc,file="TCGA_BRCA_CNV_data_for_genes_used_from_RNA_data.csv",col.names = T,row.names = T, quote = FALSE)


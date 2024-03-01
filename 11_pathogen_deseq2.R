library(DESeq2)
library(vsn)


options(error = function() {traceback(3); q(1)})

CONTROLS <- c('Healthy', 'NPC')
BASE <- '/projects/b1196/ewa_group/serniczek/data/05_pseudobulk/10_pathogen'


sanitize_name <- function(name) {
  name <- gsub(' ', '_', name, fixed = TRUE)
  name <- gsub('*', '', name, fixed = TRUE)
  name <- gsub(';', '_and', name, fixed = TRUE)
  name <- gsub('/', '_', name, fixed = TRUE)
  return(name)
}

simplify_pathogen <- function(name) {
  name <- gsub('Gram', 'G', name, fixed = TRUE)
  name <- gsub('Pseudomonas aeruginosa', 'Ps', name, fixed = TRUE)
  name <- gsub('SARS-CoV-2', 'COVID', name, fixed = TRUE)
  return(name)
}

simplify_name_for_stupid_deseq2_extra_please <- function(name) {
  name <- gsub(' ', '_', name, fixed = TRUE)
  name <- gsub('*', '', name, fixed = TRUE)
  name <- gsub(';', '_and', name, fixed = TRUE)
  name <- gsub('/', '_', name, fixed = TRUE)
  name <- gsub('Gram+', 'Gram_pos', name, fixed = TRUE)
  name <- gsub('Gram-', 'Gram_neg', name, fixed = TRUE)
  return(name)
}

process_cell_type <- function(cell_type_path) {
  cell_type_dir <- dirname(dirname(cell_type_path))
  print(sprintf('Starting %s', basename(cell_type_dir)))
  R.utils::mkdirs(sprintf('%s/%s', cell_type_dir, '_deseq2'))
  
  meta <- read.csv(cell_type_path, row.names = 1)
  
  if (nrow(meta) < 6 || length(unique(meta$sex)) < 2) {
    return()
  }
  
  controls <- intersect(unique(meta$group), CONTROLS)
  groups <- setdiff(unique(meta$group), CONTROLS)
  comparisons <- expand.grid(controls, groups, stringsAsFactors = FALSE)
  colnames(comparisons) <- c('control', 'condition')
  
  if (nrow(comparisons) == 0) {
    return()
  }
  
  counts <- read.csv(sub('meta.csv', 'count.txt', cell_type_path), sep='\t', check.names = FALSE, row.names = 1)
  
  # safe_names <- sapply(unique(meta$group), simplify_name_for_stupid_deseq2_extra_please)
  # meta$safe_group <- safe_names[meta$group]
  
  dds <- DESeqDataSetFromMatrix(counts, colData = meta, design = ~ group + sex)
  dds <- DESeq(dds, fitType = "local")
  
  transformed <- assay(vst(dds, blind = FALSE))
  write.table(transformed, sprintf('%s/%s/%s', cell_type_dir, '_deseq2', 'transformed.tsv'))
  
  pdf(sprintf('%s/%s/%s', cell_type_dir, '_deseq2', 'disp-local.pdf'), width = 6, height = 4)
  plotDispEsts(dds)
  dev.off()
  
  for (i in 1:nrow(comparisons)) {
    comparison <- comparisons[i, ]
    comparison_name <- sprintf('%s_vs_%s', comparison$control, simplify_pathogen(comparison$condition))
    comparison_dir <- sprintf('%s/%s', cell_type_dir, sanitize_name(comparison_name))
    
    # Require at least 3 samples on each side for strict testing
    if (sum(meta$group == comparison$control) < 3) {
      next
    }
    if (sum(meta$group == comparison$condition) < 3) {
      next
    }
    
    R.utils::mkdirs(comparison_dir)
    
    deseq_meta <- data.frame(control = comparison$control, condition = comparison$condition)
    write.csv(deseq_meta, sprintf('%s/%s', comparison_dir, '_deseq2.csv'))
    
    # print(comparison)
    degs <- as.data.frame(results(
      dds, 
      contrast = c('group', comparison$control, comparison$condition), 
      alpha = 0.05
    ))
    write.csv(degs, sprintf('%s/%s', comparison_dir, 'degs.csv'))
    cat('.', sep='')
  }
}


metas <- list.files(path = BASE, full.names = TRUE, recursive = TRUE, pattern = 'meta.csv')
for (meta_path in metas[3:length(metas)]) {
  process_cell_type(meta_path)
}

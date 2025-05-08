# Simulate enriched scRNA-seq data with CNAs
create_simulated_10x <- function(out_dir = "sim_data",
                                 num_cells = 100,
                                 num_cnv_cells = 20,
                                 ref_group_name = "normal",
                                 cnv_chr = "chr7",
                                 cnv_type = "gain",
                                 cnv_start_gene = 20,
                                 cnv_end_gene = 80) {
  
  dir.create(out_dir, showWarnings = FALSE)
  
  # --- 1. Gene metadata ---
  num_genes <- 2000
  gene_names <- paste0("Gene", 1:num_genes)
  chr <- rep("chr1", num_genes)
  chr[cnv_start_gene:cnv_end_gene] <- cnv_chr
  start <- seq(1, by = 1000, length.out = num_genes)
  end <- start + 999  # assume 1kb gene length
  
  gene_order <- data.frame(
    gene_name = gene_names,
    chromosome = chr,
    start = start,
    end = end
  )
  
  write.table(gene_order,
              file = file.path(out_dir, "gene_order_file.txt"),
              sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  
  # --- 2. Expression matrix ---
  cell_names <- paste0("Cell", 1:num_cells)
  counts <- matrix(rpois(num_genes * num_cells, lambda = 5),
                   nrow = num_genes, ncol = num_cells)
  rownames(counts) <- gene_names
  colnames(counts) <- cell_names
  
  # Apply CNA
  cnv_cells <- sample(cell_names, num_cnv_cells)
  if (cnv_type == "gain") {
    counts[cnv_start_gene:cnv_end_gene, cnv_cells] <- 
      counts[cnv_start_gene:cnv_end_gene, cnv_cells] * 2
  } else if (cnv_type == "loss") {
    counts[cnv_start_gene:cnv_end_gene, cnv_cells] <- 
      round(counts[cnv_start_gene:cnv_end_gene, cnv_cells] * 0.5)
  }
  
  write.table(counts,
              file = file.path(out_dir, "counts.matrix"),
              sep = "\t", quote = FALSE)
  
  # --- 3. Cell metadata ---
  cell_type <- ifelse(cell_names %in% cnv_cells, "cnv", ref_group_name)
  simulated_cnvs <- ifelse(cell_type == "cnv", cnv_type, "none")
  n_counts <- colSums(counts)
  n_genes <- colSums(counts > 0)
  
  annotation <- data.frame(
    cell_name = cell_names,
    cell_type = cell_type,
    simulated_cnvs = simulated_cnvs,
    n_genes = n_genes,
    n_counts = n_counts
  )
  
  write.table(annotation,
              file = file.path(out_dir, "annotation.txt"),
              sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  
  cat("✔ Simulation complete:", out_dir, "\n")
}


create_simulated_10x(
  out_dir = "sim_large_gain_highfreq",
  num_cells = 500,
  num_cnv_cells = 250,
  cnv_chr = "chr7",
  cnv_type = "gain",
  cnv_start_gene = 20,
  cnv_end_gene = 120
)


# Simulate a small CNA (few genes) affecting a small number of cells
create_simulated_10x(
  out_dir = "sim_small_gain_lowfreq",
  num_cells = 500,             # Total number of cells
  num_cnv_cells = 25,          # Only 5% have the CNA
  cnv_chr = "chr8",            # Simulated CNA on chr8
  cnv_type = "gain",           # Gain = overexpression
  cnv_start_gene = 100,        # Small region (e.g., 10 genes)
  cnv_end_gene = 110
)


# Simulate a medium CNA affecting ~20% of the population
create_simulated_10x(
  out_dir = "sim_medium_gain_midfreq",
  num_cells = 500,             # Total number of cells
  num_cnv_cells = 100,         # 100/500 = 20% affected
  cnv_chr = "chr12",           # Chromosome with CNA
  cnv_type = "gain",           # Gain (increased expression)
  cnv_start_gene = 300,        # Medium region (e.g., 60 genes)
  cnv_end_gene = 360
)

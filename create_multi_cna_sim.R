create_multi_cna_sim1 <- function(out_dir = "sim_mixed_cnas", num_cells = 500) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  # --- 1. Gene metadata ---
  num_genes <- 2000
  gene_names <- paste0("Gene", 1:num_genes)
  start <- seq(1, by = 1000, length.out = num_genes)
  end <- start + 999
  chr <- rep("chr1", num_genes)
  chr[501:1000] <- "chr5"
  chr[1001:1500] <- "chr9"
  
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
  
  # --- 3. Assign CNA groups ---
  set.seed(42)
  cnv_small_cells <- sample(cell_names, 25)  # 5%
  remaining_cells <- setdiff(cell_names, cnv_small_cells)
  cnv_medium_cells <- sample(remaining_cells, 100)  # 20%
  remaining_cells <- setdiff(remaining_cells, cnv_medium_cells)
  cnv_large_cells <- sample(remaining_cells, 250)  # 50%
  
  # --- 4. Apply CNAs ---
  counts[50:55, cnv_small_cells] <- counts[50:55, cnv_small_cells] * 2
  counts[600:640, cnv_medium_cells] <- counts[600:640, cnv_medium_cells] * 2
  counts[1100:1250, cnv_large_cells] <- counts[1100:1250, cnv_large_cells] * 2
  
  write.table(counts,
              file = file.path(out_dir, "counts.matrix"),
              sep = "\t", quote = FALSE,
              row.names = TRUE, col.names = TRUE)
  
  # --- 5. Cell metadata ---
  cell_type <- rep("normal", num_cells)
  names(cell_type) <- cell_names
  cell_type[cnv_small_cells] <- "small_gain"
  cell_type[cnv_medium_cells] <- "medium_gain"
  cell_type[cnv_large_cells] <- "large_gain"
  
  cell_type <- as.character(cell_type)  # avoid factor issues
  
  simulated_cnvs <- ifelse(cell_type == "normal", "none", "gain")  # binary label
  
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
  
  cat("✔ Mixed CNA simulation complete:", out_dir, "\n")
}


create_multi_cna_sim("sim_mixed_cna_goldstandard", num_cells = 500)


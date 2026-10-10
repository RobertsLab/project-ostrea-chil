# Fill the network-dependent cache for 15-GN020540-locus.Rmd from a node with
# internet access (klone login node), inside the lab R container, from code/:
#   apptainer exec --no-mount bind-paths -B /mmfs1 -B /etc/resolv.conf -B /etc/hosts <srlab R container> Rscript 15-GN020540-fetch.R
# Runs only the setup, gene, expression and haplotype-download chunks of the notebook.
library(knitr)
src <- purl("15-GN020540-locus.Rmd", output = tempfile(fileext = ".R"), quiet = TRUE,
            documentation = 1)
chunks <- split(readLines(src), cumsum(grepl("^## ----", readLines(src))))
keep <- c("packages", "parameters", "gene", "expression-counts", "expression-bams",
          "haplotypes-download")
for (ch in chunks) {
  label <- sub("^## ----(.*?)(,.*|-{2,}.*)$", "\\1", ch[1], perl = TRUE)
  if (label %in% keep) eval(parse(text = ch), envir = globalenv())
}
message("Cached: ", paste(list.files(out_dir, pattern = "09_gene_count|rnaseq_reads"), collapse = ", "))

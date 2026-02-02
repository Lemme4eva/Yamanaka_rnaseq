# ==============================================================================
# Extraktion der Top-Gene für proximale Tubuli (PTS1, PTS2, PTS3)
# aus dem Chen 2021 Mäusenieren scRNA-seq Datensatz
# ==============================================================================

# 1. DATEN LADEN
# ==============================================================================

# CSV-Datei einlesen
expr_data <- read.csv("data/reference_data/Chen2021_Mouse_Renal_Tubule_RNA-seq_Differential_Genes.csv",
                      row.names = 1,
                      check.names = FALSE)

# Kurzer Überblick
dim(expr_data)
head(expr_data[, 1:10])

# 2. SPALTEN IDENTIFIZIEREN
# ==============================================================================

# Alle Spaltennamen anschauen
colnames(expr_data)

# Spalten für die proximalen Tubuli (PTS1, PTS2, PTS3) identifizieren
pt_pattern <- c("PTS1", "PTS2", "PTS3")
pt_columns <- colnames(expr_data)[grepl(paste(pt_pattern, collapse = "|"),
                                        colnames(expr_data))]

cat("Proximale Tubuli Spalten:\n")
print(pt_columns)
cat("Anzahl PT-Spalten:", length(pt_columns), "\n\n")

# Alle anderen Segmente (zum Vergleich)
other_segments <- c("DTL", "ATL", "MTAL", "CTAL", "DCT", "CNT", "CCD", "OMCD", "IMCD")
other_columns <- colnames(expr_data)[grepl(paste(other_segments, collapse = "|"),
                                           colnames(expr_data))]

# 3. DURCHSCHNITTLICHE EXPRESSION BERECHNEN
# ==============================================================================

# Durchschnitt über alle PT-Segmente
expr_data$PT_mean <- rowMeans(expr_data[, pt_columns], na.rm = TRUE)

# Durchschnitt über alle anderen Segmente (für Spezifizität)
expr_data$Other_mean <- rowMeans(expr_data[, other_columns], na.rm = TRUE)

# Fold-change PT vs. andere Segmente (mit Pseudocount)
expr_data$PT_vs_Other_FC <- log2((expr_data$PT_mean + 1) / (expr_data$Other_mean + 1))

# 4. TOP-GENE IDENTIFIZIEREN
# ==============================================================================

# Sortieren nach PT-Expression
expr_sorted <- expr_data[order(-expr_data$PT_mean), ]

# Top 100 Gene basierend auf PT-Expression
top_100_PT <- head(expr_sorted, 100)

# Top 50 Gene
top_50_PT <- head(expr_sorted, 50)

# Optionale Filterung: Gene, die spezifisch für PT sind
# (höhere Expression in PT als in anderen Segmenten)
pt_specific <- expr_data[expr_data$PT_vs_Other_FC > 0, ]  # FC > 0 bedeutet PT > Other
pt_specific <- pt_specific[pt_specific$PT_mean > 200, ] #keep genes with minimal expression counts of 200 and more
pt_specific_sorted <- pt_specific[order(-pt_specific$PT_vs_Other_FC), ]

top_50_specific <- head(pt_specific_sorted, 50)
top_100_specific <- head(pt_specific_sorted, 100)
top_150_specific <- head(pt_specific_sorted, 150)
top_200_specific <- head(pt_specific_sorted, 200)

# 5. ERGEBNISSE EXPORTIEREN
# ==============================================================================

# Top 50 PT-Gene mit allen Details exportieren
output_top50 <- top_50_PT[, c(pt_columns, "PT_mean", "Other_mean", "PT_vs_Other_FC")]
write.csv(output_top50, "Top50_PT_genes_full_data.csv")

# Nur Gennamen und Kennzahlen
top50_summary <- data.frame(
  Gene_symbol = rownames(top_50_PT),
  PT_mean_expression = top_50_PT$PT_mean,
  Other_segments_mean = top_50_PT$Other_mean,
  PT_vs_Other_logFC = top_50_PT$PT_vs_Other_FC,
  row.names = NULL
)
#write.csv(top50_summary, "Top50_PT_genes_summary.csv", row.names = FALSE)

# Top 50 PT-spezifische Gene
top50_specific_summary <- data.frame(
  Gene_symbol = rownames(top_50_specific),
  PT_mean_expression = top_50_specific$PT_mean,
  Other_segments_mean = top_50_specific$Other_mean,
  PT_vs_Other_logFC = top_50_specific$PT_vs_Other_FC,
  row.names = NULL
)
write.csv(top50_specific_summary, "results/cell-specific analyses/Top50_PT_specific_genes_summary.csv", row.names = FALSE)


# Top 100 PT-spezifische Gene
top100_specific_summary <- data.frame(
  Gene_symbol = rownames(top_100_specific),
  PT_mean_expression = top_100_specific$PT_mean,
  Other_segments_mean = top_100_specific$Other_mean,
  PT_vs_Other_logFC = top_100_specific$PT_vs_Other_FC,
  row.names = NULL
)
write.csv(top100_specific_summary, "results/cell-specific analyses/Top100_PT_specific_genes_summary.csv", row.names = FALSE)

# Top 100 PT-spezifische Gene
top150_specific_summary <- data.frame(
  Gene_symbol = rownames(top_150_specific),
  PT_mean_expression = top_150_specific$PT_mean,
  Other_segments_mean = top_150_specific$Other_mean,
  PT_vs_Other_logFC = top_150_specific$PT_vs_Other_FC,
  row.names = NULL
)
write.csv(top150_specific_summary, "results/cell-specific analyses/Top150_PT_specific_genes_summary.csv", row.names = FALSE)

# Top 100 PT-spezifische Gene
top200_specific_summary <- data.frame(
  Gene_symbol = rownames(top_200_specific),
  PT_mean_expression = top_200_specific$PT_mean,
  Other_segments_mean = top_200_specific$Other_mean,
  PT_vs_Other_logFC = top_200_specific$PT_vs_Other_FC,
  row.names = NULL
)
write.csv(top200_specific_summary, "results/cell-specific analyses/Top200_PT_specific_genes_summary.csv", row.names = FALSE)


# 6. VISUALISIERUNG (optional)
# ==============================================================================

# Boxplot der Top 20 Gene in den verschiedenen PT-Segmenten
top_20_genes <- rownames(top_50_PT)[1:20]

# Für die Visualisierung vorbereiten
top_20_data <- expr_data[top_20_genes, pt_columns]

# Heatmap (benötigt pheatmap oder ggplot2)
if (requireNamespace("pheatmap", quietly = TRUE)) {
  pheatmap::pheatmap(log2(top_20_data + 1),
                     main = "Top 20 PT-Gene in proximalen Tubuli Segmenten",
                     cluster_rows = TRUE,
                     cluster_cols = TRUE)
}

# 7. ZUSAMMENFASSUNG AUSGEBEN
# ==============================================================================

cat("=== ZUSAMMENFASSUNG ===\n")
cat("Datensatz Dimension:", dim(expr_data), "\n")
cat("Gene insgesamt:", nrow(expr_data), "\n")
cat("Proximale Tubuli Spalten:", length(pt_columns), "\n\n")

cat("Top 5 PT-Gene (nach Expressionshöhe):\n")
print(top50_summary[1:5, ])

cat("\n\nTop 5 PT-spezifische Gene (nach PT vs. andere Segmente):\n")
print(top50_specific_summary[1:5, ])

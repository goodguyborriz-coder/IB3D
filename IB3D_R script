# ==============================================================================
# IB3D: Isolation by 3D Resistance (Čistá R verze bez Rcpp a bez alokace RAM)
# ==============================================================================

library(terra)
library(Matrix)

cat("--- SPUŠTĚNÍ ČISTÉHO 3D MODELU IB3D (Iterační low-memory verze) ---\n")

# ==============================================================================
# FÁZE 1: NAČTENÍ A PŘÍPRAVA PROSTOROVÝCH DAT
# ==============================================================================
cat("\n--- FÁZE 1: Příprava prostorových dat ---\n")

path_bathymetry <- "bathymetry.tif"
path_shipping   <- "shipping_intensity.tif"

r_bathy_raw <- rast(path_bathymetry)
r_ships_raw <- rast(path_shipping)

cat("Optimalizuji horizontální rozlišení (fact = 30)...\n")
r_bathy <- aggregate(r_bathy_raw, fact = 30, fun = "mean")
r_ships <- aggregate(r_ships_raw, fact = 30, fun = "mean")

r_ships_resampled <- resample(r_ships, r_bathy, method = "bilinear")

if (min(values(r_bathy), na.rm = TRUE) < 0) {
  r_bathy <- abs(r_bathy)
}

dz     <- 5         
z_min  <- 3         
R_base <- 1         
R_max_barrier <- 1e9 

z0     <- 5         
gamma  <- 1.5       

dim_x <- ncol(r_bathy)
dim_y <- nrow(r_bathy)
max_depth <- max(values(r_bathy), na.rm = TRUE)
dim_z <- ceiling(max_depth / dz)

cat(sprintf("Geometrie 3D mřížky: X = %d, Y = %d, Z = %d vertikálních vrstev\n", dim_x, dim_y, dim_z))

voxel_resistance <- array(R_max_barrier, dim = c(dim_y, dim_x, dim_z))
matrix_bathy <- as.matrix(r_bathy, wide = TRUE)
matrix_ships <- as.matrix(r_ships_resampled, wide = TRUE)

matrix_bathy[is.na(matrix_bathy)] <- 0
matrix_ships[is.na(matrix_ships)] <- 0

cat("Generuji 3D hydroakustické pole odporu...\n")
for (k in 1:dim_z) {
  z_current <- (k - 0.5) * dz
  layer_R <- matrix(R_max_barrier, nrow = dim_y, ncol = dim_x)
  is_water <- (matrix_bathy >= z_current) & (matrix_bathy >= z_min)
  
  if (any(is_water)) {
    attenuated_noise <- matrix_ships / ((1 + (z_current / z0))^gamma)
    layer_R[is_water] <- R_base + attenuated_noise[is_water]
  }
  voxel_resistance[, , k] <- layer_R  
}
cat("Fáze 1 úspěšně dokončena.\n")


# ==============================================================================
# FÁZE 2: LINEARIZACE 3D GRAFU (KIRCHHOFF LAPLACIAN)
# ==============================================================================
cat("\n--- FÁZE 2: Konstrukce Kirchhoffova Laplaciánu L ---\n")
alpha <- 5.0  

valid_mask <- voxel_resistance < 1e8  
nodes_idx  <- which(valid_mask)       
N_nodes    <- length(nodes_idx)

cat(sprintf("Celkový počet aktivních vodních 3D uzlů: %d\n", N_nodes))

node_ids <- array(0, dim = dim(voxel_resistance))
node_ids[nodes_idx] <- 1:N_nodes

subs <- arrayInd(nodes_idx, dim(voxel_resistance))
y <- subs[, 1]; x <- subs[, 2]; z <- subs[, 3]

max_edges <- N_nodes * 3
i_vec <- integer(max_edges); j_vec <- integer(max_edges); v_vec <- numeric(max_edges)
edge_counter <- 0

# Vazby: Směr X+
valid_neighbor <- (x + 1 <= dim_x)
if (any(valid_neighbor)) {
  idx_from <- nodes_idx[valid_neighbor]
  idx_to   <- nodes_idx[valid_neighbor] + dim(voxel_resistance)[1]
  is_ok <- valid_mask[idx_to]
  if (any(is_ok)) {
    id_from <- node_ids[idx_from[is_ok]]; id_to <- node_ids[idx_to[is_ok]]
    R_A <- voxel_resistance[idx_from[is_ok]]; R_B <- voxel_resistance[idx_to[is_ok]]
    R_vazba <- (2 * R_A * R_B) / (R_A + R_B)
    n_new <- length(id_from); rows <- (edge_counter + 1):(edge_counter + n_new)
    i_vec[rows] <- id_from; j_vec[rows] <- id_to; v_vec[rows] <- 1 / R_vazba
    edge_counter <- edge_counter + n_new
  }
}

# Vazby: Směr Y+
valid_neighbor <- (y + 1 <= dim_y)
if (any(valid_neighbor)) {
  idx_from <- nodes_idx[valid_neighbor]
  idx_to   <- nodes_idx[valid_neighbor] + 1
  is_ok <- valid_mask[idx_to]
  if (any(is_ok)) {
    id_from <- node_ids[idx_from[is_ok]]; id_to <- node_ids[idx_to[is_ok]]
    R_A <- voxel_resistance[idx_from[is_ok]]; R_B <- voxel_resistance[idx_to[is_ok]]
    R_vazba <- (2 * R_A * R_B) / (R_A + R_B)
    n_new <- length(id_from); rows <- (edge_counter + 1):(edge_counter + n_new)
    i_vec[rows] <- id_from; j_vec[rows] <- id_to; v_vec[rows] <- 1 / R_vazba
    edge_counter <- edge_counter + n_new
  }
}

# Vazby: Směr Z+ 
valid_neighbor <- (z + 1 <= dim_z)
if (any(valid_neighbor)) {
  idx_from <- nodes_idx[valid_neighbor]
  idx_to   <- nodes_idx[valid_neighbor] + (dim(voxel_resistance)[1] * dim(voxel_resistance)[2])
  is_ok <- valid_mask[idx_to]
  if (any(is_ok)) {
    id_from <- node_ids[idx_from[is_ok]]; id_to <- node_ids[idx_to[is_ok]]
    R_A <- voxel_resistance[idx_from[is_ok]]; R_B <- voxel_resistance[idx_to[is_ok]]
    R_vazba <- ((2 * R_A * R_B) / (R_A + R_B)) * alpha
    n_new <- length(id_from); rows <- (edge_counter + 1):(edge_counter + n_new)
    i_vec[rows] <- id_from; j_vec[rows] <- id_to; v_vec[rows] <- 1 / R_vazba
    edge_counter <- edge_counter + n_new
  }
}

i_vec <- i_vec[1:edge_counter]; j_vec <- j_vec[1:edge_counter]; v_vec <- v_vec[1:edge_counter]

cat("Sestavuji prostorovou matici sousednosti L...\n")
L_odchylna <- sparseMatrix(i = i_vec, j = j_vec, x = -v_vec, dims = c(N_nodes, N_nodes), symmetric = TRUE)
diag_values <- colSums(abs(L_odchylna))
L <- Diagonal(x = diag_values) + L_odchylna
cat("Fáze 2 úspěšně dokončena.\n")


# ==============================================================================
# FÁZE 3: PAMĚŤOVĚ BEZPEČNÝ VÝPOČET PO DVOJICÍCH (ČISTÉ RKO)
# ==============================================================================
# ==============================================================================
# OPTIMALIZOVANÁ FÁZE 3: BLESKOVÝ ITERAČNÍ CG SOLVER PŘES MATICOVÉ NÁSOBENÍ
# ==============================================================================

cat("\n--- Spouštím finální bleskový CG solver (63 jedinců) ---\n")

# 1. Mapování jedinců na uzly (zůstává v paměti, ale pro jistotu)
data_gen <- read.csv("genotypes_with_coordinates.csv")
N_samples <- nrow(data_gen)
sample_node_ids <- integer(N_samples)

ext_bathy <- ext(r_bathy)
res_x <- res(r_bathy)[1]; res_y <- res(r_bathy)[2]

for (i in 1:N_samples) {
  lon <- data_gen$Longitude[i]; lat <- data_gen$Latitude[i]
  col_idx <- max(1, min(dim_x, floor((lon - ext_bathy$xmin) / res_x) + 1))
  row_idx <- max(1, min(dim_y, dim_y - floor((lat - ext_bathy$ymin) / res_y)))
  
  target_node <- 0
  for (k in 1:dim_z) {
    if (valid_mask[row_idx, col_idx, k]) {
      target_node <- node_ids[row_idx, col_idx, k]
      break
    }
  }
  if (target_node == 0) {
    dist_matrix <- which(valid_mask[, , 1], arr.ind = TRUE)
    dists <- (dist_matrix[, 1] - row_idx)^2 + (dist_matrix[, 2] - col_idx)^2
    closest <- dist_matrix[which.min(dists), ]
    target_node <- node_ids[closest[1], closest[2], 1]
  }
  sample_node_ids[i] <- target_node
}

# Inicializace matice
mat_res <- matrix(0, nrow = N_samples, ncol = N_samples)
rownames(mat_res) <- data_gen$SampleID
colnames(mat_res) <- data_gen$SampleID

# Stabilizace matice (jemné uzemnění na diagonále, aby byla matice plně regulární)
L_stable <- L
diag(L_stable) <- diag(L_stable) + 1e-6

# Hlavní smyčka - běží pouze 63krát!
cat("Počítám napěťová pole jedinců přes optimalizovaný CG...\n")

for (i in 1:N_samples) {
  node_src <- sample_node_ids[i]
  
  # Injektujeme jednotkový proud do uzlu daného jedince
  b <- numeric(N_nodes)
  b[node_src] <- 1
  
  # Vysoce optimalizovaný Conjugate Gradient v R s využitím bleskového C-násobení L %*% p
  x_vec <- numeric(N_nodes)
  r <- b - as.numeric(L_stable %*% x_vec)
  p <- r
  rsold <- sum(r * r)
  
  # Max 300 iterací bohatě stačí pro konvergenci napětí
  for (iter in 1:300) {
    Ap <- as.numeric(L_stable %*% p) # TADY JE TA BLESKOVÁ MAGIE (vteřinová operace)
    alpha_cg <- rsold / sum(p * Ap)
    x_vec <- x_vec + alpha_cg * p
    r <- r - alpha_cg * Ap
    rsnew <- sum(r * r)
    if (sqrt(rsnew) < 1e-4) break
    p <- r + (rsnew / rsold) * p
    rsold <- rsnew
  }
  
  # Spočítáme rezistence pro všechny ostatní jedince z tohoto napěťového pole
  for (j in 1:N_samples) {
    node_snk <- sample_node_ids[j]
    R_eff <- x_vec[node_src] - x_vec[node_snk]
    mat_res[i, j] <- max(0, R_eff)
  }
  
  cat(sprintf("Pokrok: Jedinec %d / %d kompletně napočítán\n", i, N_samples))
}

# Symetrizace matice (odstranění drobných numerických nepřesností)
mat_res <- (mat_res + t(mat_res)) / 2

cat("\n==============================================================================\n")
cat("--- HOTOVO! VÝPOČET BYL ÚSPĚŠNĚ A RYCHLE DOKONČEN! ---\n")
cat("Výsledná párová matice 3D rezistence prostředí je v objektu 'mat_res'.\n")
cat("==============================================================================\n")
==============================================================================\n")
cat("--- VÝPOČET KOMPLETNĚ A ÚSPĚŠNĚ DOKONČEN! ---\n")
cat("Výsledná párová matice 3D rezistence prostředí je uložena v objektu 'mat_res'.\n")
cat("==============================================================================\n")
# Bezpečné uložení výsledků
write.csv(mat_res, "IB3D_final_resistance_matrix.csv", row.names = TRUE)
saveRDS(mat_res, "IB3D_final_resistance_matrix.rds")

cat("Výsledky byly bezpečně uloženy do pracovní složky!\n")
library(ggplot2)

# Převod široké matice na dlouhý formát (tidy data) pro ggplot
df_heatmap <- as.data.frame(mat_res)
df_heatmap$SampleID_A <- rownames(df_heatmap)

df_long <- reshape(df_heatmap, 
                   varying = colnames(mat_res), 
                   v.names = "Resistance",
                   timevar = "SampleID_B",
                   times = colnames(mat_res),
                   direction = "long")

# Vykreslení
ggplot(df_long, aes(x = SampleID_A, y = SampleID_B, fill = Resistance)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "3D Odpor") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Heatmapa 3D odporu prostředí mezi jedinci",
    x = "Jedinec A",
    y = "Jedinec B"
  )

# Uložení grafu v tiskové kvalitě
ggsave("IB3D_heatmapa_63_jedincu.png", width = 10, height = 9, dpi = 300)


library(org.Hs.eg.db)

# Definir los directorios desde los cuales cargar las muestras
directorios <- c(
  'C:/Users/NOE/Desktop/TFG-Example/Islet data/Mostres/samples IGT/',
  'C:/Users/NOE/Desktop/TFG-Example/Islet data/Mostres/samples ND/',
  'C:/Users/NOE/Desktop/TFG-Example/Islet data/Mostres/samples T2D/',
  'C:/Users/NOE/Desktop/TFG-Example/Islet data/Mostres/samples T3cD/'
)

# Inicializar una lista para almacenar los datos de todos los archivos
data_list <- list()

# Iterar sobre cada directorio
for (dir in directorios) {
  # Obtener la lista de archivos en el directorio actual
  files <- list.files(dir, full.names = TRUE)
  
  # Iterar sobre cada archivo en el directorio
  for (file in files) {
    # Leer el archivo CSV y agregarlo a la lista de datos
    data <- read.csv(file, sep = "\t")
    
    # Extraer el nombre del archivo sin la extensión "_htseq" o ".txt"
    nombre_archivo <- gsub("\\.txt|_htseq.*", "", basename(file), perl = TRUE)
    
    # Agregar el prefijo correspondiente al directorio
    if (grepl("samples IGT", dir)) {
      nombre_archivo <- paste0("IGT_", nombre_archivo)
    } else if (grepl("samples ND", dir)) {
      nombre_archivo <- paste0("ND_", nombre_archivo)
    } else if (grepl("samples T2D", dir)) {
      nombre_archivo <- paste0("T2D_", nombre_archivo)
    } else if (grepl("samples T3cD", dir)) {
      nombre_archivo <- paste0("T3cD_", nombre_archivo)
    }
    
    
    # Renombrar la columna "counts" con el nombre del archivo
    colnames(data)[colnames(data) == "counts"] <- nombre_archivo
    
    # Establecer la columna "Ensembl" como índice
    rownames(data) <- data$ensembl
    data$ensembl <- NULL
    
    rownames(data) <- gsub("\\.[0-9]*$", "", rownames(data))
    
    data$symbol <- mapIds(org.Hs.eg.db, 
                          keys = rownames(data),
                          keytype = 'ENSEMBL',
                          column = 'SYMBOL')
    
    data <- data[complete.cases(data), ]
    
    data <- data[!duplicated(data$symbol), ]
    
    # Establecer la columna "symbol" como índice
    rownames(data) <- data$symbol
    data$symbol <- NULL
    # Agregar los datos al data_list
    data_list[[paste0(nombre_archivo, "_data")]] <- data
  }
}

for (nombre_archivo in names(data_list)) {
  print(nombre_archivo)
}


combined_data <- NULL

# Iterar sobre cada elemento de la lista data_list y combinar los datos
for (key in names(data_list)) {
  if (is.null(combined_data)) {
    combined_data <- data_list[[key]]
  } else {
    combined_data <- cbind(combined_data, data_list[[key]])
  }
}

# Reordenar las columnas por nombre
combined_data <- combined_data[, order(colnames(combined_data))]

# Mostrar las primeras filas del dataframe combinado
head(combined_data)


)

# Instalar el paquete edgeR desde Bioconductor
BiocManager::install("edgeR")
library(edgeR)

###################################################
'''''''''''''''
levels <- c(rep("IGT", 41), rep("ND", 18), rep("T2D", 39), rep("T3cD", 35))
group <- factor(levels)

dge <- DGEList(counts = combined_data, group=group)
dge <- calcNormFactors(dge)
CPM_all <- cpm(dge)

library(dplyr)
CPM_IGT <- as.data.frame(CPM_all[, 1:41]) %>%
  tibble::rownames_to_column(var = "symbol")

# Para CPM_ND
CPM_ND <- as.data.frame(CPM_all[, 42:59]) %>%
  tibble::rownames_to_column(var = "symbol")

# Para CPM_T2D
CPM_T2D <- as.data.frame(CPM_all[, 60:98]) %>%
  tibble::rownames_to_column(var = "symbol")

# Para CPM_T3cD
CPM_T3cD <- as.data.frame(CPM_all[, 99:133]) %>%
  tibble::rownames_to_column(var = "symbol")

write.table(CPM_IGT, file = "C:/Users/NOE/Desktop/TFG-Example/Islet data/CPM Data/CPM_IGT.txt", sep = "\t", row.names = FALSE)
write.table(CPM_ND, file = "C:/Users/NOE/Desktop/TFG-Example/Islet data/CPM Data/CPM_ND.txt", sep = "\t", row.names = FALSE)
write.table(CPM_T2D, file = "C:/Users/NOE/Desktop/TFG-Example/Islet data/CPM Data/CPM_T2D.txt", sep = "\t", row.names = FALSE)
write.table(CPM_T3cD, file = "C:/Users/NOE/Desktop/TFG-Example/Islet data/CPM Data/CPM_T3cD.txt", sep = "\t", row.names = FALSE)
'''''''''''''''''''''''''''
#################################

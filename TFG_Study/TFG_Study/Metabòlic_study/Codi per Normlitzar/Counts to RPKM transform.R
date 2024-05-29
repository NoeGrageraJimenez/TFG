
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)


##########################################
#TRACTAMENT DE LES MOSTRES:

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
    
    data <- data[complete.cases(data), ] #elimina files amb valor NaN
    
    data <- data[!duplicated(data$symbol), ] #elimina duplicats
    
    # Establecer la columna "symbol" como índice

    data$symbol <- NULL
    # Agregar los datos al data_list
    data_list[[paste0(nombre_archivo, "_data")]] <- data
  }
}


################################
#CREACIÓ DATAFRAME CONJUNT:

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


##################################
#FILTRATGE DE ENSEMBLES EN EL LENGHT MAPPING:
t2g <- read.csv("C:/Users/NOE/Desktop/TFG-Example/Islet data/t2g.csv")

ensembl_of_interest <- rownames(combined_data)
ensembl_of_interest

'''''
#Mirar si estan igual ordenados el vector enseml_of_interest a los indices de combined_data:
igualmente_ordenados <- TRUE

# Comparar cada índice de combined_data con cada valor de ensembl_of_interest
for (i in seq_along(rownames(combined_data))) {
  if (rownames(combined_data)[i] != ensembl_of_interest[i]) {
    igualmente_ordenados <- FALSE
    break
  }
}

# Imprimir el resultado
if (igualmente_ordenados) {
  print("Los índices están igualmente ordenados.")
} else {
  print("Los índices no están igualmente ordenados.")
}
'''''''
t2g_filtrado <- subset(t2g, ensembl_gene_id %in% ensembl_of_interest)
orden <- match(ensembl_of_interest, t2g_filtrado$ensembl_gene_id)


# Ordenar t2g_filtrado utilizando el orden obtenido
t2g_filtrado_ordenado <- t2g_filtrado[orden, ]

'''''''''
#VER SI ESTAN IGUAL DE ORDENADOS LOS ENSEMBL DE T2G FILTERED ORDERED QUE COMBINED DATA:
indices_iguales <- TRUE

# Comparar cada índice de combined_data con cada valor de ensembl_gene_id en t2g_filtrado_ordenado
for (i in seq_along(rownames(combined_data))) {
  if (rownames(combined_data)[i] != t2g_filtrado_ordenado[i, "ensembl_gene_id"]) {
    indices_iguales <- FALSE
    break
  }
}

# Imprimir el resultado
if (indices_iguales) {
  print("Los índices están igualmente ordenados.")
} else {
  print("Los índices no están igualmente ordenados.")
}
'''''''''''''
####################
#Normalitzar les dades a fpkm:

BiocManager::install("edgeR")
library(edgeR)

x <- DGEList(counts=combined_data)
x

genelength <- t2g_filtrado_ordenado$geneLength
genelength
'''''''
#COMPROBAR QUE CADA VALOR DEL VECTOR DE LONGITUDES COINCIDE CON LA MISM POCICION DE T2G_FILTRADO_ORDENADO:
# Crear un vector para almacenar los resultados de la comparación
comparaciones <- vector("logical", length = nrow(t2g_filtrado_ordenado))

# Realizar la comparación
for (i in seq_along(comparaciones)) {
  comparaciones[i] <- t2g_filtrado_ordenado$geneLength[i] == genelength[i]
}

# Verificar si todas las comparaciones son TRUE
if (all(comparaciones)) {
  print("Todos los valores de la columna geneLength coinciden con los del vector geneLength.")
} else {
  print("Al menos un valor de la columna geneLength no coincide con el valor en la misma posición del vector geneLength.")
}
'''''



x$genes$Length <- genelength
x <- calcNormFactors(x)
rpkm <- rpkm(x)
rpkm 

rpkm_df <- data.frame(rpkm)
rpkm_df$symbol <- mapIds(org.Hs.eg.db,
                         keys = rownames(rpkm_df),
                         keytype = 'ENSEMBL',
                         column = 'SYMBOL')
rpkm_df <- rpkm_df[, c("symbol", setdiff(names(rpkm_df), "symbol"))]

#Establecer Symbol como indice:
rownames(rpkm_df) <- rpkm_df$symbol
rpkm_df <- rpkm_df[, -which(names(rpkm_df) == "symbol")]# Eliminar la columna "symbol" del dataframe

#Exportar data frame:
write.table(rpkm_df, file = "C:/Users/NOE/Desktop/TFG-Example/Islet data/RPKM_df.txt", sep = "\t", row.names = TRUE)

###############################################################################
#RPKM <- read.table("C:/Users/NOE/Desktop/TFG-Example/Islet data/RPKM_df.txt", header = TRUE, sep = "\t")

levels <- c(rep("IGT", 41), rep("ND", 18), rep("T2D", 39), rep("T3cD", 35))
group <- factor(levels)


library(dplyr)
RPKM_IGT <- as.data.frame(rpkm_df[, 1:41]) %>%
  tibble::rownames_to_column(var = "symbol")

# Para RPKM_ND
RPKM_ND <- as.data.frame(rpkm_df[, 42:59]) %>%
  tibble::rownames_to_column(var = "symbol")

# Para RPKM_T2D
RPKM_T2D <- as.data.frame(rpkm_df[, 60:98]) %>%
  tibble::rownames_to_column(var = "symbol")

# Para RPKM_T3cD
RPKM_T3cD <- as.data.frame(rpkm_df[, 99:133]) %>%
  tibble::rownames_to_column(var = "symbol")

write.table(RPKM_IGT, file = "C:/Users/NOE/Desktop/TFG-Example/Islet data/RPKM data/RPKM_IGT.txt", sep = "\t", row.names = FALSE)
write.table(RPKM_ND, file = "C:/Users/NOE/Desktop/TFG-Example/Islet data/RPKM data/RPKM_ND.txt", sep = "\t", row.names = FALSE)
write.table(RPKM_T2D, file = "C:/Users/NOE/Desktop/TFG-Example/Islet data/RPKM data/RPKM_T2D.txt", sep = "\t", row.names = FALSE)
write.table(RPKM_T3cD, file = "C:/Users/NOE/Desktop/TFG-Example/Islet data/RPKM data/RPKM_T3cD.txt", sep = "\t", row.names = FALSE)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationDbi")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGen")


library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.KnownGen)
keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)

path='C:/Users/NOE/Desktop/Dadas counts/T3cD/GSM5009263_DP066_htseq_counts.txt/GSM5009263_DP066_htseq_counts.txt';
datos <- read.table(path, header = TRUE, row.names = 1)
datos.df <- as.data.frame(datos)
datos.df
rownames(datos.df)<- gsub("\\.[0-9]*$", "", rownames(datos.df))
datos.df
datos.df$symbol <- mapIds(org.Hs.eg.db, 
                          keys = rownames(datos.df),
                          keytype = 'ENSEMBL',
                          column = 'SYMBOL') 
datos.df

#Buscar cuantos NaN hay:
num_nan <- sum(is.na(datos.df$symbol))
print(num_nan)
print(datos.df)

#Para comprobar que hay ese symbolo
#num_ABCF2 <- sum(datos.df$symbol == "ZCCHC8", na.rm = TRUE)
#print(num_ABCF2)

#Eliminar los NaN:
datos.df <- datos.df[complete.cases(datos.df), ]
print(dim(datos.df)) #Comprobar dimensiones
print(datos.df)

datos.df$ensembl <- rownames(datos.df)
datos.df
rownames(datos.df) <- NULL
head(datos.df)
datos.df <- subset(datos.df, select = -ensembl) #eliminar columna ensembl

datos.df <- datos.df[, c(2, 1)] #canviar las columnas
datos.df

nombre_dir <- "/Users/NOE/Desktop/Dadas counts/T3cD/GSM5009263_s.txt"
nombre_ar <- "GSM5009263_s.txt"
file.create(nombre_ar)

write.table(datos.df, file = paste0(nombre_dir,nombre_ar), sep = "\t", row.names = FALSE, quote=FALSE) #quote quita comillas de los nombres del arxivo creado



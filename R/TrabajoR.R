##########################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       			##################
##########################################################################################


### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
datos <- read.table("archivo-entrada.txt", header = TRUE, sep = "\t")
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)
head(datos)
tail(datos)
dim(datos)

### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación 
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)
datos$media_control <- rowMeans (datos[,2:4])
head (datos[, c("Gene.Name", "media_control")])


datos$media_condicion <- rowMeans (datos[,5:7])
head (datos[, c("Gene.Name", "media_condicion")])


### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).
hist(datos$media_control, main = "HISTOGRAMA MEDIA CONTROL", col = "skyblue", xlab = "media")
hist(datos$media_condicion, main = "HISTOGRAMA MEDIA CONDICION", col = "darkgreen", xlab = "media")


### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)
boxplot(
  list(
    Control = c(datos$C1.TPM, datos$C2.TPM, datos$C3.TPM),
    Condicion = c(datos$X1W5.TPM, datos$X1W6.TPM)
  ),
  col = c("skyblue", "darkgreen"),
  main = "boxplot control y condicion",
  ylab = "Expresión (TPM)",
  xlab = "Grupo",
  border = "black"
)

legend(
  "topleft",
  legend = c("Control", "Condición"),
  fill = c("skyblue", "darkgreen"),
  title = "Grupos",
  bty = "n"
)


### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.
shapiro.test(datos$media_control)
shapiro.test(datos$media_condicion)
## 5.2 Como p < 0.05 en Shapiro (no normal) -> test no paramétrico:

## Wilcoxon apareado (cada gen aparece en Control y en Condición)
datos$p_value <- apply(datos[, 2:7], 1, function(fila) {
  wilcox.test(as.numeric(fila[1:3]),   # Control
              as.numeric(fila[4:6]))$p.value   # Condición
})

head(datos[, c("Gene.Name", "media_control", "media_condicion", "p_value")])

### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)
# Listado de genes diferencialmente expresados
### 6.- Genes diferencialmente expresados (p-value <= 0.1)
deg <- datos[!is.na(datos$p_value) & datos$p_value <= 0.1, ]
deg[, c("Gene.Name", "media_control", "media_condicion", "p_value")]

nrow(deg)

sobre <- deg[deg$media_condicion > deg$media_control, ]
repri <- deg[deg$media_condicion < deg$media_control, ]

cat("Genes DE totales:", nrow(deg), "\n")
# me salen 468 son los P values 0,1 puros t aquellos redondeados a 0,1
cat("Sobreexpresados:", nrow(sobre), "\n")
cat("Reprimidos:", nrow(repri), "\n")

### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)

### 5.- Test de comparación por gen (Wilcoxon, grupos independientes)
shapiro.test(datos$media_control)
shapiro.test(datos$media_condicion)
# p < 0.05 -> no normal -> test no paramétrico

datos$p_value <- apply(datos[, 2:7], 1, function(fila) {
  suppressWarnings(
    wilcox.test(as.numeric(fila[1:3]), as.numeric(fila[4:6]))$p.value
  )
})

### 6.- Genes diferencialmente expresados (p-value <= 0.1)
deg <- datos[!is.na(datos$p_value) & datos$p_value <= 0.1 + 1e-9, ]
deg[, c("Gene.Name", "media_control", "media_condicion", "p_value")]

sobre <- deg[deg$media_condicion > deg$media_control, ]
repri <- deg[deg$media_condicion < deg$media_control, ]

cat("Genes DE totales:", nrow(deg), "\n")
cat("Sobreexpresados:", nrow(sobre), "\n")
cat("Reprimidos:", nrow(repri), "\n")

### 7.- Heatmap de los genes diferencialmente expresados
mat <- as.matrix(deg[, 2:7])
rownames(mat) <- deg$Gene.Name
colnames(mat) <- c("Control1", "Control2", "Control3",
                   "Condicion1", "Condicion2", "Condicion3")

mat_log <- log2(mat + 1)

# Quitar genes sin variación (darían NaN al escalar)
mat_log <- mat_log[apply(mat_log, 1, sd) > 0, , drop = FALSE]

col_columnas <- rep(c("skyblue", "darkgreen"), each = 3)

heatmap(mat_log,
        scale = "row",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        ColSideColors = col_columnas,
        margins = c(8, 8),
        labRow = if (nrow(mat_log) > 60) NA else rownames(mat_log),
        main = "Heatmap de genes diferencialmente expresados (p <= 0.1)")

legend("topright",
       legend = c("Control", "Condición"),
       fill = c("skyblue", "darkgreen"),
       bty = "n")


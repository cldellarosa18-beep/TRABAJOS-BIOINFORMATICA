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

### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)


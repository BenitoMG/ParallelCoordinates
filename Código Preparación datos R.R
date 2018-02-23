####################################################################################
# Autor: 	Benito Martín Gassol

# Introducción: este código ha sido desarrollado como parte del TFM en la
# Universidad Internacional de la Rioja para la visualización de Coordenadas
# Paralelas.

# Objetivo: el objetivo del código es la preparación de datos del fichero de datos que
# va a ser visualizado mediante Coordenadas Paralelas.

# Pasos:
#   1) Lectura del fichero y las variables de interés
#   2) Categorización por Percentiles
#   3) Categorización por Rango
#   4) Cálculo de las medidas para el algoritmo de decisión de orden de las variables
####################################################################################


####################################################################################
# 1) Inicialización
####################################################################################

# Variables a rellenar
nombre_fichero <- "datos.csv"
q <- 10 #Indica el número de intervalos en los que dividir cada variable
dir <- "C:/xampp/htdocs/UNIR_TFM_2017/"
var <- c("Tarifa","Minutos","Antiguedad","Paro","Renta")

# Lectura de datos
datos <- read.table(paste(dir, nombre_fichero,sep=""),sep = ",",header = TRUE)
datos <- apply(datos, 2, function(x) as.numeric(gsub(",",".",x)))
datos <- data.frame(datos)
cases <- nrow(datos)
color <-c("target")
var_color <- c(var,"color")
datos$color <- datos[,color]
datos <- datos[,var_color]
datos_p <- datos
datos_r <- datos
remove(datos)


####################################################################################
# 2) Categorización por percentiles
####################################################################################
for (i in 1:ncol(datos_p[,var])){
  print(var[i])
  print(unique(quantile(datos_p[,i], probs=0:q/q, na.rm = TRUE)))
  datos_p[,i] <- as.integer(cut(datos_p[,i], unique(quantile(datos_p[,i], probs=0:q/q, na.rm = TRUE)), include.lowest=TRUE))
}
datos_p$count <- 1
perc_count <- aggregate(count ~ ., datos_p[,c(var,"count")],FUN=sum)
perc_color <- aggregate(color ~ ., datos_p[,var_color],FUN=mean)
perc <- merge(perc_count, perc_color, by = var)
cutoff <- qpois(0.01,cases/q^length(var))
perc <- perc[perc$count >= cutoff,]
attach(perc)
perc <- perc[order(-Tarifa,-Minutos,-Antiguedad,-Paro,-Renta),]
detach(perc)
write.table(perc$count,paste(dir,"recuento_perc.csv",sep=""),sep=",",row.names = FALSE)
write.table(perc$color,paste(dir,"color_perc.csv",sep=""),sep=",",row.names = FALSE)
write.table(perc[,var],paste(dir,"datos_perc.csv",sep=""),sep=",",row.names = FALSE)


####################################################################################
# 3) Categorización por rango
####################################################################################
for (i in 1:ncol(datos_r[,var])){
  outliers <- quantile(datos_r[,i],probs=c(0.01,0.99), na.rm = FALSE)
  cortes <- seq(outliers[1],outliers[2],(outliers[2]-outliers[1])*1.0/q)
  cortes[1] <- min(datos_r[,i])*0.99
  cortes[q+1] <- max(datos_r[,i])*1.01
  datos_r[,i] <- as.integer(cut(as.numeric(datos_r[,i]),cortes,include.lowest = TRUE))
  remove(outliers)
  remove(cortes)
}
datos_r$count <- 1
rang_count <- aggregate(count ~ ., datos_r[,c(var,"count")],FUN=sum)
rang_color <- aggregate(color ~ ., datos_r[,var_color],FUN=mean)
rang <- merge(rang_count, rang_color, by = var)
cutoff <- qchisq(0.0001,cases/q^length(var))
rang <- rang[rang$count >= cutoff,]
hist(rang$color)
attach(rang)
rang <- rang[order(-Tarifa,-Minutos,-Antiguedad,-Paro,-Renta),]
detach(rang)
write.table(rang$count,paste(dir,"recuento_rang.csv",sep=""),sep=",",row.names = FALSE)
write.table(rang$color,paste(dir,"color_rang.csv",sep=""),sep=",",row.names = FALSE)
write.table(rang[,var],paste(dir,"datos_rang.csv",sep=""),sep=",",row.names = FALSE)


####################################################################################
# 4) Cálculo de la matriz para determinación del orden óptimo
####################################################################################
datos <- read.table(paste(dir, "datos_rang.csv",sep=""),sep = ",",header = TRUE)
datos <- apply(datos, 2, function(x) as.numeric(gsub(",",".",x)))
datos <- data.frame(datos)
cruces <- c()
varianzas <- c()
indice <- 0
for (i in 1:ncol(datos)){
  for (j in 1:ncol(datos)){
    indice <- indice + 1
    auxiliar <- 0
    var_abs <- abs(datos[,i]-datos[,j])
    varianza <- var(var_abs)*(nrow(datos)-1)/nrow(datos)
    for (k in 1:nrow(datos)){
      for (h in 1:nrow(datos)){
        if ((datos[h,i] - datos[k,i])*(datos[h,j]-datos[k,j]) < 0)
        {auxiliar <- auxiliar + 1}
      }
    }
    cruces[indice] <- auxiliar*1.0/2
    varianzas[indice] <- varianza
  }
}
matriz_cruces <- matrix(cruces, nrow = ncol(datos), ncol = ncol(datos))
matriz_varianzas <- matrix(varianzas, nrow = ncol(datos), ncol = ncol(datos))
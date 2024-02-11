# Limpiamos el espacio de trabajo
rm(list=ls())

# Cargamos bibliotecas
library(move)
library(tidyr)
library(ctmm)
library(amt)
library(sp)
library(tidyverse)

# Definimos las rutas donde se encuentran los datos
DATA_PATH <- "/home/milo/Documents/lili/taller_ecologia_movimiento/datos"
DATA_FILE_PATH <- file.path(DATA_PATH, "Movements of Australia's east coast humpback whales.csv")

OUTPUT_PATH <- "/home/milo/Documents/lili/taller_ecologia_movimiento/output"

setwd(OUTPUT_PATH)

# Cargamos los datos
datos_crudos <- read.csv(DATA_FILE_PATH)

# Eliminamos los datos faltantes (NAN)
datos_crudos <- datos_crudos %>% drop_na()

# Removemos registros duplicados
datos_crudos <- datos_crudos[!duplicated(datos_crudos$timestamp),]

# Nos quedamos con un sólo individuo
#datos_crudos <- subset(datos_crudos, datos_crudos$individual.local.identifier == 53348)

# Convertimos estos datos en un objeto MoveStack ya que incluye a varios individuos
datos_move <- move(x=datos_crudos$location.long, y=datos_crudos$location.lat,
                       time=as.POSIXct(datos_crudos$timestamp, format= "%Y-%m-%d %H:%M:%OS", tz="UTC"),
                       proj=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"),
                       data=datos_crudos, animal=datos_crudos$individual.local.identifier)

# Generamos un objeto telemetría
datos_tel <- as.telemetry(datos_move)

# Visualizamos los datos de seguimiento
plot(datos_tel)

#------------- Net squared displacement (NSD) ----------

### Generamos un data frame las coordenadas de seguimiento del individuo X53348 para obtener su NSD
df_coordenadas <- data.frame(x = datos_tel$X53348$longitude, y = datos_tel$X53348$latitude)

### Definimos un objeto make track para ejecutar la instrucción nsd
track_animal <- make_track(df_coordenadas, x,y)

### Ejecutamos la instrucción nsd y guardamos su salida
nsd_animal <- nsd(track_animal)

### Graficamos el NSD
plot(nsd_animal, main="Serie de tiempo de NSD ", xlab="Tiempo", ylab="NSD")

#------------- Kernel Density ----------

# Definimos una función que estima varios kernel density

calcula_kde <- function(animal){
  print(paste0("Ejecución para el animal ",animal@info$identity))
  print("Estimamos el modelo sin considerar la autocorrelacion")
  M.IID <- ctmm.fit(animal)
  print("Estimamos automaticamente los mejores parametros para el modelo utilizando ctmm.guess")
  m.ouf <- ctmm.guess(animal,interactive=FALSE)
  print("Seleccionamos mejor modelo de acuerdo al AIC")
  M.OUF <- ctmm.fit(animal, m.ouf)
  print("Estimamos el akde del modelo iid")
  UD0 <- akde(animal,M.IID)
  print("Estimamos el mejor modelo sin el método de pesos optimos (optimal weighting)")
  UD2 <- akde(animal,M.OUF)
  print("Estimamos el mejor modelo pero utilizando el método de pesos optimos (optimal weighting)")
  UD2w <- akde(animal,M.OUF, weights=TRUE)
  # calculate one extent for all UDs
  EXT <- extent(list(UD0,UD2,UD2w),level=0.95)

  print("Creamos plot IID AKDE")
  cairo_ps(file = paste0("iid_kde_",animal@info$identity,".eps"),onefile = FALSE,fallback_resolution=600)
  plot(animal,UD=UD0,xlim=EXT$x,ylim=EXT$y)
  title(paste0("IID KDE-",animal@info$identity))
  dev.off()

  print("Creamos plot con UD del OUF AKDE")
  cairo_ps(file = paste0("ouf_kde_",animal@info$identity,".eps"),onefile = FALSE,fallback_resolution=600)
  plot(animal,UD=UD2,xlim=EXT$x,ylim=EXT$y)
  title(paste0("OUF AKDE-",animal@info$identity))
  dev.off()

  print("Creamos plot con UD del weighted OUF AKDE")
  cairo_ps(file = paste0("w_ouf_kde_",animal@info$identity,".eps"),onefile = FALSE,fallback_resolution=600)
  plot(animal,UD=UD2w,xlim=EXT$x,ylim=EXT$y)
  title(paste0("weighted OUF AKDE-",animal@info$identity))
  dev.off()

  print("Creamos el raster con IID AKDE")
  animal_ud0 <-raster(UD0)
  writeRaster(animal_ud0, filename=paste0("iid_kde_",animal@info$identity,".tif"), overwrite=TRUE)

  print("Creamos el raster con UD del OUF AKDE")
  animal_ud2 <-raster(UD2)
  writeRaster(animal_ud2, filename=paste0("ouf_kde_",animal@info$identity,".tif"), overwrite=TRUE)

  print("Creamos el raster con UD del weighted OUF AKDE")
  animal_ud2w <-raster(UD2w)
  writeRaster(animal_ud2w, filename=paste0("w_ouf_kde_",animal@info$identity,".tif"), overwrite=TRUE)

}

### Ejecutamos la función. La función guarda los png y tiffs en el directorio OUTPUT_PATH
calcula_kde(datos_tel$X53348)

### Cargamos los tiffs
iid_kde <- raster("iid_kde_53348.tif")
ouf_kde <- raster("ouf_kde_53348.tif")
w_ouf_kde <- raster("w_ouf_kde_53348.tif")

sp::spplot(iid_kde)

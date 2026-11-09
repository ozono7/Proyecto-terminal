# ==========================================================================
# Cinética de crecimiento
# ==========================================================================
# Autor: Dr. José Juan Torres Martínez
# Docente de Ciencias Químicas
# Universidad Autónoma de San Luis Potosí (UASLP)
# México
#
# Material académico para estudiantes
#
# Descripción: Análisis de la cinética de crecimiento celular a partir de datos de densidad óptica.
# ==========================================================================


# Cargar las librerías necesarias
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

# Definir los datos
datos <- data.frame(
  Dias = 0:12,
  # -50%
  M50_1 = c(0.09, 0.102333333, 0.152666667, 0.165666667, 0.195333333, 0.257, 0.313333333, 0.354533333, 0.390333333, 0.433333333, 0.492666667, 0.563333333, 0.591333333),
  M50_2 = c(0.09, 0.106, 0.114666667, 0.174666667, 0.210333333, 0.249666667, 0.319333333, 0.359333333, 0.397666667, 0.444, 0.501666667, 0.572333333, 0.604),
  M50_3 = c(0.09, 0.106333333, 0.157666667, 0.167666667, 0.197666667, 0.249333333, 0.314, 0.354666667, 0.4, 0.442, 0.513333333, 0.568666667, 0.608),
  
  # CONTROL
  CONTROL_1 = c(0.09, 0.091666667, 0.118, 0.142, 0.201, 0.293333333, 0.398666667, 0.448333333, 0.591, 0.698, 0.795, 0.902, 0.905333333),
  CONTROL_2 = c(0.09, 0.106333333, 0.162, 0.195333333, 0.289, 0.376666667, 0.405, 0.497333333, 0.584, 0.699, 0.781333333, 0.867333333, 0.903333333),
  CONTROL_3 = c(0.09, 0.111333333, 0.164666667, 0.186, 0.259666667, 0.300666667, 0.397, 0.482666667, 0.571, 0.686, 0.782333333, 0.901666667, 0.91),
  
  # 50%
  P50_1 = c(0.09, 0.116, 0.185333333, 0.207333333, 0.243, 0.292666667, 0.376666667, 0.439666667, 0.550333333, 0.619666667, 0.689666667, 0.691, 0.688333333),
  P50_2 = c(0.09, 0.108333333, 0.15975, 0.205, 0.249333333, 0.304333333, 0.382, 0.469666667, 0.513, 0.616333333, 0.633666667, 0.684666667, 0.688666667),
  P50_3 = c(0.09, 0.107666667, 0.156666667, 0.219, 0.304, 0.369666667, 0.437, 0.442333333, 0.571333333, 0.615333333, 0.657666667, 0.685333333, 0.685333333),
  
  # 100%
  P100_1 = c(0.09, 0.116, 0.185, 0.199666667, 0.257666667, 0.303666667, 0.372333333, 0.402666667, 0.475, 0.555, 0.595666667, 0.596333333, 0.590666667),
  P100_2 = c(0.09, 0.112, 0.200666667, 0.260666667, 0.308333333, 0.314333333, 0.365666667, 0.402, 0.485333333, 0.526, 0.592333333, 0.597, 0.585333333),
  P100_3 = c(0.09, 0.116333333, 0.154333333, 0.201, 0.267666667, 0.329, 0.372, 0.403333333, 0.474333333, 0.573333333, 0.568666667, 0.596333333, 0.587)
)

# Convertir a formato largo (tidy)
datos_largo <- datos %>%
  pivot_longer(
    cols = -Dias,
    names_to = c("Tratamiento", "Replica"),
    names_sep = "_",
    values_to = "OD"
  ) %>%
  mutate(
    Tratamiento = case_when(
      Tratamiento == "M50" ~ "-50%",
      Tratamiento == "CONTROL" ~ "Control",
      Tratamiento == "P50" ~ "50%",
      Tratamiento == "P100" ~ "100%"
    ),
    Tratamiento = factor(Tratamiento, levels = c("-50%", "Control", "50%", "100%"))
  )

# Calcular medias y desviaciones estándar
datos_sumario <- datos_largo %>%
  group_by(Dias, Tratamiento) %>%
  summarise(
    Media = mean(OD),
    SD = sd(OD),
    .groups = "drop"
  ) %>%
  mutate(
    SD_sup = Media + SD,
    SD_inf = pmax(Media - SD, 0)  # Evitar valores negativos
  )

# Definir colores para cada tratamiento
colores <- c(
  "-50%" = "#4ECDC4",
  "Control" = "#FF6B6B",
  "50%" = "#45B7D1",
  "100%" = "#A593E0"
)

# Crear el gráfico 
ggplot() +
  # Área de desviación estándar (suave)
  geom_ribbon(
    data = datos_sumario,
    aes(x = Dias, ymin = SD_inf, ymax = SD_sup, fill = Tratamiento),
    alpha = 0.15
  ) +
  # Línea de la media (sin puntos)
  geom_line(
    data = datos_sumario,
    aes(x = Dias, y = Media, color = Tratamiento),
    linewidth = 1.2
  ) +
  # Puntos de réplicas (sin leyenda)
  geom_point(
    data = datos_largo,
    aes(x = Dias, y = OD, color = Tratamiento),
    size = 2,
    alpha = 0.6,
    show.legend = FALSE  # Elimina la leyenda de réplicas
  ) +
  
  # Escalas y colores
  scale_color_manual(values = colores) +
  scale_fill_manual(values = colores) +
  
  # Ejes
  scale_x_continuous(breaks = 0:12) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  
  # Etiquetas 
  labs(
    title = "Cinética de crecimiento",
    x = "Tiempo (días)",
    y = "Densidad óptica (OD 600 nm)",
    color = "Tratamiento",
    fill = "Tratamiento"
  ) +
  
  # Tema del grafico
  theme_minimal(base_size = 14) +  # Tamaño base 
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, color = "gray90")
  )

# Guardar el gráfico en alta calidad
ggsave("crecimiento.png", 
       width = 10, height = 6, dpi = 300, bg = "white")

# --------------------------------------------------------------------------
# Interpretación: Cinética de crecimiento
# --------------------------------------------------------------------------
# - Eje X: Tiempo, en días (de 0 a 12 días). Cada punto es un día de cultivo.
# - Eje Y: Densidad óptica (OD) medida a 600 nm. La densidad óptica es una medida de
#   absorbancia: cuánta luz absorbe el cultivo, que crece cuando hay más células.
#   Es adimensional, es decir, no tiene unidades: es un número sin unidad (por eso
#   el eje solo dice OD). El '600 nm' no es la unidad del eje, sino la longitud de
#   onda de la luz con la que se midió la absorbancia.
# - Lo que se aprecia: en los cuatro tratamientos la densidad óptica sube con el
#   tiempo, como una curva de crecimiento: al inicio todas parten de OD 0.09 en el
#   día 0 y luego suben. El Control llega más alto (media de aprox. 0.91 al día 12)
#   y crece más rápido a partir del día 4 o 5. El tratamiento 50% queda en medio
#   (aprox. 0.69 al día 12) y los tratamientos -50% y 100% son los más bajos
#   (aprox. 0.60 y 0.59 al día 12). Cada línea es la media de 3 réplicas, los puntos
#   son las réplicas individuales y la banda de color suave alrededor de la línea es
#   la desviación estándar.
# - Qué significa: todos los cultivos crecieron durante los 12 días, pero el Control
#   creció más y más rápido que los demás tratamientos. Que las curvas se separen a
#   partir del día 4 o 5 indica que el efecto de cada tratamiento sobre el
#   crecimiento se nota con el paso de los días, no desde el inicio.

# --------------------------------------------------------------------------
# Qué aprendimos
# --------------------------------------------------------------------------
# - Qué es la densidad óptica (OD): una medida adimensional de absorbancia que sirve
#   para seguir el crecimiento de un cultivo sin contar células una por una. El
#   número 600 nm indica la longitud de onda de la medición, no una unidad del eje.
# - Cómo leer una cinética de crecimiento: el eje X es el tiempo, el eje Y es la
#   cantidad de biomasa (aquí, OD), y la pendiente de la curva indica qué tan rápido
#   crece el cultivo en cada momento.
# - Que en un gráfico de líneas con réplicas, los puntos muestran los datos reales,
#   la línea muestra la media y la banda de color muestra la variación entre
#   réplicas (desviación estándar).
# - Que la densidad óptica se midió a 600 nm y es adimensional (no tiene unidades),
#   y que el Control creció más que los demás tratamientos (OD de aprox. 0.91 al día
#   12), con las curvas separándose a partir del día 4 o 5.

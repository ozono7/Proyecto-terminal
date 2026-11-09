# ==========================================================================
# Concentración de clorofilas por tratamiento
# ==========================================================================
# Autor: Dr. José Juan Torres Martínez
# Docente de Ciencias Químicas
# Universidad Autónoma de San Luis Potosí (UASLP)
# México
#
# Material académico para estudiantes
#
# Descripción: Análisis estadístico (prueba de Dunnett) y gráfico de la concentración de clorofilas por tratamiento.
# ==========================================================================


# 1. CARGAR PAQUETES Y DATOS
# ----------------------------------
library(tidyverse)
library(multcomp)
library(ggtext)

# Nuevos datos de clorofila calculando la media de las repeticiones
datos_clorofila <- data.frame(
  Grupo = factor(rep(c("-50", "CONTROL", "50%", "100%"), each = 3), levels = c("CONTROL", "-50", "50%", "100%")),
  Repetición = factor(rep(1:3, 4)),
  Clorofila_Total = c(5.205627, 10.661413, 9.143213,
                      5.743707, 6.494213, 2.351240,
                      7.868960, 7.517687, 7.166447,
                      5.298240, 9.244453, 9.668360)
)

# ----------------------------------
# 2. ANÁLISIS ESTADÍSTICO
# ----------------------------------
dunnett_resultados <- datos_clorofila %>%
  group_by(NULL) %>% # No agrupamos por días
  group_modify(~ {
    .x$Grupo <- relevel(.x$Grupo, ref = "CONTROL")
    modelo <- aov(Clorofila_Total ~ Grupo, data = .x)
    dunnett <- glht(modelo, linfct = mcp(Grupo = "Dunnett"))
    resumen <- summary(dunnett)
    
    data.frame(
      contrast = names(coef(dunnett)),
      estimate = as.numeric(coef(dunnett)),
      p.value = as.numeric(resumen$test$pvalues),
      stringsAsFactors = FALSE
    )
  }) %>%
  ungroup() %>%
  mutate(
    contrast = str_remove(contrast, " - CONTROL"),
    p.adj = p.value,
    p.adj.signif = case_when(
      p.adj > 0.05 ~ "ns",
      p.adj <= 0.05 & p.adj > 0.01 ~ "*",
      p.adj <= 0.01 & p.adj > 0.001 ~ "**",
      p.adj <= 0.001 ~ "***"
    )
  )

# Preparar datos para gráfico
datos_resumen_clorofila <- datos_clorofila %>%
  group_by(Grupo) %>%
  summarise(
    Media_Clorofila = mean(Clorofila_Total),
    DE_Clorofila = sd(Clorofila_Total),
    n = n(),
    .groups = "drop"
  )

resultados_grafico <- datos_resumen_clorofila %>%
  filter(Grupo != "CONTROL") %>%
  left_join(dunnett_resultados, by = c("Grupo" = "contrast")) %>%
  mutate(y.position = Media_Clorofila + DE_Clorofila + 0.8)

# ----------------------------------
# 3. GRÁFICO CON LEYENDA
# ----------------------------------
grafico_final_leyenda <- ggplot(datos_resumen_clorofila, aes(x = Grupo, y = Media_Clorofila, fill = Grupo)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = Media_Clorofila - DE_Clorofila,
                    ymax = Media_Clorofila + DE_Clorofila),
                width = 0.25,
                position = position_dodge(0.7),
                size = 0.8) +
  scale_fill_manual(values = c("CONTROL" = "#4ECDC4",
                               "-50" = "#45B7D1",
                               "50%" = "#A593E0",
                               "100%" = "#FF6B6B"),
                    name = "Grupos:") +
  geom_text(data = resultados_grafico %>% filter(!is.na(p.adj.signif)),
            aes(x = Grupo, y = y.position, label = p.adj.signif),
            size = 6, fontface = "bold", color = "black",
            position = position_dodge(width = 0.7)) +
  labs(
    title = "Concentración de Clorofila Total por Tratamiento",
    subtitle = "mg L⁻¹",
    x = "Grupo de Tratamiento",
    y = "Clorofila Total (mg L⁻¹)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
    plot.subtitle = element_text(hjust = 0.5, size = 18),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(20, 20, 20, 20),
    strip.background = element_rect(fill = "lightgray")
  )

# Mostrar gráfico
print(grafico_final_leyenda)

# Exportar en alta calidad
ggsave("clorofila_tratamientos.tiff", plot = grafico_final_leyenda,
       device = "tiff", dpi = 1200, bg = "white",
       width = 10, height = 6, units = "in", compression = "lzw")

# --------------------------------------------------------------------------
# Interpretación: Concentración de Clorofila Total por Tratamiento
# --------------------------------------------------------------------------
# - Eje X: Grupo de Tratamiento (CONTROL, -50, 50% y 100%). Es una variable
#   categórica, sin unidades: cada barra es un tratamiento distinto, no una cantidad
#   que aumente de izquierda a derecha.
# - Eje Y: Clorofila Total, en miligramos por litro (mg L⁻¹). La altura de cada
#   barra es la media de las 3 repeticiones de ese tratamiento. La línea vertical
#   que cruza cada barra es la barra de error: muestra la desviación estándar (DE),
#   es decir, cuánto variaron las repeticiones alrededor de la media. También está
#   en mg L⁻¹, igual que el eje Y.
# - Lo que se aprecia: las cuatro barras están en un rango parecido, entre aprox.
#   4.86 y 8.34 mg L⁻¹. La más baja es el CONTROL (media de aprox. 4.86 mg L⁻¹) y
#   la más alta es el tratamiento -50 (aprox. 8.34 mg L⁻¹), seguida del 100%
#   (aprox. 8.07 mg L⁻¹) y del 50% (aprox. 7.52 mg L⁻¹). Es decir, en promedio los
#   tres tratamientos quedaron por encima del control, pero las barras de error son
#   anchas, sobre todo en el CONTROL y en -50, lo que indica mucha variación entre
#   las repeticiones de un mismo tratamiento.
# - Qué significa: en promedio los tratamientos acumularon más clorofila total que
#   el control, pero la variación entre repeticiones es grande, así que la media
#   sola no basta para concluir. Hay que mirar los símbolos sobre cada barra, que
#   son el resultado de la prueba de Dunnett, que compara cada tratamiento contra el
#   CONTROL (ns = sin diferencia significativa, p > 0.05; * p ≤ 0.05; ** p ≤ 0.01;
#   *** p ≤ 0.001).

# --------------------------------------------------------------------------
# Qué aprendimos
# --------------------------------------------------------------------------
# - Cómo leer un gráfico de barras con media y desviación estándar: la altura de la
#   barra es el promedio de las repeticiones y la barra de error indica cuánto se
#   dispersaron esas repeticiones.
# - Que el eje X puede ser categórico (tratamientos) y no tener unidades, mientras
#   que el eje Y sí las tiene: aquí, miligramos por litro (mg L⁻¹).
# - Que dos barras con medias distintas no siempre son estadísticamente distintas:
#   por eso se usa la prueba de Dunnett para comparar cada tratamiento contra el
#   control antes de sacar conclusiones.
# - Que los tres tratamientos quedaron en promedio por encima del control en clorofila
#   total (el más alto fue -50%, aprox. 8.34 mg L⁻¹, y el control el más bajo, aprox.
#   4.86 mg L⁻¹), y que los símbolos sobre las barras indican si cada tratamiento
#   difiere del control según la prueba de Dunnett.

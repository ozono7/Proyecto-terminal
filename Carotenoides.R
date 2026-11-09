# ==========================================================================
# Concentración de carotenoides por tratamiento
# ==========================================================================
# Autor: Dr. José Juan Torres Martínez
# Docente de Ciencias Químicas
# Universidad Autónoma de San Luis Potosí (UASLP)
# México
#
# Material académico para estudiantes
#
# Descripción: Análisis estadístico (prueba de Dunnett) y gráfico de la concentración de carotenoides por tratamiento.
# ==========================================================================


# 1. CARGAR PAQUETES Y DATOS
# ----------------------------------
library(tidyverse)
library(multcomp)
library(ggtext)

# Nuevos datos de carotenoides calculando la media de las repeticiones
datos_carotenoides <- data.frame(
  Grupo = factor(rep(c("CONTROL", "-50%", "50%", "100%"), each = 3), levels = c("CONTROL", "-50%", "50%", "100%")),
  Repetición = factor(rep(1:3, 4)),
  Carotenoides_Total_original = c(1.607437, 1.293179, -0.050430,
                         2.116677, 0.928509, 4.309613,
                         0.011860, -0.170163, 0.735460,
                         2.969862, 2.429413, -0.455743)
)

# Valores negativos: tres mediciones dieron concentración negativa (-0.050430 en
# CONTROL, -0.170163 en 50% y -0.455743 en 100%), lo cual es físicamente imposible
# para una concentración. Son valores por debajo del límite de detección: al
# restar el blanco, una muestra sin carotenoides detectables puede dar un número
# negativo por el ruido de la medición (Analytical Methods Committee, Royal
# Society of Chemistry, 2001, "What should be done with results below the
# detection limit?", AMC Technical Brief No. 5). Para el análisis esos valores
# se registran como 0 (no detectado); los valores originales se conservan en
# Carotenoides_Total_original.
datos_carotenoides$Carotenoides_Total <- pmax(datos_carotenoides$Carotenoides_Total_original, 0)

# ----------------------------------
# 2. ANÁLISIS ESTADÍSTICO
# ----------------------------------
dunnett_resultados_carotenoides <- datos_carotenoides %>%
  group_by(NULL) %>% # No agrupamos por días
  group_modify(~ {
    .x$Grupo <- relevel(.x$Grupo, ref = "CONTROL")
    modelo <- aov(Carotenoides_Total ~ Grupo, data = .x)
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
datos_resumen_carotenoides <- datos_carotenoides %>%
  group_by(Grupo) %>%
  summarise(
    Media_Carotenoides = mean(Carotenoides_Total),
    DE_Carotenoides = sd(Carotenoides_Total),
    n = n(),
    .groups = "drop"
  )

resultados_grafico_carotenoides <- datos_resumen_carotenoides %>%
  filter(Grupo != "CONTROL") %>%
  left_join(dunnett_resultados_carotenoides, by = c("Grupo" = "contrast")) %>%
  mutate(y.position = Media_Carotenoides + DE_Carotenoides + 0.2) # Ajuste menor para la posición de las etiquetas

# ----------------------------------
# 3. GRÁFICO CON LEYENDA
# ----------------------------------
grafico_final_carotenoides <- ggplot(datos_resumen_carotenoides, aes(x = Grupo, y = Media_Carotenoides, fill = Grupo)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = Media_Carotenoides - DE_Carotenoides,
                    ymax = Media_Carotenoides + DE_Carotenoides),
                width = 0.25,
                position = position_dodge(0.7),
                size = 0.8) +
  scale_fill_manual(values = c("CONTROL" = "#4ECDC4",
                               "-50%" = "#45B7D1",
                               "50%" = "#A593E0",
                               "100%" = "#FF6B6B"),
                    name = "Grupos:") +
  geom_text(data = resultados_grafico_carotenoides %>% filter(!is.na(p.adj.signif)),
            aes(x = Grupo, y = y.position, label = p.adj.signif),
            size = 6, fontface = "bold", color = "black",
            position = position_dodge(width = 0.7)) +
  labs(
    title = "Concentración de Carotenoides por Tratamiento",
    subtitle = "mg L⁻¹",
    x = "Grupo de Tratamiento",
    y = "Carotenoides Total (mg L⁻¹)"
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
print(grafico_final_carotenoides)

# Exportar en alta calidad
ggsave("carotenoides_tratamientos.png", plot = grafico_final_carotenoides,
       dpi = 300, bg = "white",
       width = 10, height = 6, units = "in")

# --------------------------------------------------------------------------
# Interpretación: Concentración de Carotenoides por Tratamiento
# --------------------------------------------------------------------------
# - Eje X: Grupo de Tratamiento (CONTROL, -50%, 50% y 100%). Es una variable
#   categórica, sin unidades: cada barra es un tratamiento distinto, no una cantidad
#   que aumente de izquierda a derecha. Los tratamientos son niveles de un factor
#   (por ejemplo, una dosis o condición aplicada al cultivo).
# - Eje Y: Carotenoides Total, en miligramos por litro (mg L⁻¹). La altura de cada
#   barra es la media de las 3 repeticiones de ese tratamiento. La línea vertical
#   que cruza cada barra es la barra de error: muestra la desviación estándar (DE),
#   es decir, cuánto variaron las repeticiones alrededor de la media. También está
#   en mg L⁻¹, igual que el eje Y.
# - Lo que se aprecia: la barra más alta es la del tratamiento -50% (media de
#   aprox. 2.45 mg L⁻¹) y la más baja es la del 50% (media de aprox. 0.25 mg L⁻¹).
#   El CONTROL queda en aprox. 0.97 mg L⁻¹ y el 100% en aprox. 1.80 mg L⁻¹. Las
#   barras de error son grandes, sobre todo en -50% y en 100%: dentro de un mismo
#   tratamiento las repeticiones dieron valores muy distintos entre sí.
# - Tres mediciones dieron valores negativos (-0.050430 en CONTROL, -0.170163 en
#   50% y -0.455743 en 100%): no son concentraciones reales, son muestras por
#   debajo del límite de detección y en el análisis se registran como 0 (no
#   detectado), con los valores originales conservados en el script. Con ese
#   cambio las medias quedan en aprox. 0.97 (CONTROL), 2.45 (-50%), 0.25 (50%) y
#   1.80 (100%) mg L⁻¹, y la prueba de Dunnett sigue sin encontrar diferencias
#   significativas contra el control.
# - Qué significa: en promedio, el tratamiento -50% acumuló más carotenoides que el
#   control, mientras que el 50% acumuló menos. Pero como la variación entre
#   repeticiones es grande, la media sola no basta para concluir: hay que mirar los
#   símbolos sobre cada barra, que son el resultado de la prueba de Dunnett, que
#   compara cada tratamiento contra el CONTROL (ns = sin diferencia significativa,
#   p > 0.05; * p ≤ 0.05; ** p ≤ 0.01; *** p ≤ 0.001).

# --------------------------------------------------------------------------
# Qué aprendimos
# --------------------------------------------------------------------------
# - Cómo leer un gráfico de barras con media y desviación estándar: la altura de la
#   barra es el promedio de las repeticiones y la barra de error indica cuánto se
#   dispersaron esas repeticiones.
# - Que el eje X puede ser categórico (tratamientos) y no tener unidades, mientras
#   que el eje Y sí las tiene: aquí, miligramos por litro (mg L⁻¹).
# - Que una media alta con una barra de error grande es menos confiable que una media
#   parecida con una barra de error pequeña, y por eso se usa una prueba estadística
#   (Dunnett) para comparar cada tratamiento contra el control.
# - Que una concentración negativa no es un dato real: es una muestra por debajo
#   del límite de detección y en el análisis se registra como 0 (no detectado),
#   sin borrar el valor original.
# - Que el tratamiento -50% tuvo la media más alta de carotenoides (aprox. 2.45 mg L⁻¹)
#   y el 50% la más baja (aprox. 0.25 mg L⁻¹), y que los símbolos sobre las barras
#   indican si cada tratamiento difiere del control según la prueba de Dunnett.

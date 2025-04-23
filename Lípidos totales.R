# 1. CARGAR PAQUETES Y DATOS
# ----------------------------------
library(tidyverse)
library(multcomp)
library(ggtext)

# Datos originales de rendimiento de lípidos
datos_lipidos_raw <- data.frame(
  Tratamiento = c("Control", "-50%", "50%", "100%"),
  rep1 = c(22.8813559, 25.8426966, 15.7480315, 16.5354331),
  rep2 = c(22.8813559, 26.7489712, 16.4285714, 18.7969925),
  rep3 = c(21.7948718, 28.021978, 18.5185185, 16.2162162)
)

# ----------------------------------
# 2. PREPARACIÓN DE DATOS
# ----------------------------------
datos_lipidos_long <- datos_lipidos_raw %>%
  pivot_longer(-Tratamiento,
               names_to = "Repetición",
               values_to = "Rendimiento_lipidos") %>%
  mutate(
    Tratamiento = factor(Tratamiento, 
                         levels = c("Control", "-50%", "50%", "100%"),
                         labels = c("Control", "Reducción 50%", "Aumento 50%", "Aumento 100%"))
  )

# ----------------------------------
# 3. ANÁLISIS ESTADÍSTICO
# ----------------------------------
dunnett_resultados <- datos_lipidos_long %>%
  mutate(Tratamiento = relevel(Tratamiento, ref = "Control")) %>%
  {
    modelo <- aov(Rendimiento_lipidos ~ Tratamiento, data = .)
    dunnett <- glht(modelo, linfct = mcp(Tratamiento = "Dunnett"))
    resumen <- summary(dunnett)
    
    data.frame(
      contrast = names(coef(dunnett)),
      estimate = as.numeric(coef(dunnett)),
      p.value = as.numeric(resumen$test$pvalues),
      stringsAsFactors = FALSE
    )
  } %>%
  mutate(
    contrast = str_remove(contrast, " - Control"),
    p.adj = p.value,
    p.adj.signif = case_when(
      p.adj > 0.05 ~ "ns",
      p.adj <= 0.05 & p.adj > 0.01 ~ "*",
      p.adj <= 0.01 & p.adj > 0.001 ~ "**",
      p.adj <= 0.001 ~ "***"
    )
  )

# Preparar datos para gráfico
datos_resumen_lipidos <- datos_lipidos_long %>%
  group_by(Tratamiento) %>%
  summarise(
    Media_lipidos = mean(Rendimiento_lipidos),
    DE_lipidos = sd(Rendimiento_lipidos),
    n = n(),
    .groups = "drop"
  )

resultados_grafico <- datos_resumen_lipidos %>%
  filter(Tratamiento != "Control") %>%
  left_join(dunnett_resultados, by = c("Tratamiento" = "contrast")) %>%
  mutate(y.position = Media_lipidos + DE_lipidos + 1)

# ----------------------------------
# 4. GRÁFICO CON LEYENDA
# ----------------------------------
grafico_final_leyenda <- ggplot(datos_resumen_lipidos, 
                                aes(x = Tratamiento, y = Media_lipidos, fill = Tratamiento)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = Media_lipidos - DE_lipidos,
                    ymax = Media_lipidos + DE_lipidos),
                width = 0.25,
                position = position_dodge(0.7),
                size = 0.8) +
  scale_fill_manual(values = c("Control" = "#4ECDC4", 
                               "Reducción 50%" = "#45B7D1", 
                               "Aumento 50%" = "#A593E0",
                               "Aumento 100%" = "#FF6B6B"),
                    name = "Tratamientos:") +
  geom_text(data = resultados_grafico %>% filter(!is.na(p.adj.signif)),
            aes(x = Tratamiento, y = y.position, label = p.adj.signif),
            size = 6, fontface = "bold", color = "black",
            position = position_dodge(width = 0.7)) +
  labs(
    title = "Rendimiento de Lípidos Totales",
    subtitle = "Porcentaje (%)",
    x = NULL,  # Eliminamos la etiqueta del eje X
    y = "Rendimiento de Lípidos (%)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
    plot.subtitle = element_text(hjust = 0.5, size = 18),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text.x = element_blank(),  # Ocultamos las etiquetas del eje X
    axis.ticks.x = element_blank(), # Eliminamos las marcas del eje X
    axis.text.y = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# Mostrar gráfico
print(grafico_final_leyenda)

# Exportar en alta calidad
ggsave("rendimiento_lipidos.tiff", plot = grafico_final_leyenda,
       device = "tiff", dpi = 1200, bg = "white",
       width = 8, height = 6, units = "in", compression = "lzw")

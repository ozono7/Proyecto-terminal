# ----------------------------------
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

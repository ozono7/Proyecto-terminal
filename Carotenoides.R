# ----------------------------------
# 1. CARGAR PAQUETES Y DATOS
# ----------------------------------
library(tidyverse)
library(multcomp)
library(ggtext)

# Nuevos datos de carotenoides calculando la media de las repeticiones
datos_carotenoides <- data.frame(
  Grupo = factor(rep(c("CONTROL", "-50%", "50%", "100%"), each = 3), levels = c("CONTROL", "-50%", "50%", "100%")),
  Repetición = factor(rep(1:3, 4)),
  Carotenoides_Total = c(1.607437, 1.293179, -0.050430,
                         2.116677, 0.928509, 4.309613,
                         0.011860, -0.170163, 0.735460,
                         2.969862, 2.429413, -0.455743)
)

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
ggsave("carotenoides_tratamientos.tiff", plot = grafico_final_carotenoides,
       device = "tiff", dpi = 1200, bg = "white",
       width = 10, height = 6, units = "in", compression = "lzw")

# CARGAR PAQUETES
# ----------------------------------
library(tidyverse)
library(multcomp)
library(ggtext)

# ----------------------------------
# 2. PREPARACIÓN DE DATOS
# ----------------------------------
datos_lipidos_raw <- data.frame(
  Concentración = factor(c("CONTROL", "-50%", "50%", "100%"), levels = c("CONTROL", "-50%", "50%", "100%")),
  A_Y1 = c(0.192340172, 1.523091412, 1.832284179, 0.137062769),
  A_Y2 = c(0.466265326, 1.239003919, 2.414543310, 0.127472955),
  A_Y3 = c(1.699830246, 3.072862545, 2.049185725, 0.184582442),
  B_Y1 = c(0.386397923, 3.809560615, 0.012652829, 0.782929559),
  B_Y2 = c(1.998145015, 3.629901898, 0.207998763, 0.519767898),
  B_Y3 = c(0.512868957, 3.503831631, 0.272923233, 0.485387700),
  C_Y1 = c(0.757938791, 1.513072204, 1.127790707, 0.887358336),
  C_Y2 = c(0.240346493, 2.128280217, 0.232817774, 0.554892379),
  C_Y3 = c(0.264163582, 1.111530964, 0.240432372, 1.110071022)
)

# Transformación a formato largo
datos_lipidos_long <- datos_lipidos_raw %>%
  pivot_longer(
    cols = -Concentración,
    names_to = c("Grupo", "Repetición"),
    names_sep = "_",
    values_to = "Lipidos_Neutros"
  ) %>%
  mutate(
    Concentración = factor(Concentración, levels = c("CONTROL", "-50%", "50%", "100%")),
    Grupo = factor(
      case_when(
        Grupo == "A" ~ "Control",
        Grupo == "B" ~ "Tratamiento A",
        Grupo == "C" ~ "Tratamiento B"
      ),
      levels = c("Control", "Tratamiento A", "Tratamiento B")
    )
  )

# ----------------------------------
# 3. ANÁLISIS ESTADÍSTICO (Dunnett)
# ----------------------------------
dunnett_resultados_lipidos <- datos_lipidos_long %>%
  group_by(Concentración) %>%
  group_modify(~ {
    .x$Grupo <- relevel(.x$Grupo, ref = "Control")
    modelo <- aov(Lipidos_Neutros ~ Grupo, data = .x)
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
    contrast = str_remove(contrast, " - Control"),
    p.adj = p.value,
    p.adj.signif = case_when(
      p.adj > 0.05 ~ "ns",
      p.adj <= 0.05 & p.adj > 0.01 ~ "*",
      p.adj <= 0.01 & p.adj > 0.001 ~ "**",
      p.adj <= 0.001 ~ "***"
    )
  )

# Resumen de datos para gráfico
datos_resumen_lipidos <- datos_lipidos_long %>%
  group_by(Concentración, Grupo) %>%
  summarise(
    Media_Lipidos = mean(Lipidos_Neutros),
    DE_Lipidos = sd(Lipidos_Neutros),
    n = n(),
    .groups = "drop"
  )

# Preparación de resultados para gráfico
resultados_grafico_lipidos <- datos_resumen_lipidos %>%
  filter(Grupo != "Control") %>%
  left_join(dunnett_resultados_lipidos, by = c("Concentración", "Grupo" = "contrast")) %>%
  mutate(y.position = Media_Lipidos + DE_Lipidos + 0.1) # Ajuste menor en la posición y

# ----------------------------------
# 4. GRÁFICO FINAL (CON LEYENDA Y SIN ETIQUETAS EN EJE X)
# ----------------------------------
grafico_final_lipidos <- ggplot(datos_resumen_lipidos,
                                aes(x = Grupo, y = Media_Lipidos, fill = Grupo)) +
  geom_col(position = position_dodge(), width = 0.7, color = "black") +
  geom_errorbar(
    aes(ymin = Media_Lipidos - DE_Lipidos,
        ymax = Media_Lipidos + DE_Lipidos),
    width = 0.25,
    position = position_dodge(0.7),
    size = 0.8
  ) +
  facet_grid(. ~ Concentración) +
  scale_fill_manual(
    values = c(
      "Control" = "#4ECDC4",
      "Tratamiento A" = "#45B7D1",
      "Tratamiento B" = "#A593E0"
    ),
    name = "Grupos:"
  ) +
  geom_text(
    data = resultados_grafico_lipidos %>% filter(!is.na(p.adj.signif)),
    aes(y = y.position, label = p.adj.signif),
    size = 6,
    fontface = "bold",
    position = position_dodge(width = 0.7)
  ) +
  labs(
    title = "Concentración de Lípidos Neutros",
    subtitle = "mg/L",
    x = NULL,
    y = "Lípidos Neutros (mg/L)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(family = "Arial"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
    plot.subtitle = element_text(hjust = 0.5, size = 18),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text.x = element_blank(),  # Oculta etiquetas del eje X
    axis.ticks.x = element_blank(), # Elimina marcas del eje X
    axis.text.y = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "cm"),
    strip.text = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(20, 20, 20, 20),
    strip.background = element_rect(fill = "lightgray")
  )

# ----------------------------------
# 5. VISUALIZACIÓN Y EXPORTACIÓN
# ----------------------------------
print(grafico_final_lipidos)

# Exportar en alta calidad (TIFF)
ggsave(
  "lipidos_neutros_final.tiff",
  plot = grafico_final_lipidos,
  device = "tiff",
  dpi = 1200,
  bg = "white",
  width = 10,
  height = 6,
  units = "in",
  compression = "lzw"
)

# ----------------------------------
# 6. LEYENDA DE ANÁLISIS ESTADÍSTICO
# ----------------------------------
cat("
ANÁLISIS ESTADÍSTICO
---------------------------------
Método: ANOVA unidireccional + prueba de Dunnett
Comparación: Cada tratamiento vs Control

NIVELES DE SIGNIFICANCIA:
ns  = no significativo (p > 0.05)
* = p ≤ 0.05
** = p ≤ 0.01
*** = p ≤ 0.001
")

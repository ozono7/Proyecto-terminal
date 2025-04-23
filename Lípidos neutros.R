# 1. CARGAR PAQUETES Y DATOS
# ----------------------------------
library(tidyverse)
library(multcomp)
library(ggtext)

# Datos originales de concentración de lípidos neutros
datos_lipidos_neutros_raw <- data.frame(
  Tratamiento = c("Control", "-50%", "50%", "100%"),
  rep1_1 = c(0.192340172, 1.523091412, 1.832284179, 0.137062769),
  rep1_2 = c(0.466265326, 1.239003919, 2.41454331, 0.127472955),
  rep1_3 = c(0.466265326, 3.072862545, 2.049185725, 0.184582442),
  rep2_1 = c(0.386397923, 3.809560615, 0.012652829, 0.782929559),
  rep2_2 = c(0.386397923, 3.629901898, 0.207998763, 0.519767898),
  rep2_3 = c(0.512868957, 3.503831631, 0.272923233, 0.4853877),
  rep3_1 = c(0.757938791, 1.513072204, 1.127790707, 0.887358336),
  rep3_2 = c(0.240346493, 2.128280217, 0.232817774, 0.554892379),
  rep3_3 = c(0.264163582, 1.111530964, 0.240432372, 1.110071022)
)

# ----------------------------------
# 2. PREPARACIÓN DE DATOS (VERSIÓN CORREGIDA)
# ----------------------------------
datos_lipidos_neutros_long <- datos_lipidos_neutros_raw %>%
  pivot_longer(-Tratamiento,
               names_to = c("Repeticion", "Subrep"),
               names_sep = "_",
               values_to = "Concentracion_lipidos") %>%
  dplyr::select(-Subrep) %>%  # Especificamos el paquete para select
  mutate(
    Tratamiento = factor(Tratamiento, 
                         levels = c("Control", "-50%", "50%", "100%"),
                         labels = c("Control", "Reducción 50%", "Aumento 50%", "Aumento 100%"))
  )

# ----------------------------------
# 3. ANÁLISIS ESTADÍSTICO
# ----------------------------------
dunnett_resultados <- datos_lipidos_neutros_long %>%
  mutate(Tratamiento = relevel(Tratamiento, ref = "Control")) %>%
  {
    modelo <- aov(Concentracion_lipidos ~ Tratamiento, data = .)
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
datos_resumen_lipidos_neutros <- datos_lipidos_neutros_long %>%
  group_by(Tratamiento) %>%
  summarise(
    Media_lipidos = mean(Concentracion_lipidos),
    DE_lipidos = sd(Concentracion_lipidos),
    n = n(),
    .groups = "drop"
  )

resultados_grafico <- datos_resumen_lipidos_neutros %>%
  filter(Tratamiento != "Control") %>%
  left_join(dunnett_resultados, by = c("Tratamiento" = "contrast")) %>%
  mutate(y.position = Media_lipidos + DE_lipidos + 0.5)

# ----------------------------------
# 4. GRÁFICO CON LEYENDA 
# ----------------------------------
grafico_final_leyenda <- ggplot(datos_resumen_lipidos_neutros, 
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
    title = "Concentración de Lípidos Neutros",
    subtitle = "mg/L",
    x = NULL,
    y = "Concentración de Lípidos Neutros (mg/L)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
    plot.subtitle = element_text(hjust = 0.5, size = 18),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
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
ggsave("concentracion_lipidos_neutros.tiff", plot = grafico_final_leyenda,
       device = "tiff", dpi = 1200, bg = "white",
       width = 8, height = 6, units = "in", compression = "lzw")

# ----------------------------------
# 5 LEYENDA DE ANÁLISIS ESTADÍSTICO
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

# Cargar las librerías necesarias
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

# Preparación de datos
datos_biomasa <- data.frame(
  Dias = 0:12,
  B1 = c(0.0009, 0.0038, 0.0041, 0.0032, 0.0052, 0.0046, 0.0054, 0.0067, 0.0081, 0.0091, 0.0101, 0.012, 0.0108),
  B2 = c(0.0032, 0.0037, 0.0046, 0.0052, 0.0055, 0.0062, 0.0066, 0.0077, 0.0084, 0.01, 0.0111, 0.012, 0.0115),
  B3 = c(0.0044, 0.0035, 0.0036, 0.0046, 0.0037, 0.0043, 0.0053, 0.0086, 0.0081, 0.0101, 0.0116, 0.0115, 0.0142),
  OD1 = c(0.1, 0.13, 0.147666667, 0.134333333, 0.149, 0.163666667, 0.175666667,
          0.240666667, 0.236, 0.255666667, 0.269666667, 0.286, 0.287333333),
  OD2 = c(0.1, 0.128, 0.134333333, 0.174333333, 0.202666667, 0.206666667,
          0.222666667, 0.284, 0.278666667, 0.275666667, 0.293, 0.279333333, 0.261666667),
  OD3 = c(0.1, 0.099666667, 0.119, 0.119, 0.134, 0.141333333, 0.178666667,
          0.186, 0.238666667, 0.270666667, 0.315, 0.342, 0.378666667)
)

# Transformar datos a formato largo
datos_correlacion <- datos_biomasa %>%
  pivot_longer(cols = c(B1, B2, B3), names_to = "Biomasa_Rep", values_to = "Biomasa") %>%
  pivot_longer(cols = c(OD1, OD2, OD3), names_to = "OD_Rep", values_to = "OD") %>%
  filter(substr(Biomasa_Rep, 2, 2) == substr(OD_Rep, 3, 3)) %>%
  mutate(
    Réplica = factor(substr(Biomasa_Rep, 2, 2), levels = 1:3, labels = paste0("R", 1:3)),
    Biomasa = as.numeric(Biomasa),
    OD = as.numeric(OD)
  ) %>%
  dplyr::select(Dias, Biomasa, OD, Réplica)

# Calcular el modelo lineal para la correlación
modelo_correlacion <- lm(Biomasa ~ OD, data = datos_correlacion)
summary_modelo <- summary(modelo_correlacion)
r_squared <- summary_modelo$r.squared
ecuacion <- paste0("y = ", round(coef(modelo_correlacion)[2], 4), "x + ", round(coef(modelo_correlacion)[1], 4))

# Encontrar los límites para los ejes con un margen
margen_x <- (max(datos_correlacion$OD) - min(datos_correlacion$OD)) * 0.05
margen_y <- (max(datos_correlacion$Biomasa) - min(datos_correlacion$Biomasa)) * 0.05
limites_x <- c(min(datos_correlacion$OD) - margen_x, max(datos_correlacion$OD) + margen_x)
limites_y <- c(min(datos_correlacion$Biomasa) - margen_y, max(datos_correlacion$Biomasa) + margen_y)

# Gráfico de correlación
grafico <- ggplot(datos_correlacion, aes(x = OD, y = Biomasa, color = Réplica)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, fill = alpha("#FF6B6B", 0.2), color = "#FF6B6B") +
  annotate("text",
           x = min(datos_correlacion$OD) + margen_x,
           y = max(datos_correlacion$Biomasa) - margen_y,
           label = paste("Ecuación:", ecuacion, "\nR² =", round(r_squared, 4)),
           hjust = 0, vjust = 1, size = 4.5, color = "black") +
  scale_color_manual(
    name = "Réplicas",
    values = c("R1" = "#4ECDC4", "R2" = "#45B7D1", "R3" = "#A593E0")
  ) +
  scale_x_continuous(limits = limites_x, expand = c(0, 0)) +
  scale_y_continuous(limits = limites_y, expand = c(0, 0)) +
  labs(
    title = "Correlación entre Densidad Óptica y Biomasa Seca",
    x = "Densidad Óptica (OD600 nm)",
    y = "Biomasa Seca (g/L)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(linewidth = 0.4),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(20, 20, 20, 20)
  )

# Mostrar el gráfico
print(grafico)

# Exportación
ggsave("correlacion_OD_biomasa.png", plot = grafico,
       width = 15, height = 10, units = "cm",
       dpi = 1200, bg = "white")

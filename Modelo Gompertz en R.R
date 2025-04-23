# ----------------------------------------------------
# ANÁLISIS DE CRECIMIENTO MICROBIANO CON MODELO GOMPERTZ
# ----------------------------------------------------

# ====================================================
# SECCIÓN 1: CONFIGURACIÓN INICIAL
# ====================================================

# Instalar y cargar paquetes necesarios
# -----------------------------------
# pacman: Para gestión eficiente de paquetes
# ggplot2: Para visualización gráfica
# tidyr/dplyr: Para manipulación de datos
# minpack.lm: Para ajuste de modelos no lineales
# scales/broom: Para tratamiento de datos y resultados
# ggtext: Para formato de texto en gráficos
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, tidyr, dplyr, minpack.lm, scales, broom, ggtext)

# ====================================================
# SECCIÓN 2: PREPARACIÓN DE DATOS
# ====================================================

# Definición de datos experimentales
# ---------------------------------
# Datos de Densidad Óptica (OD) medidos a 750 nm
# Tres réplicas (OD1, OD2, OD3) durante 13 días (0-12)
datos <- data.frame(
  Dias = 0:12,
  OD1 = c(0.1, 0.13, 0.147666667, 0.134333333, 0.149, 0.163666667, 0.175666667, 
          0.240666667, 0.236, 0.255666667, 0.269666667, 0.286, 0.287333333),
  OD2 = c(0.1, 0.128, 0.134333333, 0.174333333, 0.202666667, 0.206666667, 
          0.222666667, 0.284, 0.278666667, 0.275666667, 0.293, 0.279333333, 0.261666667),
  OD3 = c(0.1, 0.099666667, 0.119, 0.119, 0.134, 0.141333333, 0.178666667, 
          0.186, 0.238666667, 0.270666667, 0.315, 0.342, 0.378666667)
)

# Transformación de datos a formato largo
# --------------------------------------
# pivot_longer: Convierte datos de ancho a largo
# mutate: Calcula horas (24*Dias) y log(N/N0)
# filter: Elimina valores NA si existen
datos_largos <- datos %>%
  pivot_longer(cols = starts_with("OD"), names_to = "Réplica", values_to = "OD") %>%
  mutate(Horas = Dias * 24) %>%
  filter(!is.na(OD))

# Valor inicial de OD (N0)
# ------------------------
# Se toma el valor de OD en tiempo 0 como referencia
OD0 <- 0.1  # Valor común a las tres réplicas en t=0

# Cálculo de log(N/N0)
# --------------------
# log(N/N0) representa el crecimiento logarítmico relativo
datos_largos <- datos_largos %>%
  mutate(log_NN0 = ifelse(OD > 0, log(OD / OD0), NA))

# Resumen estadístico
# -------------------
# Calcula media y desviación estándar para cada tiempo
resumen <- datos_largos %>%
  group_by(Horas) %>%
  summarise(
    Media_logNN0 = mean(log_NN0, na.rm = TRUE),
    DE = sd(log_NN0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(Media_logNN0))

# ====================================================
# SECCIÓN 3: MODELO GOMPERTZ
# ====================================================

# Función del modelo Gompertz
# ---------------------------
# Parámetros:
# t: tiempo
# A: Asíntota (máximo crecimiento posible)
# μ: Tasa de crecimiento inicial
# λ: Tiempo de latencia (retraso antes del crecimiento)
gompertz_logNN0 <- function(t, A, mu, lambda) {
  A * exp(-exp((mu * exp(1) / A) * (lambda - t) + 1))
}

# Ajuste del modelo no lineal
# ---------------------------
# nlsLM: Algoritmo de mínimos cuadrados
# Valores iniciales estimados a partir de los datos
set.seed(123)  # Para reproducibilidad
ajuste <- nlsLM(
  Media_logNN0 ~ gompertz_logNN0(Horas, A, mu, lambda),
  data = resumen,
  start = list(A = max(resumen$Media_logNN0, na.rm = TRUE) * 1.1, 
               mu = 0.005, 
               lambda = 10),
  control = nls.lm.control(maxiter = 500)
)

# Extracción de parámetros
# ------------------------
# params: Coeficientes del modelo ajustado
# conf_int: Intervalos de confianza al 95%
params <- coef(ajuste)
conf_int <- confint(ajuste)

# Cálculo de parámetros derivados
# -------------------------------
# A: Asíntota (capacidad máxima de crecimiento)
# μ: Tasa de crecimiento inicial
# λ: Tiempo de latencia
# μ_max: Tasa máxima de crecimiento (μ_max = (μ*A)/e)
# td: Tiempo de duplicación (td = ln(2)/μ_max)
A <- params["A"]
mu <- params["mu"]
lambda <- params["lambda"]
mu_max <- (mu * A) / exp(1)  # e ≈ 2.71828
td <- log(2) / mu_max

# ====================================================
# SECCIÓN 4: RESULTADOS
# ====================================================

# Mostrar resultados en consola
# -----------------------------
# Formato: Parámetro = valor (límite inferior - límite superior)
cat("\n----------------------------------------\n")
cat("RESULTADOS DEL MODELO GOMPERTZ\n")
cat("----------------------------------------\n")
cat(sprintf("A (asíntota) = %.4f (%.4f - %.4f)\n", A, conf_int["A",1], conf_int["A",2]))
cat(sprintf("μ (tasa crecimiento) = %.4f h⁻¹ (%.4f - %.4f)\n", mu, conf_int["mu",1], conf_int["mu",2]))
cat(sprintf("λ (tiempo latencia) = %.4f h (NA - %.4f)\n", lambda, conf_int["lambda",2]))
cat(sprintf("μ_max (tasa máxima) = %.4f h⁻¹ (%.4f - %.4f)\n", mu_max, 
            (conf_int["mu",1] * A)/exp(1), 
            (conf_int["mu",2] * A)/exp(1)))
cat(sprintf("td (tiempo duplicación) = %.4f h (%.4f - %.4f)\n", td, 
            log(2)/((conf_int["mu",2] * A)/exp(1)), 
            log(2)/((conf_int["mu",1] * A)/exp(1))))
cat("----------------------------------------\n")

# ====================================================
# SECCIÓN 5: VISUALIZACIÓN
# ====================================================

# Texto con parámetros para el gráfico
# ------------------------------------
# Formato similar al de la segunda imagen de referencia
param_text <- paste(
  "μ = ", sprintf("%.3f h-¹", mu), " | ",
  "λ = ", sprintf("%.1f h", lambda), " | ",
  "μ max = ", sprintf("%.1e OD/h", mu_max),
  sep = ""
)

# Construcción del gráfico principal
# ---------------------------------
grafico <- ggplot() +
  # Banda de error (desviación estándar)
  geom_ribbon(
    data = resumen,
    aes(x = Horas, ymin = Media_logNN0 - DE, ymax = Media_logNN0 + DE),
    fill = "#FF6B6B", alpha = 0.15
  ) +
  # Línea de medias experimentales
  geom_line(
    data = resumen,
    aes(x = Horas, y = Media_logNN0),
    color = "#FF6B6B", linewidth = 1, linetype = "dashed"
  ) +
  # Curva de Gompertz ajustada
  geom_line(
    data = data.frame(
      Horas = seq(0, max(resumen$Horas), length.out = 300),
      log_NN0 = gompertz_logNN0(seq(0, max(resumen$Horas), length.out = 300), A, mu, lambda)
    ),
    aes(x = Horas, y = log_NN0),
    color = "#2A363B", linewidth = 1.4
  ) +
  # Puntos experimentales
  geom_point(
    data = datos_largos,
    aes(x = Horas, y = log_NN0, color = Réplica),
    size = 3, alpha = 0.8, na.rm = TRUE
  ) +
  # Eje X (tiempo en horas)
  scale_x_continuous(
    breaks = seq(0, 288, by = 48),
    limits = c(0, 300),
    expand = c(0, 0)
  ) +
  # Eje Y (crecimiento logarítmico)
  scale_y_continuous(
    name = expression(log(frac(N, N[0]))),
    limits = c(0, max(c(resumen$Media_logNN0, datos_largos$log_NN0), na.rm = TRUE) * 1.1)
  ) +
  # Escala de colores para réplicas
  scale_color_manual(
    values = c("#4ECDC4", "#45B7D1", "#A593E0"),
    labels = c("OD1", "OD2", "OD3")
  ) +
  # Etiquetas y títulos
  labs(
    title = "Ajuste del modelo de Gompertz",
    subtitle = param_text,  # Parámetros en subtítulo
    x = "Tiempo (horas)",
    y = "log(N/N₀)",
    color = "Réplica"
  ) +
  # Estilo del gráfico
  theme_minimal(base_size = 13) +
  theme(
    text = element_text(family = "Arial"),
    panel.grid.major = element_line(color = "gray90"),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_markdown(hjust = 0.5),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(20, 20, 20, 20)
  )

# Mostrar gráfico
print(grafico)

# Guardar gráfico en archivo
# --------------------------
# Formato PNG con alta resolución (1,200 dpi)
ggsave("modelo_gompertz_final.tif", 
       plot = grafico,
       width = 18, height = 14, units = "cm",
       dpi = 1200, bg = "white")

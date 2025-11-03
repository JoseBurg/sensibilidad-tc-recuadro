# -------------------------------------------------------------------------
# ANÁLISIS VAR: SENSIBILIDAD DE EXPECTATIVAS CAMBIARIAS A CHOQUES
# -------------------------------------------------------------------------

# --- 0. Carga de Paquetes ---
# Asegurarse de tener instalados estos paquetes (install.packages("nombre"))

library(tidyverse)   # Para manipulación de datos (dplyr) y gráficos (ggplot2)
library(lubridate)   # Para manejo de fechas
library(tseries)     # Para pruebas de estacionariedad (adf.test, kpss.test) [1, 8]
library(urca)        # Para pruebas de raíz unitaria más avanzadas (ur.df, ur.kpss) 
library(vars)        # El paquete central para VAR, VARselect, irf, serial.test [3, 4, 6]

# --- 1. Carga y Preparación de Datos ---

# Asumiendo que los datos están en un dataframe llamado 'data'
# data <- read.csv("su_archivo_de_datos.csv") 

# Asegurar que la fecha esté en formato correcto
data$fecha <- as.Date(data$fecha)

# (Opcional) Convertir a un objeto de serie temporal (ts) para facilitar
# el manejo. Asumiendo que los datos inician en Junio de 2009.
data_ts <- ts(data[, -1], start = c(2009, 6), frequency = 12)


# --- 2. Análisis de Estacionariedad (Pruebas de Raíz Unitaria) ---

# I(1) Pruebas para variables en NIVELES (Sección I.C)
# Se espera p > 0.05 (No Estacionaria)
adf.test(data$exp_tc_media_12)
adf.test(data$ta_90d)
adf.test(data$ta_preferencial)
adf.test(data$dprime)
adf.test(data$ipc)

# Se espera p < 0.05 (Estacionaria) 
adf.test(data$sorpresa_tc_exp) 

# Pruebas KPSS complementarias
# Se espera p < 0.05 (No Estacionaria) [8]
kpss.test(data$exp_tc_media_12, null = "Trend")
kpss.test(data$ta_90d, null = "Trend")
kpss.test(data$ta_preferencial, null = "Trend")
kpss.test(data$dprime, null = "Trend")
kpss.test(data$ipc, null = "Trend")

# Se espera p > 0.10 (Estacionaria)
kpss.test(data$sorpresa_tc_exp, null = "Level")


# --- 3. Transformación de Variables ---

# Aplicar primera diferencia (diff) a las variables I(1) 
d.exp_tc_media_12 <- diff(data$exp_tc_media_12)
d.ta_90d <- diff(data$ta_90d)
d.ta_preferencial <- diff(data$ta_preferencial)
d.dprime <- diff(data$dprime)
d.ipc <- diff(data$ipc)

# La variable I(0) 'sorpresatc_exp' se mantiene en niveles.
# Se debe alinear su longitud, removiendo la primera observación.
sorpresa_tc_aligned <- data$sorpresa_tc_exp[-1]
N <- length(d.ipc) # Longitud final

# Verificación Opcional: Probar que las variables diferenciadas son I(0)
# Se espera p < 0.05 para todos
adf.test(d.exp_tc_media_12)
adf.test(d.ta_90d)
adf.test(d.ta_preferencial)
adf.test(d.dprime) # Esta variable no cumple el criterio 
adf.test(d.ipc)


# --- 4. Construcción del Sistema VAR ---

# Construir el sistema estacionario final usando cbind() 
# El ORDEN es crucial para la descomposición de Cholesky [29, 30]
# Orden: Exógenas -> Endógenas

var_data_stationary <- cbind(
  # d.dprime = d.dprime,
  d.ipc = d.ipc,
  d.ta_90d = d.ta_90d,
  d.ta_preferencial = d.ta_preferencial,
  sorpresatc_exp = sorpresa_tc_aligned,
  d.exp_tc_media_12 = d.exp_tc_media_12
)

# Convertir el conjunto de datos final a un objeto ts
var_data_ts <- ts(var_data_stationary, start = c(2009, 7), frequency = 12)


# --- 5. Especificación del Modelo VAR (Selección de Rezagos) ---

# Usar VARselect para determinar el número óptimo de rezagos (p) [4, 33]
# Se establece un máximo de 12 rezagos (un año)
lag_selection <- VARselect(var_data_ts, lag.max = 12, type = "const")
print(lag_selection$selection) # Muestra el p óptimo para AIC, HQ, SC, FPE [33, 34]
print(lag_selection$criteria)

# Decisión: Seleccionar el rezago (p) basado en los criterios (p.ej., 2)
# p_optimo <- lag_selection$selection # Opción parsimoniosa
p_optimo <- 2 # Asumiendo p=2 como óptimo tras diagnósticos


# --- 6. Estimación y Validación del Modelo VAR ---

# Estimación del modelo VAR(p) [5, 13]
var_model <- VAR(var_data_ts, p = p_optimo + 2)
# 1. Prueba de Autocorrelación Serial (Portmanteau Test) [39, 40, 41]
# H0: No hay autocorrelación serial. Se busca p-valor > 0.05
serial_test <- serial.test(var_model, lags.pt = p_optimo + 5, type = "PT.asymptotic")
print(serial_test)

# 2. Prueba de Estabilidad [43, 45]
# Se busca que todas las raíces (puntos) estén DENTRO del círculo unitario
stability_test <- stability(var_model)
plot(stability_test)

# 3. Prueba de Heterocedasticidad (ARCH) 
# H0: No hay efectos ARCH. Se busca p-valor > 0.05
arch_test <- arch.test(var_model, lags.multi = p_optimo)
print(arch_test)

# Si el modelo falla 'serial_test', se debe re-estimar con p+1 rezagos.


# --- 7. Análisis de Funciones de Impulso-Respuesta (FIR) ---

# Calcular la FIR (Sección III.B) [6, 37]
# Impulso: sorpresatc_exp
# Respuesta: d.exp_tc_media_12
irf_results <- irf(var_model, 
                   impulse = "sorpresatc_exp", 
                   response = "d.exp_tc_media_12", 
                   n.ahead = 24,         # Horizonte de 24 meses
                   ortho = TRUE,         # Usar Cholesky (dependiente del orden)
                   boot = TRUE,          # Calcular bandas de confianza [32]
                   ci = 0.95,            # Nivel de confianza del 95%
                   runs = 1000)          # Iteraciones de Bootstrap

# Visualización de la FIR (Sección III.C) [48, 49]
# Este es el gráfico principal que responde la pregunta
plot(irf_results, 
     ylab = "Respuesta de d.exp_tc_media_12", 
     xlab = "Meses después del choque",
     main = "Respuesta de Expectativas LP ante Choque de Sorpresa Cambiaria")

# (Opcional) Graficar todas las respuestas al choque de sorpresa
irf_all_responses <- irf(var_model, impulse = "sorpresatc_exp", n.ahead = 24, ortho = TRUE, boot = TRUE)
plot(irf_all_responses)

# (Opcional) Graficar todas las respuestas de la expectativa LP
irf_all_impulses <- irf(var_model, response = "d.exp_tc_media_12", n.ahead = 24, ortho = TRUE, boot = TRUE)
plot(irf_all_impulses)

# --- Fin del Script ---
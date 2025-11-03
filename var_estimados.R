library(bsvars)
library(bsvarSIGNs)
library(dplyr)
library(ggplot2)
library(vars)
library(urca)
library(tseries)

datos <- readxl::read_excel("./data-exp-tc.xlsx") |> 
  mutate(sorpresa_tc_exp = tc_mes - exp_tc_media_mes) |> 
  dplyr::select(-c(fecha, exp_tc_media_mes, tc_mes)) |> 
  na.omit()


data_ts <- ts(datos, start = c(2009, 6), frequency = 12)



# I(1) Pruebas para variables en NIVELES --------------
# Se espera p > 0.05 (No Estacionaria)
adf.test(datos$ta_90d)
adf.test(datos$exp_tc_media_12)
adf.test(datos$ta_preferencial)
adf.test(datos$dprime)
adf.test(datos$ipc)

# Se espera p < 0.05 (Estacionaria) 
adf.test(datos$sorpresa_tc_exp) # Estacionaria


# --- 3. Transformación de Variables ---

# Aplicar primera diferencia (diff) a las variables I(1) 
d.ta_preferencial <- diff(datos$ta_preferencial)
# d.ta_preferencial <- diff(d.ta_preferencial)

d.dprime <- diff(datos$dprime)
d.dprime <- diff(d.dprime)

d.ipc <- diff(datos$ipc)
# d.ipc <- diff(d.ipc)

d.exp_tc_media_12 <- diff(datos$exp_tc_media_12)
# d.exp_tc_media_12 <- diff(d.exp_tc_media_12)

d.ta_90d <- diff(datos$ta_90d)
# d.ta_90d <- diff(d.ta_90d)

# La variable I(0) 'sorpresatc_exp' se mantiene en niveles.
# Se debe alinear su longitud, removiendo la primera observación.
sorpresa_tc_aligned <- datos$sorpresa_tc_exp[-1]
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
  d.dprime = d.dprime,
  d.ipc = d.ipc,
  d.ta_90d = d.ta_90d,
  d.ta_preferencial = d.ta_preferencial,
  sorpresatc_exp = sorpresa_tc_aligned,
  d.exp_tc_media_12 = d.exp_tc_media_12
)


var_data_ts <- ts(var_data_stationary, start = c(2009, 8), frequency = 12)


lag_selection <- VARselect(var_data_ts, lag.max = 12, type = "const")


p_optimo <- lag_selection$selection[[1]]


# --- 6. Estimación y Validación del Modelo VAR ---

# Estimación del modelo VAR(p) [5, 13]
var_model <- VAR(var_data_ts, p = p_optimo)
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

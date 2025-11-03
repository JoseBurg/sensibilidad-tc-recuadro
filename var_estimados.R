library(bsvars)
library(bsvarSIGNs)
library(dplyr)
library(ggplot2)
library(vars)
library(urca)
library(tseries)

datos <- readxl::read_excel("./data-exp-tc.xlsx") |> 
  mutate(
    variacion_tc_esperada = (exp_tc_media_mes / lag(tc_mes)) - 1
  ) |> 
  dplyr::select(-c(fecha, exp_tc_media_mes, tc_mes)) |> 
  na.omit()


data_ts <- ts(datos, start = c(2009, 7), frequency = 12)

# Se espera p > 0.05 (No Estacionaria)
adf.test(datos$ta_90d)
adf.test(datos$exp_tc_media_12)
adf.test(datos$ta_preferencial)
adf.test(datos$dprime)
adf.test(datos$ipc)

# Se espera p < 0.05 (Estacionaria) 
adf.test(datos$variacion_tc_esperada) # Estacionaria

# --- Transformación de Variables --------
# Aplicar primera diferencia (diff) a las variables 
d.ta_preferencial <- diff(datos$ta_preferencial)

# Segunda diferencia
d.dprime <- diff(datos$dprime)
d.dprime <- diff(d.dprime)

d.ipc <- diff(datos$ipc)

d.exp_tc_media_12 <- diff(datos$exp_tc_media_12)

d.ta_90d <- diff(datos$ta_90d)

# Se debe alinear su longitud, removiendo la primera observación.
variacion_tc_esperada <- datos$variacion_tc_esperada[-1]
N <- length(d.ipc)



adf.test(d.exp_tc_media_12)
adf.test(d.ta_90d)
adf.test(d.ta_preferencial)
adf.test(d.dprime)
adf.test(d.ipc)

# Construir el sistema estacionario final usando cbind() 
# Orden: Exógenas -> Endógenas

var_data_stationary <- cbind(
  d.dprime = d.dprime,
  d.ipc = d.ipc,
  d.ta_90d = d.ta_90d,
  d.ta_preferencial = d.ta_preferencial,
  variacion_tc_esperada = variacion_tc_esperada,
  d.exp_tc_media_12 = d.exp_tc_media_12
)


var_data_ts <- ts(var_data_stationary, start = c(2009, 7), frequency = 12)

lag_selection <- VARselect(var_data_ts, lag.max = 12, type = "const")

p_optimo <- lag_selection$selection[[1]]




var_model <- VAR(var_data_ts, p = p_optimo)

serial_test <- serial.test(var_model, lags.pt = p_optimo + 5, type = "PT.asymptotic")
print(serial_test)

stability_test <- stability(var_model, type = "OLS-CUSUM")
plot(stability_test)


arch_test <- arch.test(var_model, lags.multi = p_optimo)
print(arch_test)


irf_results <- irf(var_model, 
                   impulse = "variacion_tc_esperada", 
                   response = "d.exp_tc_media_12", 
                   n.ahead = 12,         # Horizonte de 24 meses
                   ortho = FALSE,         # Usar Cholesky (dependiente del orden)
                   boot = TRUE,          # Calcular bandas de confianza [32]
                   ci = 0.95,            # Nivel de confianza del 95%
                   runs = 1000)          # Iteraciones de Bootstrap

# Visualización de la FIR 
# Este es el gráfico principal que responde la pregunta
plot(irf_results, 
     ylab = "Respuesta de d.exp_tc_media_12", 
     xlab = "Meses después del choque",
     main = "Respuesta de Expectativas LP ante Choque de Variación esperada")



irf_all_responses <- irf(var_model, impulse = "variacion_tc_esperada", n.ahead = 9, ortho = TRUE, boot = TRUE)
plot(irf_all_responses)


irf_all_impulses <- irf(var_model, response = "d.exp_tc_media_12", n.ahead = 9, ortho = TRUE, boot = TRUE)
plot(irf_all_impulses)





# Restricciones de signos -------------------------------------------------
# ---------------------------
# 1. Cargar librerías
# ---------------------------
library(vars)
library(svars)

# ---------------------------
# 2. Construir el VAR
# ---------------------------
var_model <- VAR(var_data_stationary, lag.max = 12, ic = "AIC")

summary(var_model)

# ---------------------------
# 3. Identificación estructural (Sign restrictions + structural break)
# ---------------------------

# Supón que el quiebre estructural ocurrió alrededor del año 2020 (por ejemplo, pandemia)
# Tu serie inicia en 2009M7 → el mes 2020M3 corresponde aprox. a observación 129
# (ajusta según la longitud real de tu serie)
x1 <- id.cv(var_model, SB = 129)

summary(x1)

# ---------------------------
# 4. Reordenar las columnas (según patrón de signo esperado)
# ---------------------------
# Queremos que el primer choque sea el de variación del tipo de cambio esperado (positivo)
# y que tenga un efecto positivo sobre las expectativas de largo plazo.
x1$B <- x1$B[, c(5, 6, 1, 2, 3, 4)]  # mover columna 5 (variacion_tc_esperada) al frente
x1$B[, 1] <- abs(x1$B[, 1])          # asegurar que sea un choque positivo

# ---------------------------
# 5. Impulse Response Function (IRF)
# ---------------------------
irf_result <- irf(x1, n.ahead = 12)
plot(irf_result, scales = 'free_y', main = "Respuestas a Choque de Variación Esperada del TC")

# ---------------------------
# 6. Aplicar restricción de impacto (opcional)
# ---------------------------
# Ejemplo: suponemos que d.dprime no reacciona instantáneamente al choque cambiario
restMat <- matrix(NA, ncol = 6, nrow = 6)
restMat[1, 1] <- 0  # Restringe impacto inmediato sobre la primera variable

x2 <- id.cv(var_model, SB = 129, restriction_matrix = restMat)
summary(x2)




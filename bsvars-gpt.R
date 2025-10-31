library(bsvars)
library(bsvarSIGNs)
library(dplyr)
library(ggplot2)
library(vars)
library(urca)
library(tseries)

datos <- readxl::read_excel("./data.xlsx")
# Nota: Las expectativas estan adelantada un mes, las de octubre son las de septiembre,
# esto se debe a que la encuesta inicia a principio de mes, y el observado es del 
# último día laborable

base_datos <- datos |> 
  mutate(
    remesas = log(remesas),
    across(c(
      ipc_vi, remesas  
    ),
    \(x) {x - lag(x)}
    )) |> 
  na.omit() |> 
  dplyr::select(exp_tc_mes, everything())


var_data <- base_datos %>% 
  dplyr::select(-c(periodo, diff_mn_me, remesas, tp_interbancarios))


var_ts <- ts(var_data, start = c(2010, 3), frequency = 12)

Y <- var_ts

# Determinar el número de variables (N) y el orden del rezago (p)
N <- ncol(Y)
p <- 5 # Número de rezagos (ajusta según el criterio de selección)

# Especificar la estructura base del BSVAR (incluye una constante por defecto)
specifications <- specify_bsvar$new(Y, p = p)

# 3. Definición de Restricciones de Signo
# ---
# Ejemplo de Restricción Hipotética:
# Variables: [1]=tc_venta, [2]=exp_tc_mes, [3]=..., [N]=...
# Shock de Interés: Shock 1 (Columna 1)
# Restricción: Shock 1 causa (+) en Var 1 y (-) en Var 2 en el impacto (h=0)

# Crear la matriz de restricciones N x N, inicializada a 0
sign_restrictions <- matrix(0, nrow = N, ncol = N)

# Aplicar las restricciones al Shock 1 (Columna 1)
sign_restrictions[1, 1] <- 1  # Variable 1 (+)
sign_restrictions[2, 1] <- -1 # Variable 2 (-)

# NOTA: Deja las demás celdas en 0 si no tienes restricciones para otros shocks.

# Especificar la estructura de restricciones de signo
sign_specification <- specify_signs_bsvar(
  bsvar_obj = specifications,
  sign_array = sign_restrictions,
  horizon = 1 # Restricción aplicada solo en el periodo de impacto (h=0)
)

# 4. Estimación MCMC e Identificación
# ---
# Definir parámetros de la simulación
N_draws <- 10000 # Número total de draws de MCMC
M_burnin <- 5000  # Número de draws a descartar (burn-in)

# Estimar el modelo, identificando los shocks con las restricciones de signo
cat("Iniciando la estimación MCMC con restricciones de signo...\n")
model_sign_restricted <- estimate_bsvar(
  sign_specification,
  N = N_draws,
  M = M_burnin
)
cat("Estimación finalizada.\n")


# 5. 📉 Análisis de Resultados (Funciones de Impulso-Respuesta - FIR)
# ---
# Calcular las IRFs
horizon_irf <- 24 # Horizonte de 24 periodos (e.g., meses)
cat("Calculando las funciones de impulso-respuesta (FIR) a 24 periodos...\n")
irf_sign <- compute_impulse_responses(
  model_sign_restricted,
  horizon = horizon_irf
)

# Graficar la FIR (muestra la mediana y las bandas de credibilidad)
cat("Generando gráficos de FIR.\n")
plot(irf_sign)
model_bsvar <- var_ts |> 
  specify_bsvar$new(p = 6) |> 
  estimate(S = 1000)


model_bsvar |> 
  compute_impulse_responses(12) |> 
  plot()


model_bsvar |> 
  forecast(horizon = 8) |> 
  plot()

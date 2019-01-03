## calculate R2cv of vector obs and pred, not for matrix
## R2 = max(0, 1-MSE/Var(x))
## usually, input is log(x)

R2cv <- function(obs, pred, log = T){
# first check the length of obs and pred
  if(!length(obs) == length(pred))
    print("error: the lengths of obs and pred are different.") else
    {
      # There probably NA exists, so take sum(!is.na(obs)) instead of length(obs)
      # MSE = sum((exp(obs)-exp(pred))^2, na.rm = T)/sum(!is.na(obs))
      if(log){
      MSE = mean((exp(obs)-exp(pred))^2, na.rm = T)
      var = mean((exp(obs)-mean(exp(obs), na.rm = T))^2, na.rm = T)
            R2cv = max(0, 1-MSE/var)
      } else
      {
        MSE = mean((obs-pred)^2, na.rm = T)
        var = mean((obs-mean(obs, na.rm = T))^2, na.rm = T)
        R2cv = max(0, 1-MSE/var)
      }
  }
  return(list(r2=R2cv, rmse=sqrt(MSE)))
}

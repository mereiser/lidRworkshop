mystdMetrics = function(z, i, rn, minht, above)
{
  "mode" <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  zmean = mean(z)
  hmin = if (length(z[z > minht])>0) min(z[z > minht]) else 1.37
  
  metrics = list(
    hmin = if (length(z[z > minht])>0) min(z[z > minht]) else 1.37,
    # hmax = max(z),
    hmean = mean(z[z > minht]),
    hmode = mode(z[z > minht]),
    hmedian = median(z[z > minht]),
    hsd = sd(z[z > minht]),
    hvar = var(z[z > minht]),
    hcv = (sd(z[z > minht])/mean(z[z > minht]))*100,
    hkurtosis = (sum((z[z > minht] - mean(z[z > minht]))^3)/length(z[z > minht]))/(sum((z[z > minht] - mean(z[z > minht]))^2)/length(z[z > minht]))^(3/2),
    hskewness = length(z[z > minht]) * sum((z[z > minht] - mean(z[z > minht]))^4)/(sum((z[z > minht] - mean(z[z > minht]))^2)^2),
    hmad.mean = mad(z[z > minht], center =  mean(z[z > minht])),
    hmad.med = mad(z[z > minht], center =  median(z[z > minht])),
    hmad.mode = mad(z[z > minht], center =  mode(z[z > minht])),
    hquad.mean = (mean(z[z > minht]^2))^(1/2),
    hcubic.mean = (mean(z[z > minht]^3))^(1/3),
    
    Canopy.relief.ratio = (mean(z[z > minht])-hmin)/(max(z)-hmin),
    
    Pentage.first.returns.Above.XX = length(rn[z > above & rn == 1])/length(rn[rn == 1])*100,
    # Percentage.all.returns.above.XX = length(rn[z > above])/length(rn)*100,
    All.returns.above.XX.Total.first.returns.100 = length(rn[z > above])/length(rn[rn == 1])*100,
    First.returns.above.XX = length(rn[z > above & rn == 1]),
    All.returns.above.XX = length(rn[z > above]),
    
    Percentage.first.returns.above.mean = length(rn[z > mean(z) & rn == 1])/length(rn[rn == 1])*100,
    Percentage.first.returns.above.mode = length(rn[z > mode(z) & rn == 1])/length(rn[rn == 1])*100,
    # Percentage.all.returns.above.mean = length(rn[z > mean(z)])/length(rn)*100,
    Percentage.all.returns.above.mode =  length(rn[z > mode(z)])/length(rn)*100,
    
    All.returns.above.mean.Total.first.returns.100 = length(rn[z > mean(z)])/length(rn[rn == 1])*100,
    All.returns.above.mode.Total.first.returns.100 = length(rn[z > mode(z)])/length(rn[rn == 1])*100,
    First.returns.above.mean = length(rn[z > mean(z) & rn == 1]),
    First.returns.above.mode = length(rn[z > mode(z) & rn == 1]),
    All.returns.above.mean = length(rn[z > mean(z)]),
    All.returns.above.mode = length(rn[z > mode(z)])
  )
  
  return(metrics)
}
MagnetarLightCurve <- function(t,Ep,tp){
     # Generate a light curve from the model of magnetar spin-down
     # Lsd(t) = (Ep/tp) * (1 + t/tp)^(-2)
     # Does not pass the diffusion process, nor leakage
     # t = time (after explosion in rest frame) corresponding to the calculated spin-down luminosity (day)
     # Ep = initial spin-down energy (erg)
     # tp = initial spin-down timescale (day)
     # output = total spin-down luminosity at a given time (erg/s)
     return( (Ep/(tp * 24 * 60 * 60.0)) * (1.0 + t/tp)^(-2) )
}

DiffusionMagnetarLightCurve <- function(Ep,tp,tmax,Esn,tLC,tDiff=tLC,r0=0.0,v=1e4){
     # Generate a light curve from the diffusion approximation with magnetar spin-down engine from 1 - tmax days
     # This calls function MagnetarLightCurve
     # By default, the function assumes small initial radius
     # Time step for integration is 1 day
     # Edge effect: set zero at 0 and tmax+1 days
     # Ep = initial spin-down energy (erg)
     # tp = spin-down timescale (day)
     # tmax = maximum time (after explosion in rest frame) corresponding to the calculated output luminosity (day)
     # Esn = Explosion energy (erg)
     # tLC = effective light curve timescale (i.e., = sqrt(2 * tDiff * tExp)) (day)
     # tDiff = diffusion timescale = tLC (default) (day)
     # r0 = 0 (default for assuming small initial radius) (cm)
     # v = ejecta velocity (km/s)
     # return(time grid in days, spin-down input luminosity in erg/s, spin-down diffused output luminosity in erg/s, fireball luminosity in erg/s)
     
     tgrid = c(0:tmax) # generate time grids (day)
     
     # values will be often used in the calculation
     a = (tgrid/tLC)^2 + (2.0 * r0*1e-5 * tgrid)/(v*60.0*60*24 * tLC^2) # changing units
     b = (r0*1e-5)/(v*60.0*60*24 * tLC) + tgrid/tLC
     c = exp(-a)
     d = exp(+a)
     ####
     
     LoutFireball = (Esn/(tDiff * 24.0 * 60 * 60)) * c # fireball term
     
     # magnetar term
     Linp = MagnetarLightCurve(tgrid,Ep,tp)
     integrand = d * b * Linp
     integral = c()
     integral[1] = 0.0
     LoutMagnetar = c()
     LoutMagnetar[1] = 0.0
     for (i in 2:tmax){
          integral[i] = 0.5 * (integrand[i-1] + integrand[i]) * 1.0
          LoutMagnetar[i] = LoutMagnetar[i-1] + integral[i]
     }
     LoutMagnetar[tmax+1] = 0.0
     LoutMagnetar = (2.0/tLC) * c * LoutMagnetar
     ###
     
     return(list(tgrid = tgrid[2:tmax], Linp = Linp[2:tmax], LoutMagnetar = LoutMagnetar[2:tmax], LoutFireball = LoutFireball[2:tmax]))
}

# example:
tmax = 300.0
lout = DiffusionMagnetarLightCurve(Ep=1e52,tp=50.0,tmax=tmax,Esn=1e51,tLC=80.0)
str(lout)
plot(lout$tgrid,lout$LoutMagnetar,type="l")
lines(lout$tgrid,lout$LoutFireball)
lines(lout$tgrid,lout$Linp)

     

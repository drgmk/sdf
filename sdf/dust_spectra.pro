;; make "real grain" dust spectra

;; this needs to be independent of the stellar properties, so the dust location cannot be
;; a parameter. thus use bb temperature to set radius, and then proceed. for small dust
;; the temperature is a weak function of stellar temperature, so assume a G-type star

function rgspec,bbtemp,dmin,q,wavs

  nwavs = N_ELEMENTS(wavs)
  
  ;; assume a composition (astrosil)
  ptrdielec = getdielecw('lg_amsil.dat')

  ;; assume a stellar spectrum
  starparams,'G0V',teff=tstar,lstar=lstar,mstar=mstar
  ptrsspec = kurucz(t=tstar)

  dmax = 1e5                    ; in um
  ndiams = 500
  logdmin     = alog10(dmin)
  logdmax     = alog10(dmax)
  logdiams    = logdmin+(logdmax-logdmin)*findgen(ndiams)/(ndiams-1)
  diams       = 10.0^logdiams

  ;; get the radius for this bb temp
  rs = (278.3/bbtemp)^2 * sqrt(lstar)
  nrs = 1

  ;; dust temperatures
  temps = dusttempw(diams,rs,ptrdielec=ptrdielec,lstar=lstar,tstar=tstar,ptrsspec=ptrsspec,nstep=200,/CALLC,/NOWIN)

  ;; Q_abs as function of wavelength for each size
  qabs        = fltarr(ndiams,nwavs)
  for i=0,ndiams-1 do qabs[i,*] = getqw(wavs,diams[i],ptrdielec=ptrdielec,qout=1,/CALLC)

  ;; Average QabsBnu over the size distributions
  pr = fltarr(nwavs)            ; Jy/sr
  for i=0,nwavs-1 do begin
     pr[i] = int_tabulated(logdiams,qabs[*,i]*bnuw(wavs[i],temps[*,0])* $
                           sigmadbarw(logdiams=logdiams,q=q))
  endfor

  return,pr
end

;; wavelength range slightly greater than used in sdf, readsav seems
;; to restore things with slightly different ranges
nw = 500
wavs = logarrw(min=0.05,max=1e4,n=nw)
wavs = [0.05,wavs,1e4+1]
nw = N_ELEMENTS(wavs)

nt = 50
temps = logarrw(min=10,max=1e3,n=nt)
nd = 30
dmins = logarrw(min=0.1,max=100,n=nd)
nq = 15
qs = linarrw(min=10.1/6.,max=11.9/6.,n=nq)

pr = fltarr(nw,nt,nd,nq)

for i=0,nt-1 do begin
   for j=0,nd-1 do begin
      for k=0,nq-1 do begin
         pr[*,i,j,k] = rgspec(temps[i],dmins[j],qs[k],wavs)
;         plot,wavs,pr[*,i,j,k],/xl,/yl
      endfor
   endfor
   counter,i+1,nt,'done:'
endfor

save,filename='rg_pr.xdr',pr,wavs,temps,dmins,qs

end

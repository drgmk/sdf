;; make "real grain" dust spectra

;; this needs to be independent of the stellar properties, so the dust location cannot be
;; a parameter. thus use bb temperature to set radius, and then proceed. for small dust
;; the temperature is a weak function of stellar temperature, so assume a G-type star

;; wavelength range slightly greater than used in sdf, readsav seems
;; to restore things with slightly different ranges
nwavs = 500
wavs = logarrw(min=0.05,max=1e4,n=nwavs)
wavs = [0.05,wavs,1e4+1]
nwavs = N_ELEMENTS(wavs)

;; assume a composition (astrosil)
ptrdielec = getdielecw('lg_amsil.dat')

;; assume a stellar spectrum
starparams,'G0V',teff=tstar,lstar=lstar,mstar=mstar
ptrsspec = kurucz(t=tstar)

dmax = 1e5                    ; in um
ndiams = 500

nt = 50
tbbs = logarrw(min=10,max=1e3,n=nt)
nd = 30
dmins = logarrw(min=0.1,max=100,n=nd)
nq = 15
qs = linarrw(min=10.1/6.,max=11.9/6.,n=nq)

pr = fltarr(nwavs,nt,nd,nq)

for j=0,nd-1 do begin

    logdmin     = alog10(dmins[j])
    logdmax     = alog10(dmax)
    logdiams    = logdmin+(logdmax-logdmin)*findgen(ndiams)/(ndiams-1)
    diams       = 10.0^logdiams

    ;; Q_abs as function of wavelength for each size
    qabs        = fltarr(ndiams,nwavs)
    for l=0,ndiams-1 do qabs[l,*] = getqw(wavs,diams[l],ptrdielec=ptrdielec,qout=1,/CALLC)

    for i=0,nt-1 do begin

        ;; get the radius for this bb temp
        rs = (278.3/tbbs[i])^2 * sqrt(lstar)

        ;; dust temperatures
        temps = dusttempw(diams,rs,ptrdielec=ptrdielec,lstar=lstar,tstar=tstar, $
                          ptrsspec=ptrsspec,nstep=200,/CALLC,/NOWIN)

        for k=0,nq-1 do begin

            ;; Average QabsBnu over the size distributions
            prtmp = fltarr(nwavs)            ; Jy/sr
            sdb = sigmadbarw(logdiams=logdiams,q=qs[k])
            for l=0,nwavs-1 do begin
                prtmp[l] = int_tabulated(logdiams,qabs[*,l]*bnuw(wavs[l],temps[*,0])*sdb)
            endfor

            pr[*,i,j,k] = prtmp
;            plot,wavs,pr[*,i,j,k],/xl,/yl
        endfor
    print,'iteration '+STRN(j)+':'+STRN(i)
    endfor
;    counter,j+1,nd,'done:'
endfor

save,filename='rg_pr.xdr',pr,wavs,tbbs,dmins,qs

end

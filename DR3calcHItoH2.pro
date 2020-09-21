PRO timemark

  common time,runtime,resultfolder
  print,'Start!',systime()
  runtime=string(systime())
  resultfolder=strcompress('/mnt/science1/bliu/HINSA/code/20190813/')
  FILE_MKDIR,resultfolder
  
  return
END

FUNCTION tau0,p,THINSA,Tb,TH
  Tc=3.8           ;galactic value =3.5, LG2003, LMC =3.8 calced from Hughes 2007.
  tauHI=0.7        ;average by cosmponent, D94+M97+MZ2000
  tauf=tauHI*(1-p)
  return,alog((p*Tb+(Tc-TH)*(1-tauf))/(p*Tb+(Tc-TH)*(1-tauf)-THINSA))
END

FUNCTION Kfactor,THINSA,Tb,TH
  return,tau0(0.5,THINSA,Tb,TH)/tau0(1.0,THINSA,Tb,TH)
END

PRO DR3calcHItoH2fill

timemark
common time

;HIfile='/mnt/science1/bliu/HINSA/data/HI.fits'
;HINSAfile='/mnt/science1/bliu/HINSA/code/20190813/DR3HINSAcube.fits'
;COfile='/mnt/science1/bliu/HINSA/data/LMC_MAGMA_DR3.co.base.fits'
;HINSAtaufile='/mnt/science1/bliu/HINSA/code/20190813/DR3tau0cube.fits'
;COMasksumfile='/mnt/science1/bliu/HINSA/data/gm/LMC_MAGMA_DR3.co.gm.mask_sum.fits'
;
;HI=readfits(HIfile,HIhead)
;HINSA=readfits(HINSAfile,HINSAhead)
;CO=readfits(COfile,COhead)
;HINSAtau=readfits(HINSAtaufile,HINSAtauhead)
;COmasksum=readfits(COMasksumfile,COmasksumhead)
;
;Xfactor=4E20
;Tfactor=1.      ;20190813 changed to 1, previously used 2.0
;deltz=fxpar(COhead,'CDELT3')/1.0e3
;z12CO=((dindgen(fxpar(COhead,'NAXIS3'))+0.5-fxpar(COhead,'CRPIX3'))*fxpar(COhead,'CDELT3')+fxpar(COhead,'CRVAL3'))/1.0e3
;zHI=((dindgen(fxpar(HINSAhead,'NAXIS3'))+0.5-fxpar(HINSAhead,'CRPIX3'))*fxpar(HINSAhead,'CDELT3')+fxpar(HINSAhead,'CRVAL3'))/1.0e3
;
;NHIcube=HINSAtau*0
;NH2cube=HINSAtau*0
;ratio=HINSAtau*0          ;before Kfactor correction
;ratioprime=HINSAtau*0     ;after Kfactor correction
;
;xmax=(size(HINSA))(1)-1
;
;for x=0,(size(HINSA))(1)-1 do begin
;  print,x,xmax
;  for y=0,(size(HINSA))(2)-1 do begin
;    HIspec=reform(HI(x,y,*))
;    HINSAspec=reform(HINSA(x,y,*))
;    totalHINSA=total(HINSAspec)
;    if totalHINSA gt 0 then begin
;      xyad,HINSAhead,x,y,ra,dec
;      adxy,COhead,ra,dec,xCO,yCO
;      
;      ;if COmasksum(xCO,yCO) gt 15 then begin
;      
;      COspec=reform(CO(xCO,yCO,*))
;      sum12CO=total(COspec)
;      NH2=sum12CO*deltz*Xfactor
;      
;      gsresult=gaussfit(z12CO,COspec,COgsfit,Nterms=3)
;      gsresult=gaussfit(zHI,HINSAspec,HINSAgsfit,Nterms=3)
;      TB12CO=COgsfit(0)
;      tauzero=HINSAtau(x,y)
;      Tex12CO=5.5/alog(1+5.5/TB12CO+0.82)
;      NHI=1.95E18*tauzero*(2.3548*HINSAgsfit[2])*Tex12CO
;      TH=Tfactor*Tex12CO
;      THINSA=max(HINSAspec)
;      Tb=HIspec((where(HINSAspec eq THINSA))(0))+THINSA
;      
;      HItoH2=NHI/NH2
;      
;      NHIcube(x,y)=NHI
;      NH2cube(x,y)=NH2
;      ratio(x,y)=HItoH2
;      HItoH2prime=HItoH2*Kfactor(THINSA,Tb,TH)
;      ratioprime(x,y)=HItoH2prime
;      
;      ;endif      
;    endif
;  endfor
;endfor
;  
;  save,NHIcube,filename='/mnt/science1/bliu/HINSA/code/20190813/NHIcube.sav'
;  save,NH2cube,filename='/mnt/science1/bliu/HINSA/code/20190813/NH2cube.sav'
;  save,ratio,filename='/mnt/science1/bliu/HINSA/code/20190813/ratio.sav'
;  save,ratioprime,filename='/mnt/science1/bliu/HINSA/code/20190813/ratioprime.sav';Tfactor=1

  restore,'/mnt/science1/bliu/HINSA/code/20190813/NHIcube.sav' 
  restore,'/mnt/science1/bliu/HINSA/code/20190813/NH2cube.sav' 
  restore,'/mnt/science1/bliu/HINSA/code/20190813/ratio.sav' 
  restore,'/mnt/science1/bliu/HINSA/code/20190813/ratioprime.sav'

  set_plot,'PS'
  filename='/mnt/science1/bliu/HINSA/code/20190813/'+'DR3calcHItoH2.eps'
  Device,/ENCAPSUL, BITS_PER_PIXEL=8, /COLOR, FILENAME=filename, XSIZE=26, YSIZE=20
  !P.Color = '000000'xL
  !P.Background = 'FFFFFF'xL
  loadct,39,ncolors=255,bottom=1,/silent

  ratioLG2003=[0.0012,0.00012,0.0016,0.0016,0.0023,0.00026,0.0029,0.0003,0.0036,0.00036,0.0039,0.0004,0.0041,0.0044,0.00048,0.00048,0.00055,0.00057,0.0073,0.00074,0.00081,0.00009]
  ;ratioKG2010=[0.00301,0.00279,0.00381,0.0028,0.0004,0.00014,0.00123,0.00077,0.00075,0.00229,0.00169,0.00029,0.00178,0.00204,0.00189,0.00339,0.00071,0.00079,0.00096,0.00013,0.00063,0.00191,0.00248,0.00212,0.00101,0.00139,0.00124,0.00178,0.00035,0.003,0.00166,0.00227,0.00136,0.00094,0.00066,0.00094]
  ;ratioLGKG=[0.0012,0.00012,0.0015,0.0016,0.0016,0.0023,0.00026,0.0029,0.0003,0.0036,0.00036,0.0039,0.0004,0.0041,0.0044,0.00048,0.00048,0.00055,0.00057,0.0073,0.00074,0.00081,0.00009,0.00301,0.00279,0.00381,0.0028,0.0004,0.00014,0.00123,0.00077,0.00075,0.00229,0.00169,0.00029,0.00178,0.00204,0.00189,0.00339,0.00071,0.00079,0.00096,0.00013,0.00063,0.00191,0.00248,0.00212,0.00101,0.00139,0.00124,0.00178,0.00035,0.003,0.00166,0.00227,0.00136,0.00094,0.00066,0.00094]
  ratioKG2010LOS=[0.00301,0.00123,0.00169,0.00079,0.00063,0.00191,0.00139,0.00166,0.00066]
  ratioKG2010Components=[0.00279,0.00381,0.0028,0.0004,0.00014,0.00077,0.00075,0.00229,0.00029,0.00178,0.00204,0.00189,0.00339,0.00071,0.00096,0.00013,0.00063,0.00248,0.00212,0.00101,0.00124,0.00178,0.00035,0.003,0.00227,0.00136,0.00094,0.00094]
  ratioZuo2018=[0.021]
  ratioall=[ratioLG2003,ratioKG2010LOS,ratioZuo2018]
  ratioLGKG=[ratioLG2003,ratioKG2010Components]

  mi=-6
  ma=-1
  bin=0.1
  yaxisscalefactor=5
    
;  cgplot,indgen((ma-mi)/bin)*bin+mi+0.5*bin,histogram(alog10(ratio),binsize=bin,max=ma,min=mi),psym=10,xstyle=1,yrange=[0,180],xtitle=textoidl('HINSA HI to H2 ratio[LOG(N_{HINSA}/N_{H_{2}})]'),ytitle='Counts',linestyle=1,thick=3
;  oplot,indgen((ma-mi)/bin)*bin+mi+0.5*bin,histogram(alog10(ratioprime),binsize=bin,max=ma,min=mi),psym=10,thick=3
;  bin=0.5
;  oplot,indgen((ma-mi)/bin)*bin+mi+0.5*bin,histogram(alog10(ratioLG2003),binsize=bin,max=ma,min=mi)*yaxisscalefactor,psym=10,color=250,linestyle=0,thick=3
;  bin=0.5
;  oplot,indgen((ma-mi)/bin)*bin+mi+0.5*bin,histogram(alog10(ratioKG2010Components),binsize=bin,max=ma,min=mi)*yaxisscalefactor,psym=10,color=150,linestyle=0,thick=3
;  bin=0.5
;  oplot,indgen((ma-mi)/bin)*bin+mi+0.5*bin,histogram(alog10(ratioLGKG),binsize=bin,max=ma,min=mi)*yaxisscalefactor,psym=10,color=90,linestyle=0,thick=3
;  
  cgHistoplot, alog10(ratioprime), BINSIZE=bin, maxinput=ma, mininput=mi, MAX_VALUE=180, /fill, xtitle=textoidl('HINSA HI to H2 ratio[LOG(N_{HINSA}/N_{H_{2}})]'),ytitle='Counts'
  cgHistoplot, alog10(ratio), BINSIZE=bin, maxinput=ma, mininput=mi, MAX_VALUE=180, /oplot
  bin=0.5
  cgHistoplot, alog10(ratioLG2003)*yaxisscalefactor,binsize=bin, maxinput=ma, mininput=mi, MAX_VALUE=180, /oplot, /fill, polycolor='white'
;  bin=0.5
;  oplot,indgen((ma-mi)/bin)*bin+mi+0.5*bin,histogram(alog10(ratioKG2010Components),binsize=bin,max=ma,min=mi)*yaxisscalefactor,psym=10,color=150,linestyle=0,thick=3
;  bin=0.5
;  oplot,indgen((ma-mi)/bin)*bin+mi+0.5*bin,histogram(alog10(ratioLGKG),binsize=bin,max=ma,min=mi)*yaxisscalefactor,psym=10,color=90,linestyle=0,thick=3
;  
  
  
;  print,'LMC',mean(ratio(where((Alog10(ratio) gt -6)*(Alog10(ratio) lt -1))))
;  print,'LMC corrected',mean(ratioprime(where((Alog10(ratioprime) gt -6)*(Alog10(ratioprime) lt -1))))
;  print,'LG2003',mean(ratioLG2003)
;  print,'KG2010',mean(ratioKG2010Components)
;  print,'LG+KG',mean(ratioLGKG)
;  print,'LG+KG+Zuo',mean(ratioall)
;  print,'std LGKG',STDDEV(ratioLGKG)
;  print,'std LMC corrected',STDDEV(ratioprime(where(ratioprime gt 0)))
  
;  bin=0.1
;  x=indgen((ma-mi)/bin+1)*bin+mi+0.5*bin
;  y=histogram(alog10(ratio),binsize=bin,max=ma,min=mi)
;;  print,size(x)
;;  print,size(y)
;  result=gaussfit(x,y,A,nterms=3)
;  print,A
;  
;  bin=0.1
;  x=indgen((ma-mi)/bin+1)*bin+mi+0.5*bin
;  y=histogram(alog10(ratioprime),binsize=bin,max=ma,min=mi)
;;  print,size(x)
;;  print,size(y)
;  result=gaussfit(x,y,A,nterms=3)
;  print,A
;  
;  bin=0.5
;  x=indgen((ma-mi)/bin+1)*bin+mi+0.5*bin
;  y=histogram(alog10(ratioLG2003),binsize=bin,max=ma,min=mi)
;  result=gaussfit(x,y,A,nterms=3)
;  print,A
;  
;  bin=0.5
;  x=indgen((ma-mi)/bin+1)*bin+mi+0.5*bin
;  y=histogram(alog10(ratioKG2010Components),binsize=bin,max=ma,min=mi)
;  result=gaussfit(x,y,A,nterms=3)
;  print,A
;  
;  bin=0.5
;  x=indgen((ma-mi)/bin+1)*bin+mi+0.5*bin
;  y=histogram(alog10(ratioLGKG),binsize=bin,max=ma,min=mi)
;  result=gaussfit(x,y,A,nterms=3)
;  print,A
  
  
  
  ;LG2003 reports 1.5*10^-3, KG2010 reports 1.6*10^-3
  xyouts,-5.8,170,'Black: LMC (Dashed: before optical depth correction)',charsize=1.5
  xyouts,-5.8,159,'Red: LG2003',charsize=1.5
  xyouts,-5.8,148,'Green: KG2010',charsize=1.5
  xyouts,-5.8,137,'Blue: LG2003 + KG2010',charsize=1.5
  
  Device, /close
  
;  print,size(ratioprime)
;  print,size(where(ratioprime gt 0))
;;  print,size(ratioLGKG)
;;  print,'LMC corrected:'
;  print,size(ratioprime(where((Alog10(ratioprime) gt -6)*(Alog10(ratioprime) lt -1))))
;;   print,'LMC before corrected:'
;;   print,ratio(where((Alog10(ratio) gt -6)*(Alog10(ratio) lt -1)))

print,'finished'

END

FUNCTION tau0,p,Tab
Tb=80.
Tc=3.5
TH=10.
tauHI=0.7  ;average by component, D94+M97+MZ2000
tauf=tauHI*(1-p)
return,alog((p*Tb+(Tc-TH)*(1-tauf))/(p*Tb+(Tc-TH)*(1-tauf)-Tab))
END

FUNCTION Kfactor,p,THINSA
return,tau0(p,THINSA)/tau0(1,THINSA)
END

PRO calccorrection

print,'start!'

  set_plot,'PS'
  filename='/nfshome/byliu/Data/20150911result/'+'calccorrection.eps'
  Device,/ENCAPSUL, BITS_PER_PIXEL=8, /COLOR, FILENAME=filename, XSIZE=24, YSIZE=20
  !P.Color = '000000'xL
  !P.Background = 'FFFFFF'xL
  loadct,39,ncolors=255,bottom=1,/silent

x=indgen(10000.)/10000.
cgplot,x,Kfactor(x,1),xrange=[0,1],yrange=[1,100],/YLOG,linestyle=0,xtitle='p',ytitle='K'
cgoplot,x,Kfactor(x,2),linestyle=1
cgoplot,x,Kfactor(x,3),linestyle=1
cgoplot,x,Kfactor(x,4),linestyle=1
cgoplot,x,Kfactor(x,5),linestyle=1
cgoplot,x,Kfactor(x,6),linestyle=1
cgoplot,x,Kfactor(x,7),linestyle=1
cgoplot,x,Kfactor(x,8),linestyle=1
cgoplot,x,Kfactor(x,9),linestyle=1
;cgoplot,x,Kfactor(x,0.5),linestyle=1
cgoplot,x,Kfactor(x,0.01),linestyle=0
cgoplot,x,Kfactor(x,10),linestyle=0
cgtext,0.05,4,'THINSA=0.1 K',charsize=1.5
cgtext,0.05,70,'THINSA=1.0 K',charsize=1.5
cgtext,0.48,2.4,'THINSA=10 K',charsize=1.5
cgtext,0.75,70,'Tb=80 K',charsize=1.5
cgtext,0.75,58,'Tc=3.5 K',charsize=1.5
cgtext,0.75,48,'TH=10 K',charsize=1.5
cgtext,0.75,40,'TauHI=0.7',charsize=1.5

print,Kfactor(0.5,0.1),Kfactor(0.5,1),Kfactor(0.5,10)

Device, /close

print,'Done!'

END

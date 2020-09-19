PRO timemark

  common time,runtime,resultfolder
  print,'Start!',systime()
  runtime=string(systime())
  resultfolder=strcompress('/nfshome/byliu/Data/result/'+string(runtime)+'/')
  FILE_MKDIR,resultfolder
  
  return
END

PRO readHIfits

  common HIfits,a,head,bw,freq,nx,ny,nz,crvalx,cdeltax,crpixx,crvaly,cdeltay,crpixy,crvalz,cdeltaz,crpixz,ctypey,x,y,z,zhel,vchannel
  
  fitsname='/nfshome/byliu/Data/HI.fits'
  a=mrdfits(fitsname) ;HI
  
  head=headfits(fitsname);read the header of the fits file to a vector
  
  bw= fxpar(head,'BW'); band width (MHz)
  freq= fxpar(head,'LINEFREQ'); central frequency (GHz)
  
  nx= fxpar(head,'NAXIS1'); number of elements in the first dimension
  ny= fxpar(head,'NAXIS2');
  nz= fxpar(head,'NAXIS3');
  crvalx= fxpar(head,'CRVAL1'); reference value of the first dimension
  cdeltax = fxpar(head,'CDELT1'); increasement of the first dimension
  ; in units of degree, when calculate physical
  ; scale, must changed to arcdegree
  crpixx= fxpar(head,'CRPIX1'); reference position of the first dimension
  crvaly= fxpar(head,'CRVAL2');
  cdeltay = fxpar(head,'CDELT2');
  crpixy= fxpar(head,'CRPIX2');
  crvalz= fxpar(head,'CRVAL3');
  cdeltaz = fxpar(head,'CDELT3');
  crpixz= fxpar(head,'CRPIX3');
  ctypey= fxpar(head,'CTYPE2');
  
  x=(dindgen(nx)+0.5-crpixx)*cdeltax+crvalx;
  y=(dindgen(ny)+0.5-crpixy)*cdeltay+crvaly;
  z=(dindgen(nz)+0.5-crpixz)*cdeltaz+crvalz;
  z=z/1.0e3
  zhel=z
  vchannel=abs(cdeltaz)/1.0e3; km/s?
  
  print,'HI fits read! ',systime()
  
END

PRO readCOfits

  common COfits,a12co,head12co,bw12co,freq12co,nx12co,ny12co,nz12co,crvalx12co,cdeltax12co,crpixx12co,$
    crvaly12co,cdeltay12co,crpixy12co,crvalz12co,cdeltaz12co,crpixz12co,x12co,y12co,z12co,vchannel12co
    
  fitsname12co='/nfshome/byliu/Data/MAGMADR3/Basic/LMC_MAGMA_DR3.co.base.fits'
  a12co=mrdfits(fitsname12co)
  
  head12co=headfits(fitsname12co);read the header of the fits file to a vector
  
  bw12co= fxpar(head12co,'BW'); band width (MHz)
  freq12co= fxpar(head12co,'LINEFREQ'); central frequency (GHz)
  
  nx12co= fxpar(head12co,'NAXIS1'); number of elements in the first dimension
  ny12co= fxpar(head12co,'NAXIS2');
  nz12co= fxpar(head12co,'NAXIS3');
  crvalx12co= fxpar(head12co,'CRVAL1'); reference value of the first dimension
  cdeltax12co = fxpar(head12co,'CDELT1'); increasement of the first dimension
  ; in units of degree, when calculate physical
  ; scale, must changed to arcdegree
  crpixx12co= fxpar(head12co,'CRPIX1'); reference position of the first dimension
  crvaly12co= fxpar(head12co,'CRVAL2');
  cdeltay12co = fxpar(head12co,'CDELT2');
  crpixy12co= fxpar(head12co,'CRPIX2');
  crvalz12co= fxpar(head12co,'CRVAL3');
  cdeltaz12co = fxpar(head12co,'CDELT3');
  crpixz12co= fxpar(head12co,'CRPIX3');
  
  x12co=(dindgen(nx12co)+0.5-crpixx12co)*cdeltax12co+crvalx12co;
  y12co=(dindgen(ny12co)+0.5-crpixy12co)*cdeltay12co+crvaly12co;
  z12co=(dindgen(nz12co)+0.5-crpixz12co)*cdeltaz12co+crvalz12co;
  z12co=z12co/1.0e3
  vchannel12co=abs(cdeltaz12co)/1.0e3; km/s?
  
  print,'CO fits read! ',systime()
  
END

PRO readCOmasksumfits

  common COmasksumfits,a12comasksum
    
  fitsname12comasksum='/nfshome/byliu/Data/MAGMADR3/gm/LMC_MAGMA_DR3.co.gm.mask_sum.fits'
  a12comasksum=mrdfits(fitsname12comasksum)
  
  print,'COmask fits read! ',systime()
  
END

PRO readxyHItoCOfits

  common xyHItoCOfits,xyHItoCO
    
  fitsnamexyHItoCO='/nfshome/byliu/Data/20150401result/xyHItoCO.fit'
  xyHItoCO=mrdfits(fitsnamexyHItoCO)
  
  print,'xyHItoCO fits read!',systime()
  
END

FUNCTION chi2,x
  COMMON share,specmodel,zchi,Tc,THchi,vH,sigmaHchi,vchannelchi
  tau0 = x[0]                                       ;tau0 is the optical depth of HINSA at vH?
  tau=tau0*exp(-(zchi-vH)^2/(2.0*sigmaHchi^2))
  specorigin=specmodel*exp(tau)-Tc-THchi*(exp(tau)-1.0)
  dspec=deriv(zchi,specorigin) ; derivative
  ddspec=deriv(zchi,dspec) ; second order derivative
  result=total(ddspec^2*vchannelchi)
  return,result
END

PRO calcHINSA,HIarray,COarray,z,z12CO,vchannel,HINSAarray,TH,tau0,sigmaH
  COMMON share
  zchi=z
  vchannelchi=vchannel
  spec=HIarray
  specmodel=spec
  
  gsfit=gaussfit(z12CO,COarray,result,NTerms=3)
  model12CO=(result(0)*exp(-((z12CO-result(1))/result(2))^2/2.0))

;      limit=0
;      while result(0) lt 0 do begin
;        ;print,'Negtive peak detected. Fit again.'
;        spec12CO=spec12CO-model12CO
;        result=Gsfitpart(z12CO,spec12CO,3,xmings,xmaxgs)
;        model12CO=(result(0)*exp(-((z12CO-result(1))/result(2))^2/2.0))
;        limit++
;        IF (limit eq 5) THEN BREAK
;      endwhile
;      if result(0) lt 0 then begin
;        HINSAarray=make_array((size(reform(HIarray[0,0,*])))[1],value=0.)
;      print,'Still neg peak.'
;      endif else begin
        ck=1.3806488e-16    ;cm^2 g s^-2 K^-1, Boltzmann constant
        Tc=0                ;??? ;3.5K by Lockman 1984, as cited in Li 2003
        TB0=result(0)       ;K, 12CO observed brightness temperature at line center
        T0=5.53             ;K, 12CO J=1-0 equivalent temperature
        Tbg=2.73            ;K, CMB
        Tex=T0/alog(1+1/(TB0/T0+1/(exp(T0/Tbg)-1))) ;K, 12CO excit T, formula(A1), Krco 2008   ;13.7.23
        Texcorrected=Tex*1.2
        tk=Texcorrected              ;K, kinetic temperature
        TH=tk               ;???
        vH=result(1)        ;CO & HINSA central velocity
        AvoCons=6.02214129d23 ;Avogadro constant
        mH=1.00794/AvoCons  ;g
        mCO=28.010/AvoCons  ;g
        sigmaH=sqrt((result(2))^2+2*ck*TH*(1/mH-1/mCO)/1.0d10) ;line width, using 201406 version draft's formula
;        sigma_th=sqrt(2kt/m) by http://esclab.tw/wiki/index.php/Thermal_Line_Broadening
;        sigmacaution=0
;        sigmacriteria=1.5
;        if sigmaH gt sigmacriteria then begin
;          sigmaHcopy=sigmaH
;          sigmaH=sigmacriteria
;          sigmacaution=1
;        endif
        
        sigmaHchi=sigmaH
        THchi=TH
        
        ftol = 1.0e-5
        point=AMOEBA(ftol,SCALE=[1.0],P0=[0.15],FUNCTION_VALUE=fval,FUNCTION_NAME='chi2')
        tau0=point(0)
        tau=tau0*exp(-(z-vH)^2/(2.0*sigmaH^2))
        specorigin=spec*exp(tau)-Tc-TH*(exp(tau)-1.0)
        ; =>HINSA=Sori-Sobs=(Sobs-Tex)*(exp(tau)-1)-Tc, thus HINSA profile is not central symmetric.
        HINSAarray=specorigin-spec
        
END

PRO DR3calcHINSA

timemark
readHIfits
readCOfits
readCOmasksumfits
readxyHItoCOfits

common time
common HIfits
common COfits
common COmasksumfits
common xyHItoCOfits

xyad,head12CO,0,0,a0,d0
xyad,head12CO,nx12CO-1,ny12CO-1,a1,d1
adxy,head,a0,d0,x0,y0
adxy,head,a1,d1,x1,y1
print,x0,y0,x1,y1
HINSAcube=make_array(nx,ny,nz)
THcube=make_array(nx,ny)
tau0cube=make_array(nx,ny)
sigmaHcube=make_array(nx,ny)

for xHI=fix(x0),fix(x1) do begin
  print,xHI
  for yHI=fix(y0),fix(y1) do begin
;xyad,head12CO,1334,ny12CO-141,a0,d0
;adxy,head,a0,d0,x0,y0
;
;for x=x0,x0 do begin
;  for y=y0,y0 do begin
;    xyad,head,x,y,RA,DEC
;    adxy,head12CO,RA,DEC,x12CO,y12CO
    x12CO=xyHItoCO(xHI,yHI,2)
    y12CO=xyHItoCO(xHI,yHI,3)
    if (x12CO gt 0) AND (y12CO gt 0) AND (x12CO lt nx12CO) AND (y12CO lt ny12CO) then begin
      if a12comasksum(x12CO,y12CO) ge 3 then begin
;      print,xHI,yHI
      HIarray=reform(a(xHI,yHI,*))
      COarray=reform(a12CO(x12CO,y12CO,*))
      calcHINSA,HIarray,COarray,z,z12CO,vchannel,HINSAarray,TH,tau0,sigmaH
;      plot,z,HIarray
;      oplot,z,HINSAarray
;      oplot,z12CO,COarray*5
      ;print,TH1,tau01,sigmaH1
      HINSAcube(xHI,yHI,*)=HINSAarray
      THcube(xHI,yHI)=TH
      tau0cube(xHI,yHI)=tau0
      sigmaHcube(xHI,yHI)=sigmaH
      endif
    endif
  endfor
endfor

print,'lets write'

outhead0=head
outhead1=head
sxaddpar,outhead1,'NAXIS',2
sxdelpar,outhead1,['NAXIS3','CRPIX3','CDELT3','CRVAL3','CTYPE3']

outfilename=strcompress(resultfolder+'DR3HINSAcube.fits')
mwrfits,HINSAcube,outfilename,outhead0
outfilename=strcompress(resultfolder+'DR3THcube.fits')
mwrfits,THcube,outfilename,outhead1
outfilename=strcompress(resultfolder+'DR3tau0cube.fits')
mwrfits,tau0cube,outfilename,outhead1
outfilename=strcompress(resultfolder+'DR3sigmaHcube.fits')
mwrfits,sigmaHcube,outfilename,outhead1

print,'finished ',systime()
END










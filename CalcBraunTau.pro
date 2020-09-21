PRO timemark

  common time,runtime,resultfolder
  print,'Start!',systime()
  runtime=string(systime())
  resultfolder=strcompress('/nfshome/byliu/Data/result/'+string(runtime)+'/')
  FILE_MKDIR,resultfolder

  return
END

PRO readRBHI0fits

  common RBHI0fits,aRBHI0,headRBHI0,bwRBHI0,freqRBHI0,nxRBHI0,nyRBHI0,nzRBHI0,crvalxRBHI0,cdeltaxRBHI0,crpixxRBHI0,$
    crvalyRBHI0,cdeltayRBHI0,crpixyRBHI0,crvalzRBHI0,cdeltazRBHI0,crpixzRBHI0,xRBHI0,yRBHI0,zRBHI0,vchannelRBHI0

  fitsnameRBHI0='/nfshome/byliu/Data/Braun2012/lmcnhv.fits'
  aRBHI0=mrdfits(fitsnameRBHI0)

  headRBHI0=headfits(fitsnameRBHI0);read the header of the fits file to a vector

  bwRBHI0= fxpar(headRBHI0,'BW'); band width (MHz)
  freqRBHI0= fxpar(headRBHI0,'LINEFREQ'); central frequency (GHz)

  nxRBHI0= fxpar(headRBHI0,'NAXIS1'); number of elements in the first dimension
  nyRBHI0= fxpar(headRBHI0,'NAXIS2');
  crvalxRBHI0= fxpar(headRBHI0,'CRVAL1'); reference value of the first dimension
  cdeltaxRBHI0 = fxpar(headRBHI0,'CDELT1'); increasement of the first dimension
  ; in units of degree, when calculate physical
  ; scale, must changed to arcdegree
  crpixxRBHI0= fxpar(headRBHI0,'CRPIX1'); reference position of the first dimension
  crvalyRBHI0= fxpar(headRBHI0,'CRVAL2');
  cdeltayRBHI0 = fxpar(headRBHI0,'CDELT2');
  crpixyRBHI0= fxpar(headRBHI0,'CRPIX2');

  xRBHI0=(dindgen(nxRBHI0)+0.5-crpixxRBHI0)*cdeltaxRBHI0+crvalxRBHI0;
  yRBHI0=(dindgen(nyRBHI0)+0.5-crpixyRBHI0)*cdeltayRBHI0+crvalyRBHI0;

  print,'RBHI0 fits read! ',systime()

END

PRO readRBHI1fits

  common RBHI1fits,aRBHI1,headRBHI1,bwRBHI1,freqRBHI1,nxRBHI1,nyRBHI1,nzRBHI1,crvalxRBHI1,cdeltaxRBHI1,crpixxRBHI1,$
    crvalyRBHI1,cdeltayRBHI1,crpixyRBHI1,crvalzRBHI1,cdeltazRBHI1,crpixzRBHI1,xRBHI1,yRBHI1,zRBHI1,vchannelRBHI1

  fitsnameRBHI1='/nfshome/byliu/Data/Braun2012/lmcnhcg_C60.fits'
  aRBHI1=mrdfits(fitsnameRBHI1)

  headRBHI1=headfits(fitsnameRBHI1);read the header of the fits file to a vector

  bwRBHI1= fxpar(headRBHI1,'BW'); band width (MHz)
  freqRBHI1= fxpar(headRBHI1,'LINEFREQ'); central frequency (GHz)

  nxRBHI1= fxpar(headRBHI1,'NAXIS1'); number of elements in the first dimension
  nyRBHI1= fxpar(headRBHI1,'NAXIS2');
  crvalxRBHI1= fxpar(headRBHI1,'CRVAL1'); reference value of the first dimension
  cdeltaxRBHI1 = fxpar(headRBHI1,'CDELT1'); increasement of the first dimension
  ; in units of degree, when calculate physical scale, must changed to arcdegree
  crpixxRBHI1= fxpar(headRBHI1,'CRPIX1'); reference position of the first dimension
  crvalyRBHI1= fxpar(headRBHI1,'CRVAL2');
  cdeltayRBHI1 = fxpar(headRBHI1,'CDELT2');
  crpixyRBHI1= fxpar(headRBHI1,'CRPIX2');

  xRBHI1=(dindgen(nxRBHI1)+0.5-crpixxRBHI1)*cdeltaxRBHI1+crvalxRBHI1;
  yRBHI1=(dindgen(nyRBHI1)+0.5-crpixyRBHI1)*cdeltayRBHI1+crvalyRBHI1;

  print,'RBHI1 fits read! ',systime()

END

PRO readLBYTaufits

  common LBYTaufits,aLBYTau,headLBYTau,bwLBYTau,freqLBYTau,nxLBYTau,nyLBYTau,nzLBYTau,crvalxLBYTau,cdeltaxLBYTau,crpixxLBYTau,$
    crvalyLBYTau,cdeltayLBYTau,crpixyLBYTau,crvalzLBYTau,cdeltazLBYTau,crpixzLBYTau,xLBYTau,yLBYTau,zLBYTau,vchannelLBYTau

  fitsnameLBYTau='/nfshome/byliu/Data/20150401result/3/DR3tau0cube.fits'
  aLBYTau=mrdfits(fitsnameLBYTau)

  headLBYTau=headfits(fitsnameLBYTau);read the header of the fits file to a vector

  bwLBYTau= fxpar(headLBYTau,'BW'); band width (MHz)
  freqLBYTau= fxpar(headLBYTau,'LINEFREQ'); central frequency (GHz)

  nxLBYTau= fxpar(headLBYTau,'NAXIS1'); number of elements in the first dimension
  nyLBYTau= fxpar(headLBYTau,'NAXIS2');
  crvalxLBYTau= fxpar(headLBYTau,'CRVAL1'); reference value of the first dimension
  cdeltaxLBYTau = fxpar(headLBYTau,'CDELT1'); increasement of the first dimension
  ; in units of degree, when calculate physical
  ; scale, must changed to arcdegree
  crpixxLBYTau= fxpar(headLBYTau,'CRPIX1'); reference position of the first dimension
  crvalyLBYTau= fxpar(headLBYTau,'CRVAL2');
  cdeltayLBYTau = fxpar(headLBYTau,'CDELT2');
  crpixyLBYTau= fxpar(headLBYTau,'CRPIX2');

  xLBYTau=(dindgen(nxLBYTau)+0.5-crpixxLBYTau)*cdeltaxLBYTau+crvalxLBYTau;
  yLBYTau=(dindgen(nyLBYTau)+0.5-crpixyLBYTau)*cdeltayLBYTau+crvalyLBYTau;

  print,'LBYTau fits read! ',systime()

END

PRO CalcBraunTau

timemark
readRBHI0fits
readRBHI1fits
readLBYTaufits

common time
common RBHI0fits
common RBHI1fits
common LBYTaufits

print,size(aRBHI0)
print,size(aRBHI1)
print,size(aLBYTau)
print,size(aRBHI1/aRBHI0)
print,size(aRBHI0(*))

  set_plot,'PS'
  filename='/nfshome/byliu/Data/20150911result/'+'CalcBraunTau.eps'
  Device,/ENCAPSUL, BITS_PER_PIXEL=8, /COLOR, FILENAME=filename, XSIZE=20, YSIZE=19
  !P.Color = '000000'xL
  !P.Background = 'FFFFFF'xL
  loadct,39,ncolors=255,bottom=1,/silent

cgscatter2d,aLBYTau(where((aLBYTau gt 0)*(aLBYTau lt 2)*((aRBHI1/aRBHI0) gt 1)*((aRBHI1/aRBHI0) lt 7.389))),alog((aRBHI1/aRBHI0)(where((aLBYTau gt 0)*(aLBYTau lt 2)*((aRBHI1/aRBHI0) gt 1)*((aRBHI1/aRBHI0) lt 7.389)))),fit=0,/NODISPLAY,xrange=[0,2],yrange=[0,2],xstyle=1,ystyle=1,psym=cgSYMCAT(3),BACKGROUND='white',AXISCOLOR='black',color='black',title='',xtitle='Tau',ytitle='Tau(RB)'

Device, /close

print,'finished'

END

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
  z=z/1.0e3                                ;Helio-centric
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

PRO makelog,xleft,xright,ytop,ydown

  common time
  RunLog=strcompress('/nfshome/byliu/Data/result/Log.txt')
  openw,lun,RunLog,/get_lun,/APPEND
  printf, lun, runtime
  printf, lun, 'calcxyHItoCO based on HI & 12CO map'
  printf, lun, strcompress('x12='+string(xleft)+','+string(xright)+'; y12='+string(ytop)+','+string(ydown))
  Free_lun, lun
  
END

PRO calcxyHItoCO

  common time
  common HIfits
  common COfits
  
  timemark
  readHIfits
  readCOfits
  
;  xyad,head12CO,0,0,a0,d0
;  xyad,head12CO,nx12CO-1,ny12CO-1,a1,d1
;  adxy,head,a0,d0,x0,y0
;  adxy,head,a1,d1,x1,y1
  
  print,'define studied region'
  xleft =0
  ytop  =0
  xright=nx-1
  ydown =ny-1
  makelog,xleft,xright,ytop,ydown
  
  print,'initialize the .fit file'
  xyHItoCOfits=make_array(nx,ny,4,value=0.)
  print,'size of fits',size(adxyadfits)
  outfilename=strcompress(resultfolder+'xyHItoCO.fit')
  outhead = head
  mwrfits,xyHItoCOfits,outfilename,outhead
  
  print,'modify fits header'
  sxaddpar,outhead,'NAXIS3',4
  sxaddpar,outhead,'CRPIX3',0
  sxaddpar,outhead,'CDELT3',1
  sxaddpar,outhead,'CRVAL3',0
  sxaddpar,outhead,'CTYPE3','xyHItoCO'
  modfits,outfilename,0,outhead
  
  print,'start FOR loop'
  for x=xleft,xright do begin
    for y=ny-1-ydown,ny-1-ytop do begin
      
      xyad,head,x,y,RA,DEC
      adxy,head12CO,RA,DEC,x12CO,y12CO
      
      xyHItoCOfits(x,y,*)=[RA,DEC,x12CO,y12CO]
      
    endfor
    print,x-xleft+1,' of',xright-xleft+1,'  ',systime()
    if ((x-xleft+1) mod 100) eq 0 then begin
      modfits,outfilename,xyHItoCOfits,outhead
    endif
  endfor
  
  print,'Finished! ',systime()

END
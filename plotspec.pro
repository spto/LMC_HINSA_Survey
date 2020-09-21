PRO timemark

  common time,runtime,resultfolder
  print,'Start!',systime()
  runtime=string(systime())
  resultfolder=strcompress('/mnt/science1/bliu/HINSA/code/20190813/')
  FILE_MKDIR,resultfolder

  return
END

PRO plotspec

  common time
  timemark

  HIfile='/mnt/science1/bliu/HINSA/data/HI.fits'
  HINSAfile='/mnt/science1/bliu/HINSA/code/20190813/DR3HINSAcube.fits'
  COfile='/mnt/science1/bliu/HINSA/data/LMC_MAGMA_DR3.co.base.fits'

  HI=readfits(HIfile,HIhead)
  HINSA=readfits(HINSAfile,HINSAhead)
  CO=readfits(COfile,COhead)

  deltz=fxpar(COhead,'CDELT3')/1.0e3
  z12CO=((dindgen(fxpar(COhead,'NAXIS3'))+0.5-fxpar(COhead,'CRPIX3'))*fxpar(COhead,'CDELT3')+fxpar(COhead,'CRVAL3'))/1.0e3
  z=((dindgen(fxpar(HINSAhead,'NAXIS3'))+0.5-fxpar(HINSAhead,'CRPIX3'))*fxpar(HINSAhead,'CDELT3')+fxpar(HINSAhead,'CRVAL3'))/1.0e3

  ralist=[71.841239,71.895764,72.257477,72.294002,72.373006,72.598984,72.959219,73.069755,73.212667,73.523764,73.891089,74.073412,74.676177,75.948556,76.358901,77.483187,78.337628,78.356254,78.46388,78.638794,80.554041,81.091014,81.214431,81.294496,81.47361,83.853115,83.949401,83.971081,84.623886,84.898579,84.935074,85.01205,85.318982,85.847604,86.177372,86.572708,86.757814]
  declist=[-67.195076,-67.204442,-68.604774,-68.589958,-68.504026,-69.504405,-69.354991,-66.894605,-68.064296,-69.192526,-66.471359,-66.624016,-66.133117,-67.309752,-66.898346,-68.892573,-69.38429,-67.471547,-67.12856,-68.769225,-67.961909,-68.428101,-69.672431,-69.677802,-66.23535,-67.581782,-69.218882,-69.039694,-69.036005,-69.771228,-69.625439,-69.857336,-70.925518,-69.420056,-69.470454,-69.64071,-70.769733]

  set_plot,'PS'

  for i=0,36 do begin
    ra=ralist[i]
    dec=declist[i]
    adxy,HINSAhead,ra,dec,xH,yH
    adxy,COhead,ra,dec,xCO,yCO
    xH=round(xH)
    yH=round(yH)

    spec=HI(xH,yH,*)
    spec12CO=CO(xCO,yCO,*)
    specHINSA=HINSA(xH,yH,*)
    specorigin=spec+specHINSA

    filename='/mnt/science1/bliu/HINSA/code/20190813/'+string(i+1,format='(I0)')+'.eps'
    Device,/ENCAPSUL, BITS_PER_PIXEL=8, /COLOR, FILENAME=filename, XSIZE=24, YSIZE=16
    !P.Color = '000000'xL
    !P.Background = 'FFFFFF'xL
    loadct,39,ncolors=255,bottom=1,/silent
    print,filename

    cgplot,z,spec,xrange=[200,350],yrange=[-10,100],linestyle=0,$
      thick=2,color=255,xtitle='!17v(km/s)',ytitle='!17T(K)',psym=10
    cgoplot,z12co,spec12co*10,linestyle=0,thick=1,color=240,psym=10
    cgoplot,z(where(specHINSA ge 1e-5)),specorigin(where(specHINSA ge 1e-5)),linestyle=1,thick=2,color='black',psym=10
    cgoplot,z,specHINSA,thick=2,color=80,psym=10

    figtitle=strcompress('No.'+string(i+1,format='(I0)'))
    xyouts,207,90,figtitle,charsize=2
;    xyouts,xmin+3,0.95*max(spec),'Black solid: observed HI',charsize=0.8
;    xyouts,xmin+3,0.90*max(spec),'Red solid: observed CO',charsize=0.8,color=240
;    xyouts,xmin+3,0.85*max(spec),'Yellow dash dot: gauss-fitted CO',charsize=0.8,color=190
;    xyouts,xmin+3,0.80*max(spec),'Orange dash dot: original HI',charsize=0.8,color=220
;    xyouts,xmin+3,0.75*max(spec),'Blue dash dot: HINSA',charsize=0.8,color=80

  endfor

  Device, /close

print,'finished'
END

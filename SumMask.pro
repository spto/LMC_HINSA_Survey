PRO timemark

  common time,runtime,resultfolder
  print,'Start!',systime()
  runtime=string(systime())
  resultfolder=strcompress('/nfshome/byliu/Data/result/'+string(runtime)+'/')
  FILE_MKDIR,resultfolder
  
  return
END

PRO SumMask

timemark
common time
print,resultfolder

print,'Read-in the mask'
fitsname='/nfshome/byliu/Data/MAGMADR3/gm/LMC_MAGMA_DR3.co.gm.mask.fits'
maskfit=mrdfits(fitsname)
maskhead=headfits(fitsname)
print,'Mask fits read!'

print,size(maskfit)

sumed=make_array(1600,1533)
for x=0,1599 do begin
  print,x
  for y=0,1532 do begin
      sumed(x,y)=total(maskfit(x,y,*),/NAN)
      ;if sumed(x,y) gt 0 then print,x,y,sumed(x,y)
  endfor
endfor

outfilename=strcompress(resultfolder+'LMC_MAGMA_DR3.co.gm.mask_sum.fits')
outhead=maskhead
sxaddpar,outhead,'NAXIS',2
sxdelpar,outhead,['NAXIS3','CRPIX3','CDELT3','CRVAL3','CTYPE3']
mwrfits,sumed,outfilename,outhead

print,'finished!'
END
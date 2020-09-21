PRO timemark

  common time,runtime,resultfolder
  print,'Start!',systime()
  runtime=string(systime())
  resultfolder=strcompress('/mnt/science1/bliu/HINSA/code/20190813/')
  FILE_MKDIR,resultfolder

  return
END

PRO DR3calcrings

restore,'/mnt/science1/bliu/HINSA/code/20190813/ratioprime.sav',/verbose

HINSAtaufile='/mnt/science1/bliu/HINSA/code/20190813/DR3tau0cube.fits'
HINSAtau=readfits(HINSAtaufile,HINSAtauhead)

xsize=(size(ratioprime))(1)
ysize=(size(ratioprime))(2)

print,xsize,ysize

ratioinrings=make_array(7,xsize,ysize)

ra0=79.4
dec0=-69.03
for x=0,xsize-1 do begin
  for y=0,ysize-1 do begin
    if ratioprime(x,y) ne 0 then begin
 
      xyad,HINSAtauhead,x,y,ra,dec
      ;print,ra,dec

      result=map_2points(ra0,dec0,ra,dec)
      dinkpc=result(0)/57.3*50.
      PA=result(1)+90
      ;print,dinkpc,PA

      PAring=[-3.1,-5.7,1.4,-7.1,-19.4,-21.4,-10]
      rinkpc=PAring*0
      for i=0,6 do begin
        dPA=PA-PAring[i]
        r=(i+1)*0.5
        rinkpc[i]=r*sqrt(cos(dPA)^2+0.93^2*sin(dPA)^2)
      endfor
      ;print,rinkpc
      
      ringno=(where(dinkpc lt rinkpc))(0)
      ;print,ringno,ratioprime(x,y)
      
      ratioinrings(ringno-1,x,y)=ratioprime(x,y)
      
    endif
  endfor
endfor

mi=-6
ma=-1
bin=0.1


ratio0=(ratioinrings(0,*,*))(where(ratioinrings(0,*,*) ne 0))
ratio1=(ratioinrings(1,*,*))(where(ratioinrings(1,*,*) ne 0))
ratio2=(ratioinrings(2,*,*))(where(ratioinrings(2,*,*) ne 0))
ratio3=(ratioinrings(3,*,*))(where(ratioinrings(3,*,*) ne 0))
ratio4=(ratioinrings(4,*,*))(where(ratioinrings(4,*,*) ne 0))
ratio5=(ratioinrings(5,*,*))(where(ratioinrings(5,*,*) ne 0))
ratio6=(ratioinrings(6,*,*))(where(ratioinrings(6,*,*) ne 0))
print,size(ratio0)
print,size(ratio1)
print,size(ratio2)
print,size(ratio3)
print,size(ratio4)
print,size(ratio5)
print,size(ratio6)

set_plot,'PS'
filename='/mnt/science1/bliu/HINSA/code/20190813/'+'DR3calcHItoH2_rings.eps'
Device,/ENCAPSUL, BITS_PER_PIXEL=8, /COLOR, FILENAME=filename, XSIZE=22, YSIZE=18
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct,39,ncolors=255,bottom=1,/silent

;ratioLG2003=[0.0012,0.00012,0.0015,0.0016,0.0016,0.0023,0.00026,0.0029,0.0003,0.0036,0.00036,0.0039,0.0004,0.0041,0.0044,0.00048,0.00048,0.00055,0.00057,0.0073,0.00074,0.00081,0.00009]
;ratioKG2010=[0.00301,0.00279,0.00381,0.0028,0.0004,0.00014,0.00123,0.00077,0.00075,0.00229,0.00169,0.00029,0.00178,0.00204,0.00189,0.00339,0.00071,0.00079,0.00096,0.00013,0.00063,0.00191,0.00248,0.00212,0.00101,0.00139,0.00124,0.00178,0.00035,0.003,0.00166,0.00227,0.00136,0.00094,0.00066,0.00094]
;ratioLGKG=[0.0012,0.00012,0.0015,0.0016,0.0016,0.0023,0.00026,0.0029,0.0003,0.0036,0.00036,0.0039,0.0004,0.0041,0.0044,0.00048,0.00048,0.00055,0.00057,0.0073,0.00074,0.00081,0.00009,0.00301,0.00279,0.00381,0.0028,0.0004,0.00014,0.00123,0.00077,0.00075,0.00229,0.00169,0.00029,0.00178,0.00204,0.00189,0.00339,0.00071,0.00079,0.00096,0.00013,0.00063,0.00191,0.00248,0.00212,0.00101,0.00139,0.00124,0.00178,0.00035,0.003,0.00166,0.00227,0.00136,0.00094,0.00066,0.00094]

mi=-6
ma=-1
bin=0.1

cgplot,indgen((ma-mi)/bin+1)*bin+mi,histogram(alog10(ratio0),binsize=bin,max=ma,min=mi),psym=10,thick=3,xstyle=1,yrange=[0,60],color=30,xtitle='HINSA HI to H2 ratio in rings',ytitle='counts'
oplot,indgen((ma-mi)/bin+1)*bin+mi,histogram(alog10(ratio1),binsize=bin,max=ma,min=mi),psym=10,thick=3,color=70
oplot,indgen((ma-mi)/bin+1)*bin+mi,histogram(alog10(ratio2),binsize=bin,max=ma,min=mi),psym=10,thick=3,color=110
oplot,indgen((ma-mi)/bin+1)*bin+mi,histogram(alog10(ratio3),binsize=bin,max=ma,min=mi),psym=10,thick=3,color=150
oplot,indgen((ma-mi)/bin+1)*bin+mi,histogram(alog10(ratio4),binsize=bin,max=ma,min=mi),psym=10,thick=3,color=170
oplot,indgen((ma-mi)/bin+1)*bin+mi,histogram(alog10(ratio5),binsize=bin,max=ma,min=mi),psym=10,thick=3,color=220
oplot,indgen((ma-mi)/bin+1)*bin+mi,histogram(alog10(ratio6),binsize=bin,max=ma,min=mi),psym=10,thick=3,color=254

print,size(indgen((ma-mi)/bin+1)*bin+mi)
print,size(histogram(alog10(ratio0),binsize=bin,max=ma,min=mi))

result=gaussfit(indgen((ma-mi)/bin+1)*bin+mi,histogram(alog10(ratio0),binsize=bin,max=ma,min=mi),C0,Nterms=3)
result=gaussfit(indgen((ma-mi)/bin+1)*bin+mi,histogram(alog10(ratio1),binsize=bin,max=ma,min=mi),C1,Nterms=3)
result=gaussfit(indgen((ma-mi)/bin+1)*bin+mi,histogram(alog10(ratio2),binsize=bin,max=ma,min=mi),C2,Nterms=3)
result=gaussfit(indgen((ma-mi)/bin+1)*bin+mi,histogram(alog10(ratio3),binsize=bin,max=ma,min=mi),C3,Nterms=3)
result=gaussfit(indgen((ma-mi)/bin+1)*bin+mi,histogram(alog10(ratio4),binsize=bin,max=ma,min=mi),C4,Nterms=3)
result=gaussfit(indgen((ma-mi)/bin+1)*bin+mi,histogram(alog10(ratio5),binsize=bin,max=ma,min=mi),C5,Nterms=3)
result=gaussfit(indgen((ma-mi)/bin+1)*bin+mi,histogram(alog10(ratio6),binsize=bin,max=ma,min=mi),C6,Nterms=3)

print,C0
print,C1
print,C2
print,C3
print,C4
print,C5
print,C6

xyouts,-5.8,57,'0.0 < a < 0.5 kpc',charsize=1.3,charthick=4,color=30
xyouts,-5.8,55,'0.5 < a < 1.0 kpc',charsize=1.3,charthick=4,color=70
xyouts,-5.8,53,'1.0 < a < 1.5 kpc',charsize=1.3,charthick=4,color=110
xyouts,-5.8,51,'1.5 < a < 2.0 kpc',charsize=1.3,charthick=4,color=150
xyouts,-5.8,49,'2.0 < a < 2.5 kpc',charsize=1.3,charthick=4,color=170
xyouts,-5.8,47,'2.5 < a < 3.0 kpc',charsize=1.3,charthick=4,color=220
xyouts,-5.8,45,'3.0 < a < 3.5 kpc',charsize=1.3,charthick=4,color=254

Device, /close

print,'finished'
END
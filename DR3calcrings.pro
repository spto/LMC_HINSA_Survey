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

print,'finished'
END

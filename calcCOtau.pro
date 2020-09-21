PRO deriveCOtau,TB12CO,TB13CO,tau12CO,iter,Tex12CO

;Zeng, Mao, Pei 2006
;assume 12CO optically thick

Tex12CO1=5.5/alog(1+5.5/(TB12CO+0.82))
tau13CO=-1*alog((1-TB13CO/5.3*((exp(5.3/Tex12CO1-1))^(-1)-0.16)^(-1)))
ratio12to13=50
tau12CO=ratio12to13*tau13CO


Tex12CO=0
for i=1,iter do begin
Tex12CO=5.5/alog(1+1/(TB12CO/(5.5*(1-exp(-1*tau12CO)))+0.15))  ;smaller tau12CO, larger Tex
tau13CO=-1*alog((1-TB13CO/5.3*((exp(5.3/Tex12CO-1))^(-1)-0.16)^(-1)))
tau12CO=ratio12to13*tau13CO
endfor

print,Tex12CO1,Tex12CO;,TB12CO,TB13CO,tau13CO,tau12CO

END

PRO calcCOtau2019

;;readin the files:DR3 12CO, 13CO
;;to select:
;;/mnt/science1/bliu/HINSA/data/gm/LMC_MAGMA_DR3.co.gm.mask_sum.fits : >5 (there is 12CO)
;;/mnt/science1/bliu/HINSA/data/13CO_snrpk/*** : >4 (13CO strong enough)
;;to calculate:
;;/mnt/science1/bliu/HINSA/data/LMC_MAGMA_DR3.co.base.fits
;;/mnt/science1/bliu/HINSA/data/13CO_maps/***
;;
;;scan the 13CO files
;;if CO.gm.mask_sum > 5 && snrpk > 4, fit the 13CO, fit the 12CO, calc tau
;;
;;standard CO analysis?
;;save the result on a 12CO-sized array
;;make histogram of CO tau
;;
;COMasksumfile='/mnt/science1/bliu/HINSA/data/gm/LMC_MAGMA_DR3.co.gm.mask_sum.fits'
;list13snrpk=file_search('/mnt/science1/bliu/HINSA/data/13CO_snrpk/','*')
;COspecfile='/mnt/science1/bliu/HINSA/data/LMC_MAGMA_DR3.co.base.fits'
;list13map=file_search('/mnt/science1/bliu/HINSA/data/13CO_maps/','*')
;;
;COmasksum=readfits(COMasksumfile,COmasksumhead)
;COspec=readfits(COspecfile,COspechead)
;snrpk13CO1=readfits(list13snrpk(0),snrpk13COhead1)
;snrpk13CO2=readfits(list13snrpk(1),snrpk13COhead2)
;snrpk13CO3=readfits(list13snrpk(2),snrpk13COhead3)
;snrpk13CO4=readfits(list13snrpk(3),snrpk13COhead4)
;snrpk13CO5=readfits(list13snrpk(4),snrpk13COhead5)
;snrpk13CO6=readfits(list13snrpk(5),snrpk13COhead6)
;map13CO1=readfits(list13map(0),map13COhead1)
;map13CO2=readfits(list13map(1),map13COhead2)
;map13CO3=readfits(list13map(2),map13COhead3)
;map13CO4=readfits(list13map(3),map13COhead4)
;map13CO5=readfits(list13map(4),map13COhead5)
;map13CO6=readfits(list13map(5),map13COhead6)
;;
;print,size(snrpk13CO1)
;print,size(map13CO1)
;;
;map13COhead=[ptr_new(map13COhead1,/No_copy),ptr_new(map13COhead2,/No_copy),$
;  ptr_new(map13COhead3,/No_copy),ptr_new(map13COhead4,/No_copy),$
;  ptr_new(map13COhead5,/No_copy),ptr_new(map13COhead6,/No_copy)]
;map13CO=[ptr_new(map13CO1,/No_copy),ptr_new(map13CO2,/No_copy),$
;  ptr_new(map13CO3,/No_copy),ptr_new(map13CO4,/No_copy),$
;  ptr_new(map13CO5,/No_copy),ptr_new(map13CO6,/No_copy)]
;snrpk13CO=[ptr_new(snrpk13CO1,/No_copy),ptr_new(snrpk13CO2,/No_copy),$
;  ptr_new(snrpk13CO3,/No_copy),ptr_new(snrpk13CO4,/No_copy),$
;  ptr_new(snrpk13CO5,/No_copy),ptr_new(snrpk13CO6,/No_copy)]
;;
;n=0
;tau0result=[0]
;tau1result=[0]
;tau5result=[0]
;tau10result=[0]
;tau100result=[0]
;for i=0,5 do begin
;  xyad,*map13COhead(i),0,0,a0,d0
;  adxy,COspechead,a0,d0,x0,y0
;  z12CO=((dindgen(fxpar(COspechead,'NAXIS3'))+0.5-fxpar(COspechead,'CRPIX3'))*fxpar(COspechead,'CDELT3')+fxpar(COspechead,'CRVAL3'))/1.0e3
;  z13CO=((dindgen(fxpar(*map13COhead(i),'NAXIS3'))+0.5-fxpar(*map13COhead(i),'CRPIX3'))*fxpar(*map13COhead(i),'CDELT3')+fxpar(*map13COhead(i),'CRVAL3'))/1.0e3
;  for x=0,(size(*map13CO(i)))(1)-1 do begin
;    for y=0,(size(*map13CO(i)))(2)-1 do begin
;      if (*snrpk13CO(i))(x,y) ge 3 && COmasksum(x+x0,y+y0) ge 5 then begin
;        spec12CO=reform(COspec(x+x0,y+y0,*))
;        spec13CO=reform((*map13CO(i))(x,y,*))
;        result12CO=gaussfit(z12CO,spec12CO,gspara12CO,nterms=3)
;        result13CO=gaussfit(z13CO,spec13CO,gspara13CO,nterms=3)
;        if gspara13CO(1) gt gspara12CO(1)-3*gspara12CO(2) && gspara13CO(1) lt gspara12CO(1)+3*gspara12CO(2) && gspara12CO(2) lt 4 then begin
;;          window,n mod 4
;;          plot,z12CO,spec12CO
;;          oplot,z13CO,spec13CO,color=255
;;          print,gspara12CO,gspara13CO
;          deriveCOtau,gspara12CO(0),gspara13CO(0),tau_0,0,Tex12CO0
;;          deriveCOtau,gspara12CO(0),gspara13CO(0),tau1,1,Tex12CO1
;;          deriveCOtau,gspara12CO(0),gspara13CO(0),tau5,5,Tex12CO5
;;          deriveCOtau,gspara12CO(0),gspara13CO(0),tau10,10,Tex12CO10
;          deriveCOtau,gspara12CO(0),gspara13CO(0),tau100,100,Tex12CO100
;          tau0result=[tau0result,tau_0]
;;          tau1result=[tau1result,tau1]
;;          tau5result=[tau5result,tau5]
;;          tau10result=[tau10result,tau10]
;          tau100result=[tau100result,tau100]
;          n++
;        endif
;      endif
;    endfor
;  endfor
;endfor
;;
;cgdisplay,800,600
;tau0result=tau0result(1:(size(tau0result))(1)-1)
;;tau1result=tau1result(1:(size(tau1result))(1)-1)
;;tau5result=tau5result(1:(size(tau5result))(1)-1)
;;tau10result=tau10result(1:(size(tau10result))(1)-1)
;tau100result=tau100result(1:(size(tau100result))(1)-1)
;;print,min(tauresult),max(tauresult)
;;
;save,tau0result,filename='/mnt/science1/bliu/HINSA/data/2019result/tau0.sav'
;save,tau100result,filename='/mnt/science1/bliu/HINSA/data/2019result/tau100.sav'

  restore,'/mnt/science1/bliu/HINSA/data/2019result/tau0.sav'
  restore,'/mnt/science1/bliu/HINSA/data/2019result/tau100.sav'

  set_plot,'PS'
  filename='/mnt/science1/bliu/HINSA/code/20190813/'+'calcCOtau.eps'
  Device,/ENCAPSUL, BITS_PER_PIXEL=8, /COLOR, FILENAME=filename, XSIZE=26, YSIZE=20
  !P.Color = '000000'xL
  !P.Background = 'FFFFFF'xL
  loadct,39,ncolors=255,bottom=1,/silent

  mi=0
  ma=20
  bin=0.5
  cgplot,indgen((ma-mi)/bin)*bin+mi,histogram(tau100result,binsize=bin,max=ma,min=mi),yrange=[0,60],psym=10,xstyle=1,xtitle=textoidl('\tau(^{12}CO)'),ytitle='counts',charsize=1.5
  ;oplot,indgen((ma-mi)/bin)*bin+mi,histogram(tau1result,binsize=bin,max=ma,min=mi),psym=10,color=150
  ;oplot,indgen((ma-mi)/bin)*bin+mi,histogram(tau5result,binsize=bin,max=ma,min=mi),psym=10,color=180
  ;oplot,indgen((ma-mi)/bin)*bin+mi,histogram(tau10result,binsize=bin,max=ma,min=mi),psym=10,color=210
  ;oplot,indgen((ma-mi)/bin)*bin+mi,histogram(tau100result,binsize=bin,max=ma,min=mi),psym=10,color=240

result0=where(tau100result gt 0,count0)
result1=where(tau100result(result0) lt 1,count1)
result2=where(tau100result(result0) lt 2,count2)
;print,'number of pixels :',n,',number of pixels with tau12CO > 0 :',count0
print,'number and fraction of pixels where tau12CO < 1:',count1,count1*1.0/count0
print,'number and fraction of pixels where tau12CO < 2:',count2,count2*1.0/count0
print,'average tau12CO',total(tau100result(result0))/count0
print,'MEDIAN tau12CO',MEDIAN(tau100result(result0))

  Device, /close

  print,'finished'

END

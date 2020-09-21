PRO DR3printtable

xarray=[1378,1477,1428,1269,1147,1088,579,571,613,677,731,906,666,1447,1396,671,673,1073,728,913,686,593,1491,949,1441,1101,1096,723,1246,1551,1412,1428,884,1555,1489,1486,901,1395,1096,952,916,595]
yarray=[1601,1453,1018,1404,1131,1046,781,985,1029,971,1074,998,761,988,1052,955,997,1157,1106,1221,1104,1018,1144,1304,1249,1387,1448,1366,1479,1386,1510,1535,1609,1387,1147,1163,997,1056,1451,1300,1220,1021]

tau0file='/mnt/science1/bliu/HINSA/code/20190813/DR3tau0cube.fits'
THfile='/mnt/science1/bliu/HINSA/code/20190813/DR3THcube.fits'
sigmaHfile='/mnt/science1/bliu/HINSA/code/20190813/DR3sigmaHcube.fits'

tau0cube=readfits(tau0file,tau0head)
THcube=readfits(THfile,THhead)
sigmaHcube=readfits(sigmaHfile,sigmaHhead)

restore,'/mnt/science1/bliu/HINSA/code/20190813/NHIcube.sav',/verbose
restore,'/mnt/science1/bliu/HINSA/code/20190813/NH2cube.sav',/verbose
restore,'/mnt/science1/bliu/HINSA/code/20190813/ratioprime.sav',/verbose

print,size(ratioprime)

for m=0,(size(xarray))(1)-1 do begin
    x=xarray(m)
    y=yarray(m)
    xyad,tau0head,x,y,ra,dec
    tau0=tau0cube(x,y)
    TH=THcube(x,y)
    sigmaH=sigmaHcube(x,y)
    NHI=NHIcube(x,y)
    NH2=NH2cube(x,y)
    ratio=ratioprime(x,y)
    print,m+1,ra,dec,tau0,TH,sigmaH,NHI,NH2,ratio
endfor

print,'finished'
END
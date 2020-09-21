PRO DR3HItoH2_ringsplot

radius=[0.25,0.75,1.25,1.75,2.25,2.75,3.25]
ratio=[-3.05316,-2.90043,-2.89712,-2.91288,-2.96242,-2.75711,-2.93326]
sigma=[0.468707,0.402636,0.46669,0.350117,0.370984,0.37023,0.360335]

set_plot,'PS'
filename='/mnt/science1/bliu/HINSA/code/20190813/'+'DR3HItoH2_rings.eps'
Device,/ENCAPSUL, BITS_PER_PIXEL=8, /COLOR, FILENAME=filename, XSIZE=24, YSIZE=16
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct,39,ncolors=255,bottom=1,/silent

xtitle='Radius(kpc)'
ytitle=textoidl('HINSA HI to H2 ratio[LOG(N_{HINSA}/N_{H_{2}})]')

cgPlot, radius, ratio, Color='black', PSym=-16, SymColor='black', $
  SymSize=1.0, Thick=thick, XTitle=xtitle, YTitle=ytitle, $
  Position=position, YRange=[-6,0], XRange=[0,3.5], xstyle=1, YStyle=1, $
  ERR_YLow=sigma, ERR_YHigh=sigma, ERR_Color='black'

Device, /close
print,'finished'

END

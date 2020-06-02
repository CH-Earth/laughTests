pro colbeck1976contour_testSUMMA

; define plotting parameters
window, 0, xs=1200, ys=1000, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=5
!P.COLOR=0
erase, color=255
!P.MULTI=[0,2,3,0,0]

; define the date format
dummy = label_date(date_format=['%H!C%I'])

; define the summa path
root = '/Users/mac414/summaTestCases_2.1/'

; define the file path
file_path = root + 'output/syntheticTestCases/colbeck1976/'

; define the path to the graphics file
gpath = root + 'plottingTestCases/zFigures/'

; define the name of the graphics file
gname = 'contourPlot_colbeck1976.png'

; define the plot range
vmin = 0.00d
vmax = 0.2d
 
; *****
; * ANALYTICAL SOLUTION...
; ************************

; define constants
gravity      =       9.80616d  ; acceleration of gravity              (m s-2)
Cp_air       =    1005.d       ; specific heat of air                 (J kg-1 K-1)
Cp_ice       =    2114.d       ; specific heat of ice                 (J kg-1 K-1)
Cp_soil      =     850.d       ; specific heat of soil                (J kg-1 K-1)
Cp_water     =    4181.d       ; specific heat of liquid water        (J kg-1 K-1)
dynVisc      =       0.001781d ; dynamic viscosity of water           (kg m-1 s-1)
Tfreeze      =     273.16d     ; temperature at freezing              (K)
LH_fus       =  333700.0d      ; latent heat of fusion                (J kg-1)
LH_vap       = 2501000.0d      ; latent heat of vaporization          (J kg-1)
LH_sub       = 2834700.0d      ; latent heat of sublimation           (J kg-1)
iden_air     =       1.293d    ; intrinsic density of air             (kg m-3)
iden_ice     =     917.0d      ; intrinsic density of ice             (kg m-3)
iden_water   =    1000.0d      ; intrinsic density of liquid water    (kg m-3)
secprday     =   86400.d       ; number of seconds in a day

; define domain
depth       = 1.d   ; snow depth (meters)
snowDensity = 300.d ; snow density (kg m-3)

; define flux
rainfall  = 10.d^(-5.d)   ; rainfall (m s-1)
duration  = 10800.d       ; duration of rainfall flux (seconds)
startRain = 0.d           ; start of rainfall

; define exponent in flux parameterization
n = 3.d

; define initial relative water content
Sw = [0.07d, 0.d, 0.d]

; define porosity
porosityVec = (snowDensity - iden_ice) / (iden_water*Sw - iden_ice)
print, 'porosity = ', porosityVec

; define the volumetric fraction of liquid water
volFracLiq = Sw*porosityVec
print, 'volFracLiq = ', volFracLiq

; define the volumetric fraction of ice
volFracIce = 1.d - porosityVec
print, 'volFracIce = ', volFracIce

; compute fraction of liquid water
volFracTot = volFracIce*iden_ice/iden_water + volFracLiq
fracLiquid = volFracLiq/volFracTot

; compute temperature based on the fraction of liquid water (K)
fc_param = 50.d
ripeTemp = Tfreeze - ((1.d/fracliquid[0] - 1.d)/fc_param^2.d)^(0.5d)
print, 'ripeTemp = ', ripeTemp

; refine residual liquid water content
residLiqVec = Sw[0]*porosityVec

; define experiments
expName   = ['Ripe snow', 'Refrozen snow', 'Fresh snow']    ; experiment name
grainSize = [ 2.0d,       2.0d,  0.2d] / 1000.d             ; grain size (m)
snowTemp  = [ripeTemp, -5.d + Tfreeze, -5.d + Tfreeze]      ; temperature (K)
thetaLiq  = volFracLiq                                      ; initial volumetric fraction of liquid water (-)

; loop through experiments
for iExp=0,n_elements(volFracLiq)-1 do begin

 ; define the filename
 file_name = file_path + 'colbeck1976-exp'+strtrim(iExp+1,2)+'_output_timestep.nc'
 
 ; print the experiment name
 print, expName[iExp]
 
 ; get porosity and residual liquid water content
 porosity = porosityVec[iExp]
 residLiq = residLiqVec[iExp]
 
 ; compute hydraulic conductivity (m s-1)
 permeability = grainSize[iExp]^2.d * 0.077d*exp(-7.8d*snowDensity/iden_water)  ; permeability (m2)
 conductivity = permeability*iden_water*gravity/dynVisc   ; hydraulic conductivity (m s-1)
 print, 'ratio           = ', permeability/grainSize[iExp]^2.d
 print, 'grainsize^2     = ', grainSize[iExp]^2.d
 print, 'permeability    = ', permeability
 print, 'conductivity    = ', conductivity
 
 ; define the volumetric fraction of water required to maintain the speed of the wetting front
 satRequire   = (rainfall/conductivity)^(1.d/n)
 thetaRequire = residLiq + satRequire*(porosity - residLiq) - thetaLiq[iExp]
 print, 'residLiq        = ', residLiq
 print, 'porosity        = ', porosity
 print, 'satRequire      = ', satRequire
 print, 'thetaRequire    = ', thetaRequire
 
 ; define the level of saturation necessary to satisfy thermal requirements
 heatRequire = (snowTemp[iExp] - Tfreeze)*(porosity - 1.d) * iden_ice*cp_ice/(iden_water*LH_fus)
 print, 'heatRequire     = ', heatRequire
 
 ; compute the time required for the wetting front to reach the bottom of the snowpack
 watRequire      = depth*(thetaRequire + heatRequire)  ; water required so the entire profile = thetaRequire (m)
 timeWetting     = watRequire/rainfall
 celerityWetting = depth/timeWetting
 print, 'watRequire      = ', watRequire
 print, 'timeWetting     = ', timeWetting
 print, 'timeWetting (h) = ', timeWetting/3600.d
 print, 'celerityWetting = ', celerityWetting
 
 ; compute the celerity (m s-1)
 celerity = (n/(porosity - residLiq)) * conductivity^(1.d/n) * rainfall^((n - 1.d)/n)
 print, 'celerity        = ', celerity
 
 ; compute the lag time (s)
 timeLag    = depth/celerity
 timeDrying = duration + timeLag
 print, 'lag time        = ', timeLag
 
 ; check if the snowpack reaches steady state
 if(timeWetting lt timeDrying)then begin
 
  ; no shock exists
  timeShock  = 10.d * 3600.d  ; set the shock time to the end of the simulation (10 hours)
  timeBottom = timeWetting    ; water reaches the bottom of the snowpack at timeWetting
 
 ; check if the drying front overtakes the wetting front (kinematic shock)
 endif else begin
 
  ; define the depth of the kinematic shock (the depth that the lines intersect)
  ; time is equal, so: depth/celerityWetting + startRain = depth/celerity + duration
  ; then: depth/celerityWetting - depth/celerity = duration - startRain
  ; so
  depthShock = (duration - startRain) / (1.d/celerityWetting - 1.d/celerity)
  print, 'depthShock = ', depthShock
 
  ; get the time when the wetting front stops propagating at its maximum value
  timeShock  = depthShock/celerity + duration
  print, 'timeShock  = ', timeShock/3600.d
 
  ; define short-cut variables
  a = conductivity*(depthShock*(porosity - residLiq)/(n*conductivity))^(n/(n - 1.d))
  b = (porosity - residLiq)*( ((porosity - residLiq)/(n*conductivity))^(1.d/(n - 1.d)) )
  c = residLiq + heatRequire
 
  ; time that shock reaches the bottom of the snowpack
  tmp1       = a*(1.d - n) - b*((n - 1.d)/n)*( depth^(n/(n - 1.d)) - depthShock^(n/(n - 1.d)) )
  tmp2       = a*(1.d - n)*(timeShock - duration)^(1.d/(1.d - n)) + c*(depth - depthShock)
  tb         = (tmp2/tmp1)^(1.d - n)
  timeBottom = tb + duration
  print, 'tb = ', tb + duration
  print, 'tb (hours) = ', (tb + duration)/3600.d

 endelse  ; if kinematic shock
 
 ; *****
 ; * COMPUTE ANALYTICAL SOLUTIONS FOR THE TIME SERIES OF VOLUMETRIC LIQUID WATER PROFILE...
 ; ****************************************************************************************

 ; define the volumetric liquid water profile for multiple times
 nTime = 60*10    ; every minute for 10 hours
 xTime = 3600.d*10.d*(dindgen(nTime)+1.d)/double(nTime)  ; time in seconds

 ; define depths
 nDepth  = 100
 yDepth  = depth*(dindgen(nDepth)+0.5d)/double(nDepth)
 yHeight = depth*(dindgen(nDepth+1))/double(nDepth)

 ; define an array for the solution
 volFracLiq = dblarr(nDepth,nTime)

 ; loop through times
 for iTime=0,nTime-1 do begin

  ; identify specific cases
  case 1 of

   ; case 1: wetting front while rain is still falling
   xTime[iTime] lt timeWetting and xTime[iTime] lt duration: begin
    depthWetting = (xTime[iTime] - startRain)*celerityWetting
    depthDrying  = 0.d
   end

   ; case 2: wetting front while rain is still falling
   xTime[iTime] ge timeWetting and xTime[iTime] lt duration: begin
    depthWetting = depth 
    depthDrying  = 0.d
   end

   ; case 3: drying front above wetting front
   xTime[iTime] ge duration and xTime[iTime] lt min([timeShock,timeDrying]): begin
    depthWetting = (xTime[iTime] - startRain)*celerityWetting
    depthDrying  = (xTime[iTime] - duration)*celerity
   end

   ; case 4: propagation of the kinematic shock
   xTime[iTime] ge timeShock and xTime[iTime] lt timeBottom: begin
    depthWetting = call_function('getDepth', xTime[iTime], duration, timeShock, depthShock, depth, a, b, c, n)
    depthDrying  = depthWetting
   end

   ; case 5: water reached the bottom of the snowpack
   xTime[iTime] ge timeBottom and xTime[iTime] ge duration: begin
    depthWetting = depth
    depthDrying  = depth
   end

   ; check that we found everything
   else: stop, 'could not find case'

  endcase

  ; define areas of the snowpack
  iDry = where(yDepth gt depthWetting,                           nDry)  ; no change from intial state
  iSat = where(yDepth le depthWetting and yDepth gt depthDrying, nSat)  ; maximum possible wetting
  iMid = where(yDepth lt depthDrying,                            nMid)  ; drying phase
  if(nDry+nSat+nMid ne nDepth)then stop, 'problem identifying different areas of the snowpack'

  ; define volumetric fraction of liquid water
  if(nDry gt 0)then volFracLiq[iDry,iTime] = thetaLiq[iExp]                ; initial state
  if(nSat gt 0)then volFracLiq[iSat,iTime] = thetaRequire + thetaLiq[iExp] ; maximum possible wetting
  if(nMid gt 0)then begin
   satProfile             = ( (porosity - residLiq)*yDepth[iMid]/(n*conductivity*(xTime[iTime] - duration)) )^(1.d/(n - 1.d))
   volFracLiq[iMid,iTime] = residLiq + (porosity - residLiq)*satProfile
  endif

 endfor   ; looping through time

 ; *****
 ; * NUMERICAL SOLUTION...
 ; ***********************
 
 ; define the HRU
 iHRU=0

 ; define the variable name
 summaVarName = 'mLayerVolFracLiq'
 
 ; open files
 ncFileID = ncdf_open(file_name, /nowrite)
 
  ; get time units
  ivar_id = ncdf_varid(ncFileID,'time')
  ncdf_attget, ncFileID, ivar_id, 'units', bunits
  cunits = string(bunits)
 
  ; extract the units "words"
  tunit_words = strsplit(string(cunits),' ',/extract)
  tunit_idate = fix(strsplit(tunit_words[2],'-',/extract))
  tunit_ihour = fix(strsplit(tunit_words[3],':',/extract))
  bjulian     = julday(tunit_idate[1],tunit_idate[2],tunit_idate[0],tunit_ihour[0],tunit_ihour[1],tunit_ihour[2])
 
  ; get the offset in days
  if(strtrim(tunit_words[0],2) eq 'seconds') then aoff=1.d/86400.d else stop, 'unknown time units'
 
  ; extract the time vector
  ncdf_varget, ncFileID, ivar_id, atime
  djulian_summa = bjulian + atime*aoff
 
  ; extract the number of snow layers
  ivar_id = ncdf_varid(ncFileID,'nSnow')
  ncdf_varget, ncFileID, ivar_id, xVar
  nSnow=reform(xVar)
 
  ; extract the total number of layers
  ivar_id = ncdf_varid(ncFileID,'nLayers')
  ncdf_varget, ncFileID, ivar_id, xVar
  nLayers=reform(xVar)
 
  ; get the number time elements
  nTime = n_elements(djulian_summa)
 
  ; get the matrix for the desired variable
  summaVar = fltarr(max(nSnow),nTime)
 
  ; get the matrix for the coordinate variable
  snowHeight = fltarr(max(nSnow)+1,nTime)
 
  ; loop through time
  for iTime=0,nTime-1 do begin
 
   ; get the height
   ivar_id = ncdf_varid(ncFileID,'iLayerHeight')
   ncdf_varget, ncFileID, ivar_id, xVar, offset=[iHRU,0,iTime], count=[1,nSnow[itime]+1,1]
   snowHeight[0:nSnow[iTime],itime] = reform(xVar)
 
   ; get the desired variable
   ivar_id = ncdf_varid(ncFileID,summaVarName)
   ncdf_varget, ncFileID, ivar_id, xVar, offset=[iHRU,0,iTime], count=[1,nSnow[itime],1] 
   summaVar[*,iTime] = reform(xVar)
 
  endfor  ; looping through time
 
  ; extract the outflow from the snowpack
  ivar_id = ncdf_varid(ncFileID,'scalarRainPlusMelt')
  ncdf_varget, ncFileID, ivar_id, xVar
  summaFlux=reform(xVar)
 
 ncdf_close, ncFileID

 ; *****
 ; * COMMON PLOTTING INFORMATION...
 ; ********************************

 ; define x title
 if(iExp eq 2)then xtit='Time (hours)' else xtit=' '

 ; define margins
 if(iExp eq 0)then ymar=[0,4]
 if(iExp eq 1)then ymar=[2,2]
 if(iExp eq 2)then ymar=[4,0]
 
 ; *****
 ; * PLOT ANALYTICAL CALCULATIONS...
 ; *********************************

 ; define plot title
 if(iExp eq 0)then ptit='Analytical solution' else ptit=' ' 

 ; make a base plot
 plot, indgen(5), xrange=[0,10], yrange=[1,0], xstyle=9, ystyle=1, $
  xticklen=(-0.02), xtitle=xtit, xmargin=[10,4], ymargin=ymar, ytitle='Depth (m)', title=ptit, /nodata

 ; make a hovmuller plot
 variableHeight=0
 time_hovmuller, xTime/3600.d, yHeight, volFracLiq, vmin, vmax, variableHeight

 ; re-draw the axes
 plots, [0,10], [0,0]
 plots, [10,10], [1,0]

 ; label the plots
 xyouts, -3, 0.5, '!17'+expName[iExp]+'!3', alignment=0.5, orientation=90, charsize=3

 ; *****
 ; * PLOT NUMERICAL SIMULATIONS...
 ; *******************************

 ; define plot title
 if(iExp eq 0)then ptit='Model simulations' else ptit=' ' 

 ; define time
 time0 = min(djulian_summa) - (djulian_summa[1] - djulian_summa[0])
 time1 = max(djulian_summa)
 
 ; make a base plot
 plot, indgen(5), xrange=[time0,time1], yrange=[0,-1], xstyle=13, ystyle=13, $
  xticks=1, xtickname=[' ',' '], xmargin=[0,14], ymargin=ymar, title=ptit, /nodata
 
 ; make a hovmuller plot
 variableHeight=1
 time_hovmuller, djulian_summa, snowHeight, summaVar, vmin, vmax, variableHeight
 
 ; re-draw the axes
 plots, [time0,time1], [-1,-1]
 plots, [time1,time1], [-1,0]
 axis, yaxis=0, yrange=[1,0], /save
 axis, xaxis=0, xrange=[0,10], xticklen=(-0.02), xtitle=xtit, /save

endfor  ; looping through experiments

; plot colorbar
nbins   = 201
ccolors = 40.d + 180.d*dindgen(nbins)/double(nbins-1)
clabels = vmin + (vmax - vmin)*dindgen(nbins)/double(nbins-1)
colorbar, 0.855, 0.90, 0.1, 0.9, clabels, ccolors, every=20, charsize=2, /nobox, /norm
xyouts, 0.97, 0.5, 'Volumetric liquid water content (-)', alignment=0.5, orientation=90, charsize=3, /normal
 
; write figure
write_png, gpath+gname, tvrd(true=1)
 
 
stop
end

; ========================================================================================
; ========================================================================================
; ========================================================================================
; ========================================================================================
; ========================================================================================

; *****
; * make a hovmuller diagram for the day-time...
; **********************************************

pro time_hovmuller, dTime, simHeight, simVar, vmin, vmax, variableHeight

; get the dates
xSmall = 1.d-6
caldat, dTime-xSmall, im, id, iyyy, ih, imin, asec

; get the time step
dt = dtime[1] - dtime[0]

; loop through the time series
for iTime=0,n_elements(dTime)-1 do begin

 ; get the julian day
 djulian = julday(im[itime], id[itime], iyyy[itime], ih[itime], imin[itime])

 ; get the height
 if(variableHeight eq 1)then height=reform(simHeight[*,itime]) else height=simHeight

 ; loop through the depths
 for iDepth=0,n_elements(simVar[*,itime])-1 do begin 

  ; identify the color
  if(simVar[iDepth,iTime] gt -9998.)then begin
   icolor = ( (simVar[iDepth,iTime] - vmin) / (vmax - vmin) )*180.d + 40.d
   if(simVar[iDepth,iTime] lt vmin)then icolor=40
   if(simVar[iDepth,iTime] gt vmax)then icolor=250
  endif else begin
   icolor = 255
  endelse

  ; plot data
  x0 = djulian
  x1 = djulian + dt
  y0 = height[iDepth]
  y1 = height[iDepth+1]
  polyfill, [x0,x1,x1,x0], [y0,y0,y1,y1], color=icolor

 endfor  ; layers
endfor  ; time

end

; *****
; * compute depth for a given time after the kinematic shock...
; *************************************************************

; use bi-section method to compute depth for a given time after the kinematic shock
function getDepth, xTime, timeCenter, timeShock, depthShock, depthTotal, a, b, c, n

; return value:
; depth of the kinematic shock at time xTime

; inputs:
; xTime      = desired time
; timeCenter = center time of the characteristics, the time that rainfall ends (timeShock > timeEnd)
; timeShock  = time of the kinematic shock (xTime > timeShock)
; depthShock = depth of the kinematic shock
; depthTotal = total depth of snow
; a, b, c    = constants
; n          = exponent

; define the print flag
print=0  ; 1=print

; define the number of iterations
niter = 100

; define the convergence criteria
convTol = 0.000001d

; define brackets
d0 = depthTotal
d1 = depthShock

; iterate
for iter=1,niter do begin

 ; get midpoint
 depthTrial = 0.5d*(d0 + d1)

 ; trial solution
 tmp1 = a*(1.d - n) - b*((n - 1.d)/n)*( depthTrial^(n/(n - 1.d)) - depthShock^(n/(n - 1.d)) )
 tmp2 = a*(1.d - n)*(timeShock - timeCenter)^(1.d/(1.d - n)) + c*(depthTrial - depthShock)
 tt   = (tmp2/tmp1)^(1.d - n) + timeCenter

 ; function evaluation
 if(abs(tt - xTime) lt convTol)then break
 if(iter eq niter)then stop, 'convergence problem in getDepth'

 ; print
 if(print eq 1)then print, 'iter, d0, d1, depthTrial, tt, xTime = ', iter, d0, d1, depthTrial, tt, xTime, format='(a,1x,i4,1x,10(f15.8,1x))'

 ; update bounds
 if(xTime lt tt)then d0=depthTrial
 if(xTime gt tt)then d1=depthTrial

endfor  ; iterating

; return depth
return, depthTrial

end




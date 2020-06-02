pro colbeck1976timeseries_testSUMMA

; define plotting parameters
window, 0, xs=1000, ys=1000, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=4
!P.COLOR=0
erase, color=255
!P.MULTI=[0,1,3,0,0]

; define the summa path
root = '/Users/mac414/summaTestCases_2.1/'

; define the file path
file_path = root + 'output/syntheticTestCases/colbeck1976/'

; define the path to the graphics file
gpath = root + 'plottingTestCases/zFigures/'

; define the name of the graphics file
gname = 'timeseriesPlot_colbeck1976.png'

; define the date format
dummy = label_date(date_format=['%H!C%I'])

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
for iExp=0,2 do begin

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
 
  ; define the time the water reaches the bottom of the snowpack
  timeShock = timeDrying
 
 ; check if the drying front overtakes the wetting front (kinematic shock)
 endif else begin
 
  ; define the depth of the kinematic shock (the depth that the lines intersect)
  ; time is equal, so: depth/celerityWetting + startRain = depth/celerity + duration
  ; then: depth/celerityWetting - depth/celerity = duration - startRain
  ; so
  depthWet = (duration - startRain) / (1.d/celerityWetting - 1.d/celerity)
  print, 'depthWet = ', depthWet
 
  ; get the time when the wetting front stops propagating at its maximum value
  timeWet  = depthWet/celerity + duration
  print, 'timeWet  = ', timeWet/3600.d
 
  ; define short-cut variables
  a = conductivity*(depthWet*(porosity - residLiq)/(n*conductivity))^(n/(n - 1.d))
  b = (porosity - residLiq)*( ((porosity - residLiq)/(n*conductivity))^(1.d/(n - 1.d)) )
 
  ; time that shock reaches the bottom of the snowpack
  tmp1      = a*(1.d - n) - b*((n - 1.d)/n)*( depth^(n/(n - 1.d)) - depthWet^(n/(n - 1.d)) )
  tmp2      = a*(1.d - n)*(timeWet - duration)^(1.d/(1.d - n)) + (residLiq + heatRequire)*(depth - depthWet)
  tb        = (tmp2/tmp1)^(1.d - n)
  timeShock = tb + duration
  print, 'tb = ', tb + duration
  print, 'tb (hours) = ', (tb + duration)/3600.d
 
 endelse  ; if kinematic shock
 
 ; define the falling limb of the hydrograph
 ntime = 100
 time_end = 10.d * 3600.d
 timeFall = timeShock + (time_end - timeShock)*dindgen(ntime)/double(ntime-1)
 fluxFall = conductivity*( (porosity - residLiq)*depth/(n*conductivity*(timeFall - duration)) )^(n/(n - 1.d))
 print, 'timeFall = ', timeFall
 print, 'fluxFall = ', fluxFall
 
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
 ; * PLOT THE HYDROGRAPH...
 ; ************************
 
 ; define time
 time0 = min(djulian_summa) - (djulian_summa[1] - djulian_summa[0])
 time1 = max(djulian_summa)

 ; define x title
 if(iExp eq 2)then xtit='Time (hours)' else xtit=' '

 ; define y margins
 if(iExp eq 0)then ymar=[3,1]
 if(iExp eq 1)then ymar=[4,0]
 if(iExp eq 2)then ymar=[5,-1]
 
 ; make a base plot
 plot, indgen(5), xrange=[time0,time1], yrange=[0,0.011], xstyle=1, ystyle=1, $
  xticks=1, xtickname=[' ', ' '], xtitle=xtit, ytitle='Flux (mm/s)', ymargin=ymar, /nodata
 
 ; plot the summa simulations
 oplot, djulian_summa, summaFlux*1000.d, color=40, thick=3
 
 ; plot another x axis
 axis, xaxis=0, xrange=[0,10], xticklen=(-0.02), /save 
 
 ; plot the experiment
 xyouts, 9.75, 0.0085, expName[iExp], charsize=3, alignment=1

 ; plot the rising hydrograph (when the snowpack reaches steady state)
 if(timeWetting lt timeDrying)then begin
  plots, [0,timeWetting]/3600.d, [0,0]*1000.d, color=80, thick=2, linestyle=2
  plots, [timeWetting,timeWetting]/3600.d, [0,rainfall]*1000.d, color=80, thick=2, linestyle=2
  plots, [timeWetting,timeDrying]/3600.d, [rainfall,rainfall]*1000.d, color=80, thick=2, linestyle=2
 
 ; plot the rising hydrograph (when the snowpack does not reach steady state)
 endif else begin
  plots, [timeShock,timeShock]/3600.d, [0,fluxFall[0]]*1000.d, color=80, thick=2, linestyle=2
 endelse

 ; plot the falling hydrograph
 oplot, timeFall/3600.d, fluxFall*1000.d, color=80, thick=2, linestyle=2

endfor  ; looping through experiments

; make a legend
plots, [0.5,1.5], [0.009, 0.009], color=80, thick=2, linestyle=2
plots, [0.5,1.5], [0.008, 0.008], color=40, thick=3
xyouts, 1.55, 0.009-0.00025, 'Analytical solution', charsize=2
xyouts, 1.55, 0.008-0.00025, 'Numerical simulation', charsize=2

; write figure
write_png, gpath+gname, tvrd(true=1)


stop
end


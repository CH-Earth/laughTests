pro miller1998_testSUMMA

; define plotting parameters
window, 0, xs=1000, ys=500, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=2
!P.COLOR=0
erase, color=255
!P.MULTI=1   ; [0,1,2]

; define the HRU
iHRU=0

; define the maximum soil depth
maxDepth = 10.d

; define the time step
dTime = 900.d

; define the desired time (units of seconds)
xTime = [0.18d, 2.25d, 1.d] * 86400.d

; define the soil depth
sDepth = [10.d,5.d,2.d]

; define the summa path
root = '/Users/mac414/summaTestCases_2.1/'

; define the path to the graphics file
gpath = root + 'plottingTestCases/zFigures/'

; define the name of the graphics file
gname = 'syntheticTestCase_miller1998.png'

; define the file path
fpath = root + 'output/syntheticTestCases/miller1998/'

; define experiment names
expName = ['millerSand','millerLoam','millerClay']

; define the plot name
plotName = ['sand','loam','clay loam']

; define file suffix
fSuff = '_output_timestep.nc'

; define colors
icol = [80,210,150]

; loop through variables (liquid and matric)
for ivar=1,1 do begin

 ; define variables
 if(ivar eq 0)then cvar='mLayerVolFracLiq'
 if(ivar eq 1)then cvar='mLayerMatricHead'

 ; make a base plot
 if(ivar eq 0) then begin
  plot, indgen(5), yrange=[0,maxDepth], xrange=[0.05,0.45], xstyle=1, ystyle=1, $
   ytitle='Height above profile bottom (m)', xtitle='Volumetric liquid water content (-)', /nodata
 endif else begin
  plot, indgen(5), yrange=[0,maxDepth], xrange=[-5,0.5], xstyle=1, ystyle=1, $
   ytitle='Height above profile bottom (m)', xtitle='Pressure head (m)', /nodata
 endelse

 ; loop through files
 for ifile=0,n_elements(expName)-1 do begin
  
  ; define file names
  filenm = expName[ifile] + fSuff

  ; define the desired time step
  nTime = floor(xTime[ifile]/dTime + 0.5d)

  ; open file
  nc_file = ncdf_open(fpath+filenm, /nowrite)

  ; extract the time vector
  ivar_id = ncdf_varid(nc_file,'time')
  ncdf_varget, nc_file, ivar_id, atime

  ; get the number of soil layers
  ivar_id = ncdf_varid(nc_file,'nSoil')
  ncdf_varget, nc_file, ivar_id, nSoil, offset=[iHRU,ntime], count=[1,1]

  ; get the mid-point of each layer
  ivar_id = ncdf_varid(nc_file,'mLayerHeight')
  ncdf_varget, nc_file, ivar_id, mLayerHeight, offset=[iHRU,0,ntime], count=[1,nSoil,1]

  ; get the desired variable for all layers
  ivar_id = ncdf_varid(nc_file,cvar)
  ncdf_varget, nc_file, ivar_id, avar, offset=[iHRU,0,ntime], count=[1,nSoil,1]

  ; plot the data
  ;oplot, avar[0,*], reverse(reform(mLayerHeight[0,*])), color=icol[ifile], psym=sym(1), symsize=1
  oplot, avar[0,*], reverse(reform(mLayerHeight[0,*])), color=icol[ifile], thick=3
  print, 'time = ', long(ntime+1)*long(dTime)
  print, reform(avar[0,*])

  ; plot the legend
  ypos = 9.5d - 2.d*double(ifile+1)/double(n_elements(expName))
  plots, [-4.5,-4], [ypos, ypos], color=icol[ifile], thick=3
  xyouts, -3.95, ypos-0.1, plotName[ifile]

  ; close NetCDF file
  ncdf_close, nc_file

 endfor  ; loop through files
endfor  ; loop through variables

; make a figure
write_png, gpath+gname, tvrd(true=1)

stop
end

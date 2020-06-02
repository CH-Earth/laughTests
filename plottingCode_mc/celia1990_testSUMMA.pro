pro celia1990_testSUMMA

; define plotting parameters
window, 0, xs=1000, ys=700, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=2
!P.COLOR=0
erase, color=255
!P.MULTI=1   ; [0,1,2]

; define the summa path
root = '/Users/mac414/summaTestCases_2.1/'

; define the output file path
fpath = root + 'output/syntheticTestCases/celia1990/'

; define the output file name
fname = 'celia1990_output_timestep.nc'

; define the path to the graphics file
gpath = root + 'plottingTestCases/zFigures/'

; define the name of the graphics file
gname = 'syntheticTestCase_celia1990.png'

; define the HRU index
iHRU=0

; define desired times
desiredTimes=[10,32,49]

; loop through variables (liquid and matric)
for ivar=1,1 do begin

 ; define variables
 if(ivar eq 0)then cvar='mLayerVolFracLiq'
 if(ivar eq 1)then cvar='mLayerMatricHead'

 ; make a base plot
 if(ivar eq 0) then begin
  plot, indgen(5), yrange=[0.6,0], xrange=[0.1,0.22], xstyle=1, ystyle=1, $
   ytitle='Depth (m)', xtitle='Volumetric liquid water content (-)', /nodata
 endif else begin
  plot, indgen(5), yrange=[0.6,0], xrange=[-10,0], xstyle=1, ystyle=1, $
   ytitle='Depth (m)', xtitle='Pressure head (m)', /nodata
 endelse

 ; open file
 nc_file = ncdf_open(fpath+fname, /nowrite)

 ; loop through time
 for jtime=0,n_elements(desiredTimes)-1 do begin

  ; define time index
  itime = desiredTimes[jtime]

  ; define color
  icol = fix((float(itime+1)/52.)*200.) + 50

  ; get the number of soil layers
  ivar_id = ncdf_varid(nc_file,'nSoil')
  ncdf_varget, nc_file, ivar_id, nSoil, offset=[iHRU,itime], count=[1,1]

  ; get the mid-point of each layer
  ivar_id = ncdf_varid(nc_file,'mLayerHeight')
  ncdf_varget, nc_file, ivar_id, mLayerHeight, offset=[iHRU,0,itime], count=[1,nSoil[0],1]

  ; get the desired variable for all layers
  ivar_id = ncdf_varid(nc_file,cvar)
  ncdf_varget, nc_file, ivar_id, avar, offset=[iHRU,0,itime], count=[1,nSoil[0],1]

  ; plot the data
  oplot, avar[0,*], mLayerHeight[0,*], color=icol
  oplot, avar[0,*], mLayerHeight[0,*], color=icol, psym=sym(1)

  ; plot the legend
  ypos = 0.1d*double(jtime+1)/double(n_elements(desiredTimes))
  plots, [-9.5,-8.5], [ypos, ypos], color=icol
  plots, -9, ypos, color=icol, psym=sym(1)
  xyouts, -8.4, ypos+0.005, strtrim(long(itime+1)*1800L,2)+' s'

  ; print
  print, 'time = ', long(itime+1)*1800L, ypos
  print, reform(avar[0,*])

 endfor ; loop through time

 ; close NetCDF file
 ncdf_close, nc_file

endfor  ; loop through variables

; make a figure
write_png, gpath+gname, tvrd(true=1)

stop
end

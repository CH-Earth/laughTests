pro mizoguchi1990_testSUMMA

; define plotting parameters
window, 1, xs=1400, ys=1000, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=5
!P.COLOR=0
erase, color=255
!P.MULTI=[0,3,1,0,0]

; define the summa path
root = '/Users/mac414/summaTestCases_2.1/'

; define the path to the validation data
vpath = root + 'testCases_data/validationData/'

; define the file of the validation data
vfile = 'mizoguchiLabData.txt'

aData = fltarr(4,20)

; read in the mizoguchi data
openr, in_unit, vpath+vfile, /get_lun
 readf, in_unit, aData
free_lun, in_unit

; define constants
Tfreeze = 273.16
iden_ice= 917.

; define the path to the graphics file
gpath = root + 'plottingTestCases/zFigures/'

; define the name of the graphics file
gname = 'syntheticTestCase_mizoguchi1990.png'

; define file path
file_path = root + 'output/syntheticTestCases/mizoguchi1990/'

; define file prefix
file_pref = 'mizoguchi1990_output_timestep'

; define suffixes
cSuffix = ['']

; define the number of HRUs
nHRU=1

; define time steps for output
ixDesire = [720,1440,3000]

; define the y title
ytit = ['Depth (m)', ' ', ' ']

; define the y margin
xmar0 = [10,6,2]
xmar1 = [ 2,6,10]

; define y legend
yleg = ['12 hours','24 hours','50 hours']

for iSuffix=0,0 do begin

 ; *****
 ; * PLOT THE MODEL SIMULATIONS...
 ; *******************************

 ; declare file
 filenm = file_path + file_pref + cSuffix[iSuffix]+'.nc'

 ; open file
 nc_file = ncdf_open(filenm, /nowrite)
 
  ; get time units
  ivar_id = ncdf_varid(nc_file,'time')
  ncdf_attget, nc_file, ivar_id, 'units', bunits
  cunits = string(bunits)
 
  ; extract the units "words"
  tunit_words = strsplit(string(cunits),' ',/extract)
  tunit_idate = fix(strsplit(tunit_words[2],'-',/extract))
  tunit_ihour = fix(strsplit(tunit_words[3],':',/extract))
  bjulian     = julday(tunit_idate[1],tunit_idate[2],tunit_idate[0],tunit_ihour[0],tunit_ihour[1],tunit_ihour[2])
 
  ; get the offset in days
  if(strtrim(tunit_words[0],2) eq 'seconds') then aoff=1.d/86400.d else stop, 'unknown time units'
 
  ; extract the time vector
  ncdf_varget, nc_file, ivar_id, atime
  djulian = bjulian + atime*aoff
 
  ; define the date format
  dummy = label_date(date_format=['%D!C%M'])
 
  ; get the number of time elements
  ntime = n_elements(djulian)-1

  ; loop through the HRUs
  for iHRU=0,nHRU-1 do begin

   ; get the number of soil layers
   ivar_id = ncdf_varid(nc_file,'nSoil')
   ncdf_varget, nc_file, ivar_id, nSoil, offset=[iHRU,0], count=[1,ntime]
 
   ; loop through time
   for jtime=0,n_elements(ixDesire)-1 do begin

    ; make a base plot
    plot, indgen(5), xrange=[0.,0.5], yrange=[0.205,-0.005], xstyle=1, ystyle=1, $
     xtitle='Total volumetric water (-)', ytitle = ytit[jtime], xmargin=[xmar0[jtime],xmar1[jtime]], $
     /nodata

    ; make the plot title
    xyouts, 0.025, 0.005, yleg[jtime], charsize=3

    ; make a legend
    plots, [0.025, 0.125], [0.015, 0.015], color=60, thick=2
    xyouts, 0.127, 0.0165, 'SUMMA', charsize=2

    plots, 0.075, 0.02, psym=sym(1), symsize=2, color=250
    xyouts, 0.127, 0.0215, 'Observed', charsize=2

 
    ; get the time index
    itime = ixDesire[jtime]-1

    ; check time
    print, 'elapsed time = ', (60.d + atime[itime]-atime[0])/3600.d

    ; get the height of the layer mid-points
    ivar_id = ncdf_varid(nc_file,'mLayerHeight')
    ncdf_varget, nc_file, ivar_id, mLayerHeight, offset=[iHRU,0,itime], count=[1,nSoil[itime],1]
 
    ; get the volumetric liquid water content (-)
    ivar_id = ncdf_varid(nc_file, 'mLayerVolFracLiq')
    ncdf_varget, nc_file, ivar_id, mLayerVolFracLiq, offset=[iHRU,0,itime], count=[1,nSoil[itime],1]
 
    ; get the volumetric ice content (-)
    ivar_id = ncdf_varid(nc_file, 'mLayerVolFracIce')
    ncdf_varget, nc_file, ivar_id, mLayerVolFracIce, offset=[iHRU,0,itime], count=[1,nSoil[itime],1] 

    ; get the total water content
    volFracTotal = mLayerVolFracLiq + mLayerVolFracIce

    ; define the color
    icol = 160.*float(jTime)/float(n_elements(ixTime)-1) + 80 

    ; plot the model simulations
    oplot, volFracTotal, mLayerHeight, color=60, thick=2

    ; plot the observations
    oplot, aData[jtime+1,*], -aData[0,*]/1000.d, psym=sym(1), symsize=2, color=250

    endfor  ; looping through time

  endfor  ; looping through HRUs

 ; close the NetCDF file
 ncdf_close, nc_file

endfor ; loop through experiments

; make a figure
write_png, gpath+gname, tvrd(true=1)


stop
end

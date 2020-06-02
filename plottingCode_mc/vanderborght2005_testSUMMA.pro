; *****
; * COLLECTION OF FUNCTIONS....
; *****************************

; function to compute hydraulic conductivity from matric head
function hydCond_psi, psi, ksat, alpha, n, m
 if(psi lt 0.d)then begin
  k = ksat*( ( (1.d - (psi*alpha)^(n-1.d) * (1.d + (psi*alpha)^n)^(-m))^2.d ) / ( (1.d + (psi*alpha)^n)^(m/2.d) ) )
 endif else begin
  k = ksat
 endelse
 return, k
end

; function to compute error in hydraulic conductivity
function solve_psi, psi
 common vGn, qIn, ksat, alpha, n, m
 return, qIn - call_function('hydCond_psi', psi, ksat, alpha, n, m)
end

; function to compute change in x given change in matric head
function integrate_psi, psi0, psi1, qIn, ksat, alpha, n, m
 nInt = 100
 dPsi = (psi1 - psi0)/double(nInt)
 psi  = psi0 + (psi1 - psi0) * (dindgen(nInt)+0.5d)/double(nInt)

 delX = 0.d
 for i=0,nInt-1 do begin
  k    = call_function('hydCond_psi', psi[i], ksat, alpha, n, m)
  delX = delX - dPsi / (qIn/k - 1.d)
  ;delX = - dPsi / (qIn/k - 1.d)
  ;print, 'psi, k, dPsi, delX = ', psi[i], k, dPsi, delX
 endfor
 return, delX
end

; *****
; * MAIN PROCEDURE...
; *******************

; plot the vanderbroght test cases
pro vanderborght2005_testSUMMA

; define common parameters
common vGn, qIn, ksat, alpha, n, m

; define plotting parameters
window, 0, xs=2000, ys=800, retain=2
device, decomposed=0
LOADCT, 39
!P.BACKGROUND=255
!P.CHARSIZE=6
!P.COLOR=0
erase, color=255
!P.MULTI=[0,3,1,0,0]

; define the summa path
root = '/Users/mac414/summaTestCases_2.1/'

; define the path to the graphics file
gpath = root + 'plottingTestCases/zFigures/'

; define the name of the graphics file
gname = 'syntheticTestCase_vanderborght2005.png'

; define file path
fpath = root + 'output/syntheticTestCases/vanderborght2005/'

; define number of experiments
nExp = 3

; define the number of layers
nLayers=100

; define inflow
qIn = 0.5d ; cm d-1

; define the soil types
upperType = ['loam','sand','clay']
lowerType = ['sand','loam','sand']

; define the y titles
ytitle=['depth (m)', ' ', ' ']

; define the x margin
xmar0 = [8, 6, 4]
xmar1 = [0, 2, 4]

; initialize psi for the matric head calculations
psi_init = -25.d

; define the number of values for the curve
nCurve = 1000

; loop thru experiments
for iExp=0,nExp-1 do begin

 ; get the parameters for the lower type
 get_params, lowerType[iExp], ksat, alpha, n, m

 ; compute the matric head based on steady state conditions (outflow = inflow)
 psiLower = newton(psi_init, 'solve_psi')

 ; get the parameters for the upper type
 get_params, upperType[iExp], ksat, alpha, n, m

 ; compute the matric head based on inflow
 psiUpper = newton(psi_init, 'solve_psi')
 
 ; define psi boundaries
 psi = psiLower + (psiUpper - psiLower) * (dindgen(nCurve)+1.d)/double(nCurve)

 ; define xCoord
 xCoord = dblarr(nCurve)
 xCoord[0] = 50.d

 ; define parameters
 nTrial   = 1000
 delX_max = 2.d

 ; loop through points on the curve
 for iCurve=1,nCurve-1 do begin
  
  ; initialize psi limits
  psi0_base = psi[iCurve-1]
  psi1_base = psi[iCurve]

  psi0 = psi0_base
  dPsi = psi1_base - psi0

  ; initialize x coordinate
  xCoord[iCurve] = xCoord[iCurve-1]

  ; loop to constrain delX
  for iTrial=0,nTrial-1 do begin

   ; update psi
   psi1 = psi0 + dPsi

   ; constrain dPsi
   if(dPsi gt 0.d)then begin
    if(psi1 gt psi1_base)then begin
     psi1 = psi1_base
     dPsi = psi1 - psi0
    endif
   endif else begin
    if(psi1 lt psi1_base)then begin
     psi1 = psi1_base
     dPsi = psi1 - psi0
    endif
   endelse

   ; compute delX
   delX = call_function('integrate_psi', psi0, psi1, qIn, ksat, alpha, n, m)

   ; check delX
   if(abs(delX) gt delX_max)then begin
    dPsi = 0.3d * dPsi
    continue
   endif else begin
    psi0 = psi1
    xCoord[iCurve] = xCoord[iCurve] + delX
   endelse

   ; check completion
   if(abs(psi0-psi1_base) lt 1.d-12)then break

  endfor  ; trial

 endfor  ; looping through points on the curve

 ; define file names
 newFile = fPath + 'vanderborght2005_exp' + strtrim(iExp+1, 2) + '_output_timestep.nc'

 ; make the base plot
 plot, indgen(5), xrange=[-0.5,0], yrange=[1,0], xstyle=1, ystyle=1, $
  xtitle='matric head (m)', ytitle = ytitle[iExp], xmargin=[xmar0[iExp],xmar1[iExp]], /nodata
 
 ; open file
 nc_file = ncdf_open(newFile, /nowrite) 

  ; read in the depth
  ivar_id = ncdf_varid(nc_file,'mLayerHeight')
  ncdf_varget, nc_file, ivar_id, depth
  depth = reform(depth)

  ; read in the depth
  ivar_id = ncdf_varid(nc_file,'mLayerMatricHead')
  ncdf_varget, nc_file, ivar_id, matricHead
  matricHead = reform(matricHead)

 ; close netcdf file
 ncdf_close, nc_file

 ; plot data
 oplot, matricHead, depth, color=250, thick=2 

 ; plot the boundary between layers
 plots, [-0.5,0], [0.5,0.5], linestyle=2

 ; define the layers
 xyouts, -0.40, 0.45, upperType[iExp], charsize=3
 xyouts, -0.40, 0.57, lowerType[iExp], charsize=3

 ; overplot the analytical solution
 oplot, psi/100.d, xcoord/100.d, linestyle=2, thick=2
 oplot, [psiLower, psiLower]/100.d, [1,0.5], linestyle=2, thick=2

endfor  ; looping through variables

; save figure
write_png, gpath+gname, tvrd(true=1)

stop
end

; define the parameters
pro get_params, soilType, ksat, alpha, n, m 

case soilType of

 'sand': begin
   ksat  = 1000.d    ; cm d-1
   alpha =   -0.15d  ; cm-1
   n     =    3.d
 end

 'loam': begin
   ksat  =   50.d    ; cm d-1
   alpha =   -0.04d  ; cm-1
   n     =    1.6d
 end
  
 'clay': begin
   ksat  =   10.d    ; cm d-1
   alpha =   -0.01d  ; cm-1
   n     =    1.1d
 end
 
 else: print, 'soilType has an illegal value'

endcase

; compute m
m = 1.d - 1.d/n

end

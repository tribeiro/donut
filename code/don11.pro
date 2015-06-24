; Measuring low-order aberrations from defocused images
; May 3, 2006
;------------------------------------------------
pro readpar, file
; Read/save parameters in a structure

@donut.common

 line = ''
 tmps = ''
 sname= 'donpar'

 res = findfile(file, count=c)

  if (c eq 0) then begin
    print, 'File ',file,' is not found, exiting'
    return
  endif
  openr, unit, file, /get_lun
  while not eof(unit) do begin
    readf, unit, line
    if strpos(line,';') eq -1 then tmps = tmps + ' ' + line $
       else tmps = tmps + ' ' + strmid(line,0,strpos(line,';'))
    if strpos(line,'}') ne -1 then tmps = strcompress(strtrim(tmps,2))	
  endwhile      
  res = execute(sname+' = '+tmps)  ; execute command, create structure
  close, unit 
  free_lun, unit  
  print, 'Parameters are read from ', file

  init
end
;-------------------------------------------------------
pro init
; Pre-compute the parameters and save them in the COMMOM block

@donut.common

   ngrid = donpar.ngrid
   d = donpar.d
   eps = donpar.eps
   lambda = donpar.lambda
   pixel = donpar.pixel
   sflag = 0 ; fast seeing blur, set to 1 for slow calc.


  ; Find reasonable limits on the grid parameters
  asperpix =  206265.*(1e-6*lambda)/D ; maximum fine-pixel size
 
  ccdpix = float(pixel)
  k = floor(alog10(ccdpix/asperpix)/alog10(2.)) +1
  npixperpix = 2^k
  fovpix = 2*ngrid/npixperpix     ; CCD field size
  asperpix = ccdpix/npixperpix
  size = 206266.*(1e-6*lambda)/asperpix
  Rpix = ngrid/size*D

  print, 'Rebinning factor: ', npixperpix
  print, 'Grid pixel: ',asperpix,' arcsec'
  print, 'Grid size: ', size,' m'
  print, 'CCD pixel:  ',ccdpix,' arcsec'
  print, 'CCD field:  ',fovpix*ccdpix,' arcsec'
  print, 'CCD format: ', fovpix

  r = shift(dist(2*ngrid),ngrid,ngrid) ;distance from grid center, pixs
  print,'INIT:',Rpix,Rpix*eps
  inside = where( (r le Rpix) and (r ge Rpix*eps) )
  pupil = fltarr(2*ngrid,2*ngrid)    ; zero array
  pupil(inside) = 1
  n = n_elements(inside)

  x = (findgen(2*ngrid) - ngrid) # replicate(1.,2*ngrid)
  theta = atan(transpose(x),x)
  theta(ngrid,ngrid)=0.
  zgrid = fltarr(n,2)
  zgrid(*,0) = r(inside)/Rpix
  zgrid(*,1) = theta(inside)
  ;print,zgrid(*,0)
end
;-------------------------------------------------------
function getimage, z
; z is the Zernike vector in microns, starting from Z=2 (tip)
; z[0] is seeing in arcseconds
@donut.common

COMMON imagedata, uampl, filter2, seeing

;  t0 = systime(1)
  fact = 2.*!pi/donpar.lambda
  nzer = n_elements(z)
  phase = zgrid(*,0)*0.  ; empty array for phase
  for j=1, nzer-1 do phase += fact*z[j]*zernike_estim(j+1,zgrid)  
print ,'GETIMAGE[zgrid]:',phase[0]
  tmp = fltarr(ngrid*2,ngrid*2)
  uampl = complex(tmp, tmp)
  uampl(inside) = complex(cos(phase), sin(phase))

  seeing = z[0]

;---------  compute the image ----------------------
  imh = abs(shift(fft(shift(uampl,ngrid,ngrid),/inverse),ngrid,ngrid))^2
  if (sflag gt 0) then begin ; exact seing blur
    filter2 = exp(-2.*!pi^2*(seeing/2.35/asperpix/2/ngrid)^2*r^2) ; unbinned image
    imh = abs(fft(shift(fft(imh),ngrid,ngrid)*filter2))
    impix = rebin(imh,fovpix,fovpix) ; rebinning into CCD pixels
  endif else begin 
    rr = shift(dist(fovpix),fovpix/2,fovpix/2) 
    filter2 = exp(-2.*!pi^2*(seeing/2.35/ccdpix/fovpix)^2*rr^2) ; binned image
    impix = rebin(imh,fovpix,fovpix) ; rebinning into CCD pixels
    impix = abs(fft(shift(fft(impix),fovpix/2,fovpix/2)*filter2)) ; Seeing blur
  endelse
;  tvscl, congrid(impix,2*ngrid,2*ngrid)
  return, impix/total(impix)  
end
;-------------------------------------------------------
function newimage, a, jzer
; a is the amplitude change of aberration (micron), Jzer is the Zernike number 
; (1=seeing, 4=focus etc.)

@donut.common
COMMON imagedata, uampl, filter2, seeing

  newampl = uampl
  if (jzer gt 1) then begin  ; Change Zernike coefficient 
    newphase =  2.*!pi/donpar.lambda*a*zernike_estim(jzer,zgrid)  
    newampl(inside) *= complex(cos(newphase), sin(newphase))
    filter = filter2
  endif else begin  ; new seeing
    newseeing = seeing + a
    if (sflag gt 0) then begin
       filter = exp(-2.*!pi^2*(newseeing/2.35/asperpix/2/ngrid)^2*r^2) ; unbinned image
    endif else begin 
     rr = shift(dist(fovpix),fovpix/2,fovpix/2) 
     filter = exp(-2.*!pi^2*(newseeing/2.35/ccdpix/fovpix)^2*rr^2) ; binned image
    endelse
  endelse
;---------  compute the image ----------------------
  imh = abs(shift(fft(shift(newampl,ngrid,ngrid),/inverse),ngrid,ngrid))^2
  if (sflag gt 0) then begin ; exact seing blur
    imh = abs(fft(shift(fft(imh),ngrid,ngrid)*filter))
    impix = rebin(imh,fovpix,fovpix) ; rebinning into CCD pixels
  endif else begin 
    impix = rebin(imh,fovpix,fovpix) ; rebinning into CCD pixels
    impix = abs(fft(shift(fft(impix),fovpix/2,fovpix/2)*filter)) ; Seeing blur
  endelse
;  tvscl, congrid(impix,2*ngrid,2*ngrid)
  return, impix/total(impix)  
end
;-------------------------------------------------------
pro getmom, impix1, zestim
; returns the vector of Zernike aberrations in microns
@donut.common

 n = ngrid/npixperpix
 xx = (findgen(2*n)-n)#replicate(1,2*n)
 yy = transpose(xx)

  thresh = max(impix1)*donpar.thresh
  impix = impix1
  impix (where(impix le thresh)) = 0.

 imh0 = total(impix)
print,'GETMOM[xx*impix]:',total(xx*impix)
 xc = total(xx*impix)/imh0 ;
 yc = total(yy*impix)/imh0 ;
 mxx = total(impix*(xx-xc)^2)/imh0
 myy = total(impix*(yy-yc)^2)/imh0
 mxy = total(impix*(xx-xc)*(yy-yc))/imh0

  scale = npixperpix/(Ngrid/Rpix)

print, 'GETMOM:',scale,xc,yc
  a2 = scale*(xc+0.5)*!pi*0.5
  a3 = scale*(yc+0.5)*!pi*0.5
  print, 'GETMOM:',a2,a3
; stop
  a4 = scale*sqrt((mxx + myy)*0.5)/1.102
  a4 = sqrt( (a4^2 - (0.5/2.35)^2) > 0) ; subtract 0.5arcsec seeing
;  a4 = scale*sqrt((mxx + myy)*0.5)/1.102/1.35
  a5 = scale*mxy*(mxx*myy)^(-0.25)/1.45
  a6 = scale*(mxx - myy)*0.5*(mxx*myy)^(-0.25)/1.45
  zestim = [0.,a2,a3,a4,a5,a6]*donpar.lambda/(2.*!pi) ; estimated Zernike aberrations
  zestim[0] = 0.5 

end
;-------------------------------------------------------
pro displ, image
  tvscl, congrid(image,256,256)
end
;-------------------------------------------------------
pro find, impix, zres, nzer, chi2, model
; Nzer is the highest Zernike number, zres is the  result

@donut.common

  nzer = nzer > 6
  n1 = n_elements(zres)
  nzer1 = n1 > nzer
  z0 = fltarr(nzer1)
  z0[0:n1-1] = zres
       
  impixnorm = impix/total(impix)
  impixmax = max(impix)

  xi = fltarr(nzer)
  for j=1,nzer-1 do xi[j] = donpar.lambda/(2.*!pi)*0.5/(long(sqrt(8L*(j+1L)-6L)-1L)/2L) ; 0.5/n amplitudes
  xi[0] = 0.1
;  xi = xi*0.5 ; Added Dec. 7, 2005
;  xi = xi*2. ; Added Dec. 7, 2005

   indonut = where(impixnorm gt impixmax*donpar.thresh) 
;   indonut = where(impixnorm ge 0.) 
;  sig = sqrt(sigpix(indonut)) ; pixel rms due to noise  

  im = impixnorm(indonut)
  n = n_elements(indonut)
  chi2old = n^2
  chi2norm = total(im^2)

  ncycle = 20
  thresh0 = 0.01  ; initial SVD threshold
  norm = max(impix)
  thresh = thresh0  ; SVD inversion threshold, initial
  print, 'Z  ',indgen(nzer)+1, format='(A4,100I8)'  
;  print, 'um ', z0[0:nzer-1], format='(A4,100F8.3)'
  lambda = 1. ; for L-M method

print,"IMPIX:",size(impix)
print,"INDONUT:",size(indonut)
print,"IM:",size(im)
for k=0,ncycle-1 do begin
  model = getimage(z0)
  print,"MODEL",size(model)
  im0 = model(indonut)
  print,"IM0:",size(im0)
  chi2 = sqrt(total((im0 - im)^2)/chi2norm )
  print, 'Cycle: ',k+1, '  RMS= ',chi2*100, ' percent' 
  print, 'um ', z0[0:nzer-1], format='(A4,100F8.3)'

  thresh = thresh*0.5 ; > 1e-3
;   thresh = chi2*0.1  ; Dec 7, 2005

;  if ((chi2 ge chi2old*0.99) and (k gt 3)) then break
  if (chi2 lt 1e-4) then break


  ; do not degrade aberrations
  if (chi2 le chi2old) then begin
     zres=z0 
     lambda = lambda*0.1
     if ((chi2 ge chi2old*0.99) and (k gt 3)) then break
     chi2old = chi2 
  endif else  begin
     z0 = zres
     thresh = thresh0
     lambda = lambda*10.
     print, 'Divergence... Now LM parameter = ', lambda
  endelse

  
  if (k mod 2 eq 0) then begin
    imat = fltarr(n,nzer)
    print, 'Computing the interaction matrix...'
    print, 'IMAT:',size(imat)
    for j=0,nzer-1 do imat(*,j) =  ((newimage(xi[j],j+1))(indonut) - im0)/xi[j]
;    for j=0,nzer-1 do imat(*,j) =  ((newimage(xi[j],j+1))(indonut) - im0)/xi[j]/sig
    tmat = (transpose(imat) # imat)
;    for j=0,nzer-1 do tmat[j,j] *= (1.+lambda) ; L-M correction
    svd_invert,  tmat, tmp, thresh
    invmat = tmp # transpose(imat)
  endif

  dif = im - im0  
;  dif = (im - im0)/sig 
  dz = invmat # dif  
  z0[0:nzer-1] += 0.7*dz
;  z0[0:nzer-1] += dz
  z0[0] = z0[0] > 0.2
  z0[0] = z0[0] < 1.5

  d1 = min(dif) & d2 = max(dif)
  ; display the image (left: input, right: model)
  tmp = fltarr(512,256)
  tmp(0:255,*) = congrid(impix, 256,256)
  tmp(256:511,*) = congrid(model, 256,256 )
;  tmp(256:511,*) = congrid((impix-model-d1)/(d2-d1)*impixmax, 256,256 )
;  tmp(256:511,*) = congrid((impix-model+0.5*impixmax), 256,256 )
  tvscl, tmp
;  stop

endfor

;  zres = z0[0:nzer-1] ; Do not truncate to permit further fitting

print, 'Fitting done!'

end
;-------------------------------------------------------
pro fit, impix, immod, zres, efoc, chi2
; preliminary and final fitting

@donut.common

 getmom, impix, zres

 nzer = donpar.nzer
 if (donpar.static ne '') then readz, z0, donpar.static else z0 = fltarr(nzer)

 z0[0:5] = zres
 if (efoc lt 0) then z0[3:5] *= -1.
 zres = z0

 find, impix, zres, nzer, chi2, immod

end
;-------------------------------------------------------
pro writepar, file
; Write parameters into a file

@donut.common

  openw, unit, file, /get_lun
  
  n = n_tags(donpar)
  parnames = tag_names(donpar)

  printf, unit, '{donpar, '
  for i=0,n-1 do  begin
  if (i eq n-1) then sep = "}" else sep = ","
     if (datatype(donpar.(i)) eq 'STR') then $
     printf, unit, parnames[i],": '",donpar.(i),"'"+sep else $
     printf, unit, parnames[i],": ",donpar.(i),sep
  endfor
 
   close, unit 
   free_lun, unit  
  print, 'Parameters are saved in ',file
end
;------------------------------------------------
pro savez, z, file
; Save Zernike vector in ASCII file

 openw, unit, file,/get_lun
 fmt='(F8.3,1X,I2)' 
 for i=0,n_elements(z)-1 do printf, unit, z[i], i+1, format=fmt
 close, unit 
 free_lun, unit 
  print, 'Zernike vector is saved in ',file 
end
;------------------------------------------------
pro readz, z, file
; Read Zernike vector from ASCII file

  res = findfile(file, count=c)
  if (c eq 0) then begin
    print, 'File ',file,' is not found, exiting'
    return
  endif

 openr, unit, file,/get_lun
 tmp = fltarr(300)
 i = 0 & x = 0.0
 while not eof(unit) do begin
  readf, unit, x
  tmp[i] = x
  i +=1
 endwhile
 close, unit 
 free_lun, unit 
  if (i eq 0) then z=-1 else z = tmp[0:i-1]
  print, i,' Zernike coefficients are retrieved from ',file 

end
;------------------------------------------------
pro saveres, resfile, z, chi2, imfile

@donut.common

  openw, unit, resfile,/get_lun, /append
  printf, unit, imfile, donpar.xc, donpar.yc, flux, chi2, z,  $
format='(A20,2I6,E10.3,F8.4,100F8.3)'
  close, unit 
  free_lun, unit 
  print, 'Results are saved!'

end
;------------------------------------------------
; image extraction
; ----- function extract, filename, xc, yc, nccd
function extract, img, xc, yc, nccd

@donut.common


  ix1 = xc-nccd >0 & ix2 = xc+nccd < n_elements(img(*,0))-1
  iy1 = yc-nccd >0 & iy2 = yc+nccd < n_elements(img(0,*))-1

  img1 = img[ix1:ix2,iy1:iy2] ; cut out the required part
  img1 = img1 - min(img1)  ; subtract background
  itot = total(img1)

  ; find the center-of-gravity
   nx = n_elements(img1(*,0)) & ny = n_elements(img1(0,*))  
   xx = (findgen(nx)-nx/2)#replicate(1,ny)
   ix = total(img1*xx)/itot + 2

   yy = replicate(1,nx) # (findgen(ny)-ny/2)
   iy = total(img1*yy)/itot +2

  ix = fix(ix)+ nx/2 & iy = fix(iy)+ ny/2

  ix1 = ix-nccd/2 >0 & ix2 = ix1+nccd < nx-1
  iy1 = iy-nccd/2 >0 & iy2 = iy1+nccd < ny-1

  if (ix2-ix1 lt nccd-1) or (iy2-iy1 lt nccd-1) then begin
      print, 'Image is cut on one side!'
      return, -1
  endif

  print, ix1,ix2,iy1,iy2
  impix = img1[ix1:ix2,iy1:iy2]
  i = (sort(impix))[fix(0.1*fovpix^2)]
  print,sort(impix)
  print, 'EXTRACT[i]:',i
  backgr = impix[i] ; 10% quantile of pixel distribution
;  stop
    print, 'EXTRACT[BACKGND]:',backgr
  impix = float(impix) - backgr
  flux = total(impix)
  print, 'Total flux, ADU: ', flux
  impix = impix/flux
  sigpix = (impix > 0)*flux*donpar.eadu + donpar.ron^2  ; variance in each pixel 
;  stop

  print,'IMPIX SIZE:',size(impix)
  return, impix
end
;------------------------------------------------
;------------------------------------------------


;------ end of the program ----------------------------  

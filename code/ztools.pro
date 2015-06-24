; Toolbox of Zernike polynomials
;--------------------------------------------------
FUNCTION cova_zern1,jmax
; Compute the  Noll covariance matrix of Zernike polynomials

c = fltarr(jmax-1,jmax-1)
FOR j=2,jmax DO BEGIN
FOR jp=2,jmax DO BEGIN 
  zern_num,j,n=n,m=m
  zern_num,jp,n=np,m=mp
;k_zz = 7.2e-3 *(!Dpi)^(8./3.)* gamma(14D0/3D0) ; gives 2.242434, c_11 = 0.44815
k_zz = 2.242434

  K = k_zz * (-1)^((n+np-2D0*m)/2D0) *sqrt((n+1)*(np+1))
  IF (m EQ mp) AND ((j*jp/2D0 NE long(j*jp/2D0)) OR ((j/2D0 EQ long(j/2D0)) AND $
                                                   (jp/2D0 EQ long(jp/2D0))) $
                  OR (m EQ 0)) THEN C[j-2,jp-2] = K * gamma((n+np-5D0/3D0)/2D0) / $
          (gamma((n-np+17D0/3D0)/2D0) * gamma((np-n+17D0/3D0)/2D0) * $
           gamma((n+np+23D0/3D0)/2D0)) ELSE C[j-2,jp-2] = 0D0
ENDFOR
ENDFOR
return,c
END
;-------------------------------------------------------------------
PRO Zern_num, num, M = m, N = n, INFO = info
;Description
;  ZERN_NUM computes the azimuthal degree M and the radial order N
;  for the sequential Zernike polynomial NUM
;  INFO  =  result display
    on_error, 2
    IF n_elements(num) NE 1L THEN BEGIN
        message, "SYNTAX: ZERN_NUM, NUM, M=m, N=n"
    END
    IF ((j = long(num))) NE num OR j LT 1L THEN BEGIN
        message, "ERROR: NUM must be an integer greater or equal 1"
    END
    n = long(sqrt(8L*j-7L)-1L)/2L
    IF n MOD 2 THEN BEGIN                  ; odd n
        m = 1L+2L*((j-1L-(n*(n+1L))/2L)/2L)
    END ELSE BEGIN                         ; even n
        m = 2L*((j-(n*(n+1L))/2L)/2L)
    END
    IF keyword_set(info) THEN print, j, n, m
END
;-------------------------------------------------------------------
FUNCTION Fact, n   ;factorial
    on_error, 2
    IF (n_elements(n) NE 1) THEN $
      message, "Syntax: n! = FACT(n) (n must be a scalar integer in [0, 170])"
    IF n LE 18 THEN BEGIN
        f =  1D
        FOR i = 2L, n DO f = f*i
        return, f
    END ELSE BEGIN
        return, exp(total(alog(1+dindgen(n))))
    END
END
;-------------------------------------------------------------------
FUNCTION zernike_estim, mode, grid
; Input: grid = set of points (polar co.) on which the Zernike must be evaluated
; output: vector of Zernike values on the grid
zern_num, mode, M = m, N = n
p = (mode MOD 2)
R=0d0

For J=0,(n-m)/2 Do Begin
   print,"-->",J
   S=J
   R=R+(-1.)^J*Fact(n-J)/(Fact(S)*Fact((n+m)/2-J)$_
                          *Fact((n-m)/2-J))*grid(*,0)^(n-2*J)
Endfor

print, 'ZERNIKE_ESTIM:',mode,(n-m)/2,size(R)
print, 'ZERNIKE_ESTIM[GRID0]:',grid[0,0],grid[1,0],grid[2,0],grid[3,0],grid[4,0],grid[5,0]
print, 'ZERNIKE_ESTIM[GRID1]:',grid[0,1],grid[1,1],grid[2,1],grid[3,1],grid[4,1],grid[5,1]

IF (m EQ 0) Then ZZ=Sqrt(n+1d0)*R

IF (m NE 0) Then Begin
   IF (p EQ 0) Then ZZ=Sqrt(n+1d0)*Sqrt(2d0)*Cos(m*grid(*,1))*R
   IF (p GT 0) Then ZZ=Sqrt(n+1d0)*Sqrt(2d0)*Sin(m*grid(*,1))*R
EndIF
return,zz
END
;----------------------------------------------------------------
;------------------------------------------------------
PRO   svd_invert, matrix, inv_matrix, threshold
; Returns the SVD-inverted matrix 

svdc, matrix, ws, u, v , /double
;svdc, matrix, ws, u, v
;--------- BEWARE! Double-precision is actually needed, but...
;---------- SVDC is bugged in IDL! --------------------------
;print,matrix
ww = max(ws)
n = n_elements(ws)
invw = identity(n)
ncount = 0
FOR i=0, n-1 DO BEGIN 
  IF ws(i) LT ww*threshold THEN BEGIN 
     invw(i,i) = 0.
     ncount = ncount+1
     print ,'SVD_INVERT: Value ',i,'=',ws(i)," rejected (threshold=",threshold,")."
  ENDIF  ELSE invw(i,i)=1./ws(i)
ENDFOR

 print, ncount, ' singular values rejected in inversion'
inv_matrix =  (v ##  invw) ## transpose(u)

END 
;------------------------------------------------------
; Zern_deriv,j
;
;
; This function calculates the x and y derivative coefficients needed to compute 
; the derivative of the jth Zernike polynomial.
; (d/dx)Zj=SUM_j' gammax_j' Z_j'
; (d/dy)Zj=SUM_j' gammay_j' Z_j'
; gammax and gammay is the output vector gamma=[2,j]
;
; Date : 9 December 1999
; Written by Elise Viard, eviard@eso.org


FUNCTION zern_derivx,j

    zern_num, j, M = m, N = n
;   zern_degree, j, n, m
   gam = fltarr(j)
   for j2=1,j do begin
      zern_num, j2, M = m2, N = n2
;      zern_degree, j2, n2, m2
      IF ((m-m2)^2 eq 1) THEN BEGIN
         IF (m NE 0) AND (m2 NE 0) THEN BEGIN
            IF ((j mod(2) EQ 0) AND (j2 mod(2) EQ 0)) OR $
             ((j mod(2) ne 0) AND (j2 mod(2) ne 0)) THEN $
             gam(j2-1) = sqrt((n+1)*(n2+1)) ELSE gam(j2-1) = 0 
         ENDIF ELSE IF ((m EQ 0) AND (j2 mod(2) EQ 0)) THEN $
          gam(j2-1) = sqrt(2.*(n+1)*(n2+1)) ELSE $
          IF ((m2 EQ 0) AND (j mod(2) EQ 0)) THEN gam(j2-1) = sqrt(2.*(n+1)*(n2+1))$
         ELSE  gam(j2-1) = 0 
      ENDIF ELSE gam(j2-1) = 0
   ENDFOR 
   return,gam
END          

;-------------------------------------------------
FUNCTION zern_derivy,j,n,m
    zern_num, j, M = m, N = n
;   zern_degree, j, n, m
   gam = fltarr(j)
   for j2=1,j do begin
      zern_num, j2, M = m2, N = n2
;      zern_degree, j2, n2, m2
      IF ((m-m2)^2 eq 1) THEN BEGIN
         IF (m NE 0) AND (m2 NE 0) THEN BEGIN
            IF ((j mod(2) EQ 0) AND (j2 mod(2) NE 0)) OR $
             ((j2 mod(2) EQ 0) AND (j mod(2) NE 0)) THEN BEGIN
              IF m2 EQ (m+1) AND (j mod(2) NE 0) THEN sig = -1 $
               ELSE IF m2 EQ (m-1) AND (j mod(2) EQ 0) THEN  sig = -1 $
                ELSE sig = 1
               gam(j2-1) =  sig*sqrt((n+1)*(n2+1)) 
            ENDIF ELSE  gam(j2-1) = 0
         ENDIF ELSE IF ((m EQ 0) AND (j2 mod(2) NE 0)) THEN $
          gam(j2-1) = sqrt(2.*(n+1)*(n2+1)) ELSE $
          IF ((m2 EQ 0) AND (j mod(2) NE 0)) THEN gam(j2-1) = sqrt(2.*(n+1)*(n2+1)) $
         ELSE  gam(j2-1) = 0 
       ENDIF ELSE gam(j2-1) = 0
   ENDFOR 
   return,gam
END

;--------------------------------------------
FUNCTION zern_deriv,j
gam = fltarr(2,j)
gam[0,*] = zern_derivx(j) 
gam[1,*] = zern_derivy(j)

return,gam
END 
;-------------------------------------------------------
function getftzer, Jzer
; Compute the Fourier Transform of Zernike mode

; ngrid = 128 ; grid half-size, pixels
; Rpix = 100 ; pupil radius in pixels

 x = (findgen(2*ngrid) - ngrid) # replicate(1.,2*ngrid)
 y = transpose(x)
 theta = atan(y,x)
 theta(ngrid,ngrid)=0.
  
  Zern_num, Jzer, M = m, N = n
  f =  shift(dist(2*ngrid),ngrid,ngrid)/(2*ngrid)*Rpix & f(ngrid,ngrid) = 1e-3

  ftmod = sqrt(n+1D0)*beselj(2*!dpi*f,n+1)/(!dpi*f)   
  IF (m EQ 0) Then  zz = ftmod*complex(0,1D0)^(n/2)
  IF (m NE 0) Then Begin
    iF ((Jzer MOD 2) EQ 0) then fact=sqrt(2.)*Cos(m*theta) $ 
                           else fact=sqrt(2.)*sin(m*theta)
    zz = ftmod*fact*(-1)^((n-m/2))*complex(0,1D0)^m
  EndIF 

;    tvscl, zz
;    tvscl, imaginary(zz), 2*ngrid, 0
; aber = shift(fft(shift(zz,ngrid,ngrid)), ngrid, ngrid) ; get the mode shape
    return, zz

end
;-------------------------------------------------------

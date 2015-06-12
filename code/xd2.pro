; DONUT GUI, Apr 29 2005

;----------------------------------------------------------
;----------------------------------------------------------
pro xdedit_event, ev

@donut.common

common paredit, whichpar, parfile

  widget_control,ev.id, get_uvalue=uvalue 
  if (uvalue ne 'Paredit') then return

 type = tag_names(ev, /structure)
 case type of
'WIDGET_BUTTON': begin
    widget_control,ev.id, get_value=value  
    case value of
   'OK': widget_control, ev.top, /destroy
   'Save param': writepar, donpar.parfile
   else:
   endcase
   end  
;---------------------------
else: begin
    sel = where(ev.id eq whichpar)
    if (sel ne -1) then begin
      widget_control,ev.id, get_value=value  
      value = value[0]  
      donpar.(sel) = value
    endif 
   end
endcase

end
;----------------------------------------------------------
pro xd_event, ev

@donut.common

common paredit, whichpar, parfile

  common xcontrol, impix,immod,zres,z0,chi2,iter,xc,yc,efoc,ws0,ws1,img,label1,label2,label3,button3,button4,button5,button6,button7,button8,zval


  type = tag_names(ev, /structure)
;  print, 'Event type: ', '*'+type+'*'

  case type of
;--------------------
'WIDGET_BUTTON': begin
    widget_control,ev.id, get_value=value  
    case value of
   'Exit': widget_control, ev.top, /destroy
   'Edit parameters': begin
       n = n_tags(donpar)
       whichpar = intarr(n)
       bas0 = widget_base(title='Parameter editor',/column,uvalue='Paredit')
       parname = tag_names(donpar)
       for i=0,n-1 do begin
  case datatype(donpar.(i)) of
'INT':   field = cw_field(bas0,title=parname[i], value = donpar.(i), /all_events, /int,uvalue='Paredit') 
'LON':   field = cw_field(bas0,title=parname[i], value = donpar.(i), /all_events, /lon,uvalue='Paredit') 
'FLO':   field = cw_field(bas0,title=parname[i], value = donpar.(i), /all_events, /flo,uvalue='Paredit') 
'STR':   field = cw_field(bas0,title=parname[i], value = donpar.(i), /all_events, /str,uvalue='Paredit') 

  endcase
        whichpar[i] = field
       endfor
       but_s = widget_button(bas0, value='Save param',uvalue='Paredit')
       but_d = widget_button(bas0, value='OK',uvalue='Paredit')
       widget_control, bas0, /realize, group_leader=base
       xmanager, 'xdedit', bas0, group_leader=base
     end
  'File': begin
        bas0 = widget_base(title='File selector',/column)      
        filesel = cw_filesel(bas0,path=donpar.datadir) 
        widget_control, bas0, /realize
        xmanager, 'xd', bas0
    end
  'Extract': begin
    nx = n_elements(img(*,0)) & ny = n_elements(img(0,*))
    bas0 = widget_base(title='Image display')
    draw1 = widget_draw(bas0,xsize=nx,ysize=ny,x_scroll_size=512,y_scroll_size=512,/scroll,/button_events)
    widget_control, bas0, /realize
    widget_control, draw1, get_value=ws1
    print, 'Image display window is ', ws1
    wset, ws1
     imax = max(img)
     tvscl, (img < imax*0.1)
    xmanager, 'xd', bas0   
   end
  'Fit': begin
     widget_control, /hourglass
     fit, impix, immod, zres, donpar.efoc, chi2
     widget_control, button5, sensitive=1
     widget_control, button6, sensitive=1
     widget_control, button7, sensitive=1
     widget_control, label3, set_value = 'RMS: '+string(100*chi2, format='(F6.2)')+'%'
     widget_control, zval[0], set_value = string(zres[0],format='(F8.3)')
     for i=1,(donpar.nzer-3) < 8 do widget_control, zval[i], set_value = string(zres[i+2],format='(F8.3)')
   end 
  'Fit again': begin
     widget_control, /hourglass
     find, impix, zres, donpar.nzer, chi2, immod
     widget_control, label3, set_value = 'RMS: '+string(100*chi2, format='(F6.2)')+'%'
     widget_control, zval[0], set_value = string(zres[0],format='(F8.3)')
     for i=1,(donpar.nzer-3)<8 do widget_control, zval[i], set_value = string(zres[i+2],format='(F8.3)')   
   end
  'Save result': begin
     saveres, donpar.result, zres, chi2, donpar.imfile
   end
  'Save static': begin
      savez, zres, donpar.static
   end
  'Save image': begin
     tmp = fltarr(512,256)
     tmp(0:255,*) = congrid(impix, 256,256)
     tmp(256:511,*) = congrid(immod, 256,256 )
     set_plot, 'ps'
     device, filename='fit.ps', bits=8
     tvscl, tmp
     device, /close
     set_plot, 'x'      
     print, 'Image is saved in fit.ps'
   end
  ;---
   else: begin
    ;  toggle button of EFOC
    widget_control,ev.id, get_uvalue=uvalue  
    if (uvalue eq 'efoc') then begin
       donpar.efoc *= -1  
       widget_control, label2, set_value = 'EFOC: '+string(donpar.efoc)
     endif
   end
   endcase
  end
  ;--------------------
'FILESEL_EVENT': begin
case ev.done of
0: begin
   donpar.imfile=''
   img = -1   
  end
1: begin
    widget_control,ev.id, get_value=value    
    slash = strpos(value, '/', /reverse_search)
    if (slash gt -1) then begin 
       file = strmid(value,slash+1)
       donpar.datadir = strmid(value,0,slash+1)
    endif
    donpar.imfile = file
    print, 'Selected file: ', donpar.imfile 
    widget_control, ev.top, /destroy
    widget_control, label1, set_value = file
    widget_control, button3, sensitive=1
    img = readfits(donpar.datadir+donpar.imfile)
    if (donpar.xc gt 0) then begin
        impix = extract( img, donpar.xc, donpar.yc, fovpix)
        tvscl, congrid(impix,256,256)
        widget_control, button4, sensitive=1
    endif
     widget_control, button5, sensitive=0
     widget_control, button6, sensitive=0
     widget_control, button7, sensitive=0
     chi2 = 0. & iter= 0
     widget_control, label3, set_value = 'RMS: '+string(100*chi2, format='(F6.2)')+'%'
   end
2:  widget_control, ev.top, /destroy
  endcase

end
;---------------------------
'WIDGET_DRAW': begin
  if (ev.type eq 0) then begin 
    xc = ev.x &  yc = ev.y
    donpar.xc = xc & donpar.yc = yc
    impix = extract(img,xc,yc, fovpix)
    widget_control, ev.top, /destroy
    wset, ws0
    tvscl, congrid(impix,256,256)
    widget_control, button4, sensitive=1 ; fit enable
  endif
end
;---------------------------
else:
endcase
;----------------------
 
end

;------------------------------------ 
pro xd, pfile=pfile

@donut.common

common paredit, whichpar, parfile

       common xcontrol, impix,immod,zres,z0,chi2,iter,xc,yc,efoc,ws0,ws1,img,label1,label2,label3,button3,button4,button5,button6,button7,button8,zval

  if  keyword_set(pfile) then parfile = pfile else parfile = 'donut.par'

  readpar, parfile
  loadct, 3

  iter = 0 & chi2 = 0.
  zres = fltarr(donpar.nzer)

  base = widget_base(title='DONUT GUI',/column)
;  label = widget_label(base,value='--------------- DONUT aberration  retriever ---------------------')
  ; ---- main buttons
  base1 =  widget_base(base, /row)
  button1 = widget_button(base1, value='File')
  button2 = widget_button(base1, value='Edit parameters')
  button3 = widget_button(base1, value='Extract',sensitive=1)
  button4 = widget_button(base1, value='Fit',sensitive=1)
  button5 = widget_button(base1, value='Fit again',sensitive=0)
  button6 = widget_button(base1, value='Save result',sensitive=0)
  button9 = widget_button(base1, value='Exit')

  draw2 = widget_draw(base,xsize=512,ysize=256)


; ---   parameters
  base2 =  widget_base(base, /row)
  label1 = widget_label(base2, value='                         ',/frame)
  label2 =  widget_button(base2, value='EFOC: '+string(donpar.efoc),uvalue='efoc')
  label3 =  widget_label(base2, value='RMS: '+string(100*chi2, format='(F6.2)')+'%',/frame )
  button7 = widget_button(base2, value='Save static',sensitive=0)
  button8 = widget_button(base2, value='Save image')

  base3 =  widget_base(base, /row)
  zl = intarr(9) 
  zl[0] =  widget_label(base3, value='Seeing,"',xsize=60)
  for i=1,8 do zl[i] = widget_label(base3, value='Z'+strtrim(string(i+3)),xsize=60)

  base4 =  widget_base(base, /row)
  zval = intarr(9)
  for i=0,8 do zval[i] = widget_label(base4, value='-',xsize=60)

  widget_control, base, /realize
  widget_control, draw2, get_value=ws0
  wset, ws0

;  --- input data if defined -----------------
  if (donpar.imfile gt '') then begin  
    img = readfits(donpar.datadir+donpar.imfile)
    widget_control, label1, set_value=donpar.imfile
    if (donpar.xc gt 0) then begin
        impix = extract( img, donpar.xc, donpar.yc, fovpix)
        tvscl, congrid(impix,256,256)

    endif
  endif


  xmanager, 'xd', base, /no_block


end

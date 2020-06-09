



subroutine spl2pp (mw, mh, rknot, zknot, copy, breakr, lr, breakz, &
 lz, coef)
!**********************************************************************
!**                                                                  **
!**                                                                  **
!**                                                                  **
!**********************************************************************
! translates to pp representation


use orbit_parameters_mod

implicit real (8)(a - h, o - z)  
save


! parameter (nw = 129, nh = 129, krord = 4, kzord = 4)  

parameter (krord = 4, kzord = 4)  

!      parameter  (lr0=nw-krord+1,lz0=nh-kzord+1)
!      dimension rknot(nw+krord),zknot(nh+kzord)
!      dimension copy(mw,mh),coef(krord,lr0,kzord,lz0)
!      dimension breakr(lr0+1),breakz(lz0+1)
!      dimension work4(krord,nw,nh), work5(nh,krord,lr0)
!     *         ,work6(kzord,kzord,nw,krord)

dimension rknot (mw + krord), zknot (mh + kzord)  
dimension copy (mw, mh), coef (krord, mw - krord+1, kzord, mh - &
 kzord+1)

dimension breakr (mw - krord+2), breakz (mh - kzord+2)  
!     dimension work4(krord,nw,nh), work5(mh,krord,mw-krord+1)
!     *         ,work6(kzord,kzord,nw,krord)
dimension work4 (krord, nw, nh), work5 (nh, krord, nw - krord+1), &
 work6 (kzord, kzord, nw, krord)
equivalence (work4, work6)  
call bspp2d (rknot, copy, mw, krord, mh, work4, breakr, work5, lr)  
ndum = lr * krord  
call bspp2d (zknot, work5, mh, kzord, ndum, work6, breakz, coef, &
 lz)
return  
end subroutine spl2pp

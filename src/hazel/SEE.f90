module SEE
use vars
use maths
use allen
implicit none
contains

!------------------------------------------------------------
! Fill and solve the SEE equations
!------------------------------------------------------------
    subroutine fill_SEE(in_params, in_fixed, component)
    type(variable_parameters) :: in_params
    type(fixed_parameters) :: in_fixed
    integer :: nt, component, error
    real(kind=8) :: nbar, w, gei(0:2)
    integer :: n, lu2, ll2, i
    integer :: ju2, jup2, ku, ku2, qu, qu2, iru, j, jl2, jlp2, ql, ql2, irl, kl, kl2, kr, kr2
    integer :: irp, qp, qp2, kp, kp2, js2, jt2, ir, j2, jp2, k, k2, q, q2, nterm, l2
    real(kind=8) :: x0, x1, x2, x3, x4, x5, x6, x7, sign, sign0
    real(kind=8) :: y0, y1, y2, y3, y4, j10
    real(kind=8) :: flarmor, bcoeff, d, thb, chb
    integer, allocatable :: indx(:)


        if (verbose_mode == 1) then
            print *, 'Number of transitions : ', atom%ntran
        endif
        
        if (.not.allocated(aesto)) allocate(aesto(atom%ntran))
        if (.not.allocated(ntlsto)) allocate(ntlsto(atom%ntran))
        if (.not.allocated(ntusto)) allocate(ntusto(atom%ntran))
        
        if (.not.allocated(SEE_A)) allocate(SEE_A(nrhos,nrhos))
        if (.not.allocated(SEE_mag_A)) allocate(SEE_mag_A(nrhos,nrhos))
        if (.not.allocated(SEE_b)) allocate(SEE_b(nrhos))
        
        SEE_A = 0.d0
        SEE_mag_A = 0.d0
        SEE_b = 0.d0

        do nt = 1, atom%ntran
            
! Use Allen's tables to calculate the anisotropy and the value of nbar
            nbar = nbar_allen(atom%wavelength(nt), in_fixed, in_params, atom%reduction_factor(nt) * in_fixed%nbarExternal(nt))
            w = omega_allen(atom%wavelength(nt), in_fixed, in_params, atom%reduction_factor_omega(nt) * in_fixed%omegaExternal(nt))

! Neglect the influence of anisotropy
            if (in_fixed%use_atomic_pol == -1) then
                w = 0.d0                
            endif
                
            aesto(nt) = atom%ae(nt)
            ntlsto(nt) = atom%nterml(nt)
            ntusto(nt) = atom%ntermu(nt)
            lu2 = lsto2(atom%ntermu(nt))
            ll2 = lsto2(atom%nterml(nt))
            gei(0) = nbar
            gei(1) = nbar * atom%j10(nt)
            gei(2) = nbar * w / sqrt(2.d0)

            
!----------------------------------------------------------------------
!----- Relaxation rates due to spontaneous emission
!----------------------------------------------------------------------
            do i = 1, nrhos
            if (ntab(i) == atom%ntermu(nt)) then
                SEE_A(i,i) = SEE_A(i,i) - atom%ae(nt)
                endif
            enddo
            
!----------------------------------------------------------------------
!----- Transfer rates due to absorption
!----------------------------------------------------------------------
            do i = 1, nrhos
            
            if (ntab(i) == atom%ntermu(nt)) then
                ju2 = j2tab(i)
                jup2 = jp2tab(i)
                ku = ktab(i)
                ku2 = ku*2
                qu = qtab(i)
                qu2 = qu*2
                iru = irtab(i)
                do j = 1, nrhos
                    if (ntab(j) == atom%nterml(nt)) then
                        jl2 = j2tab(j)
                        if (iabs(ju2-jl2) > 2 .or. ju2+jl2 == 0) then
                            else
                            jlp2 = jp2tab(j)
                            if (abs(jup2-jlp2) > 2 .or. jup2+jlp2 == 0) then
                                else
                                ql = qtab(j)
                                ql2 = ql*2
                                if (qu2 == ql2) then
                                    irl = irtab(j)
                                    if (iru ==  irl) then
                                        kl = ktab(j)
                                        kl2 = kl*2
                                        do kr = 0, 2, 1
                                            kr2 = kr*2
                                            if (kr > ku+kl .or. kr < abs(ku-kl)) then
                                                else
                                                x1 = dfloat(lu2+1)
                                                x2 = dsqrt(dfloat(3*(ju2+1)*(jup2+1)*(jl2+1)*(jlp2+1)*(ku2+1)*(kl2+1)*(kr2+1)))
                                                x3 = 1.d0
                                                if (mod(kl2+ql2+jlp2-jl2,4) /= 0) then
                                                        x3 = -1.d0
                                                    endif
                                                x4 = w9js(ju2,jl2,2,jup2,jlp2,2,ku2,kl2,kr2)
                                                x5 = w6js(lu2,ll2,2,jl2,ju2,is2)
                                                x6 = w6js(lu2,ll2,2,jlp2,jup2,is2)
                                                x7 = w3js(ku2,kl2,kr2,-qu2,ql2,0)
                                                x0 = x1*x2*x3*x4*x5*x6*x7*atom%ae(nt)*gei(kr)
                                                SEE_A(i,j) = SEE_A(i,j) + x0
                                                endif
                                            enddo
                                        endif
                                    endif
                                endif
                            endif
                        endif
                    enddo
                
                
!----------------------------------------------------------------------
!----- We repeat the loop for rho^K_Q(j,j') with j>j'
!----- and for rho^K_Q(j,j') with j=j' and Q<0
!----------------------------------------------------------------------
                
                do j = 1, nrhos
                    if (ntab(j) == atom%nterml(nt)) then
                        if (j2tab(j) == jp2tab(j) .and. qtab(j) == 0) then
                            else
                            jl2 = jp2tab(j)
                            if (abs(ju2-jl2) > 2 .or. ju2+jl2 == 0) then
                                else
                                    jlp2 = j2tab(j)
                                    if (abs(jup2-jlp2) > 2 .or. jup2+jlp2 == 0) then
                                    else
                                    ql = -qtab(j)
                                    ql2 = ql*2
                                    if (qu2 == ql2) then
                                        irl = irtab(j)
                                        if (iru ==  irl) then
                                            kl = ktab(j)
                                            kl2 = kl*2
                                            do kr = 0, 2, 1
                                                kr2 = kr*2
                                                if (kr > ku+kl .or. kr < abs(ku-kl)) then
                                                    else
                                                    x1 = dfloat(lu2+1)
                                                    x2 = dsqrt(dfloat(3*(ju2+1)*(jup2+1)*(jl2+1)*(jlp2+1)*(ku2+1)*(kl2+1)*(kr2+1)))
                                                    x3 = 1.d0
                                                    if (mod(kl2+ql2+jlp2-jl2,4) /= 0) then
                                                            x3 = -1.d0
                                                        endif
                                                    x4 = w9js(ju2,jl2,2,jup2,jlp2,2,ku2,kl2,kr2)
                                                    x5 = w6js(lu2,ll2,2,jl2,ju2,is2)
                                                    x6 = w6js(lu2,ll2,2,jlp2,jup2,is2)
                                                    x7 = w3js(ku2,kl2,kr2,-qu2,ql2,0)
                                                    x0 = x1*x2*x3*x4*x5*x6*x7*atom%ae(nt)*gei(kr)
                                                        sign = 1.d0
                                                        if (mod(jl2-jlp2-ql2,4) /= 0) sign = -1.d0
                                                        if (irl == 1) sign0 = sign
                                                        if (irl == 2) sign0 = -sign
                                                    SEE_A(i,j) = SEE_A(i,j) + sign0 * x0
                                                    endif
                                                enddo
                                            endif
                                        endif
                                    endif
                                endif
                            endif
                        endif
                    enddo               
                endif                           
            enddo
            
!----------------------------------------------------------------------
!----- Transfer rates due to spontaneous emission
!----------------------------------------------------------------------         
            
            do i = 1, nrhos
            if (ntab(i) == atom%nterml(nt)) then
                jl2 = j2tab(i)
                jlp2 = jp2tab(i)
                kl = ktab(i)
                kl2 = kl*2
                ql = qtab(i)
                ql2 = ql*2
                irl = irtab(i)
                do j = 1, nrhos
                    if (ntab(j) == atom%ntermu(nt)) then
                    ju2 = j2tab(j)
                        if (iabs(ju2-jl2) > 2 .or. ju2+jl2 == 0) then
                            else
                            jup2 = jp2tab(j)
                            if (iabs(jup2-jlp2) > 2 .or. jup2+jlp2 == 0) then
                                else
                                iru = irtab(j)
                                if (iru == irl) then
                                    ku = ktab(j)
                                    if (kl ==  ku) then
                                        qu = qtab(j)
                                        if (ql == qu) then
                                            x1 = dfloat(lu2+1)
                                            x2 = dsqrt(dfloat((jl2+1)*(jlp2+1)*(ju2+1)*(jup2+1)))
                                            x3 = 1.d0
                                            if (mod(2+kl2+jlp2+jup2,4) /= 0) then
                                                    x3 = -1.d0
                                                endif
                                            x4 = w6js(jl2,jlp2,kl2,jup2,ju2,2)
                                            x5 = w6js(lu2,ll2,2,jl2,ju2,is2)
                                            x6 = w6js(lu2,ll2,2,jlp2,jup2,is2)
                                            x0 = x1*x2*x3*x4*x5*x6*atom%ae(nt)
                                            SEE_A(i,j) = SEE_A(i,j) + x0
                                            endif
                                        endif
                                    endif
                                endif
                            endif
                        endif
                    enddo
                

!----------------------------------------------------------------------
!----- We repeat the loop for rho^K_Q(j,j') with j>j'
!----- and for rho^K_Q(j,j') with j=j' and Q<0
!----------------------------------------------------------------------
                    do j = 1, nrhos
                    if (ntab(j) ==  atom%ntermu(nt)) then
                        if (j2tab(j) == jp2tab(j) .and. qtab(j) == 0) then
                        else
                                ju2 = jp2tab(j)
                            if (iabs(ju2-jl2) > 2 .or. ju2+jl2 == 0) then
                                else
                                jup2 = j2tab(j)
                                if (iabs(jup2-jlp2) > 2 .or. jup2+jlp2 == 0) then
                                    else
                                    iru = irtab(j)
                                    if (iru == irl) then
                                        ku = ktab(j)
                                            if (kl == ku) then
                                            qu = -qtab(j)
                                            qu2 = qu*2
                                            if (ql == qu) then
                                                x1 = dfloat(lu2+1)
                                                x2 = dsqrt(dfloat((jl2+1)*(jlp2+1)*(ju2+1)*(jup2+1)))
                                                x3 = 1.d0
                                                if (mod(2+kl2+jlp2+jup2,4) /= 0) then
                                                        x3 = -1.d0
                                                    endif
                                                x4 = w6js(jl2,jlp2,kl2,jup2,ju2,2)
                                                x5 = w6js(lu2,ll2,2,jl2,ju2,is2)
                                                x6 = w6js(lu2,ll2,2,jlp2,jup2,is2)
                                                x0 = x1*x2*x3*x4*x5*x6*atom%ae(nt)
                                                sign = 1.d0
                                                if (mod(ju2-jup2-qu2,4) /= 0) sign = -1.d0
                                                if (iru == 1) sign0 = sign
                                                if (iru == 2) sign0 = -sign
                                                SEE_A(i,j) = SEE_A(i,j) + sign0*x0
                                                endif
                                            endif
                                        endif
                                    endif
                                endif
                            endif
                        endif
                    enddo
                endif
            enddo
            
!----------------------------------------------------------------------
!----- Relaxation rates due to absorption
!----------------------------------------------------------------------
        do i = 1, nrhos
            if (ntab(i) == atom%nterml(nt)) then
                j2 = j2tab(i)
                jp2 = jp2tab(i)
                k = ktab(i)
                k2 = k*2
                q = qtab(i)
                q2 = q*2
                ir = irtab(i)
                
                    do j = 1, nrhos
                    if (ntab(j) == atom%nterml(nt)) then
                        js2 = j2tab(j)
                        jt2 = jp2tab(j)
                        if (js2 /= j2 .and. jt2 /= jp2) then
                            else
                            kp = ktab(j)
                            kp2 = kp*2
                            qp = qtab(j)
                            if (q == qp) then
                                qp2 = qp*2
                                irp = irtab(j)
                                if (ir == irp) then
                                    do kr = 0, 2, 1
                                        kr2 = kr*2
                                        if (kr <= ll2) then
                                            if (kr > k+kp .or. kr < iabs(k-kp)) then
                                                else
                                                x1 = dfloat(lu2+1)
                                                x2 = dsqrt(dfloat(3*(k2+1)*(kp2+1)*(kr2+1)))
                                                x3 = 1.d0
                                                if (mod(2+lu2-is2+j2+qp2,4) /= 0) then
                                                        x3 = -1.d0
                                                    endif
                                                x4 = w6js(ll2,ll2,kr2,2,2,lu2)
                                                x5 = w3js(k2,kp2,kr2,q2,-qp2,0)
                                                x6 = x1*x2*x3*x4*x5
                                                y0 = 0.d0
                                                
!--------------------------------------------------------------------
!----- This is calculated only if j=js
!--------------------------------------------------------------------
                                                if (j2 == js2) then
                                                    if (kr2 > jp2+jt2 .or. kr2 < iabs(jp2-jt2)) then
                                                        else
                                                        y1 = dsqrt(dfloat((jp2+1)*(jt2+1)))
                                                        y2 = w6js(ll2,ll2,kr2,jt2,jp2,is2)
                                                        y3 = w6js(k2,kp2,kr2,jt2,jp2,j2)
                                                        y0 = y1*y2*y3
                                                        endif
                                                    endif
                                                                                        
!---------------------------------------------------------------------
!----- This is calculated only if j=jt
!---------------------------------------------------------------------
                                                if (jp2 == jt2) then
                                                    if (kr2 > j2+js2 .or. kr2 < iabs(j2-js2)) then
                                                        else
                                                        y1 = dsqrt(dfloat((j2+1)*(js2+1)))
                                                        y2 = 1.d0
                                                        if (mod(js2-jp2+k2+kp2+kr2,4) /= 0) then
                                                                y2 = -1.d0
                                                            endif
                                                        y3 = w6js(ll2,ll2,kr2,js2,j2,is2)
                                                        y4 = w6js(k2,kp2,kr2,js2,j2,jp2)
                                                        y0 = y0+y1*y2*y3*y4
                                                        endif
                                                    endif

                                                x0 = 0.5d0*x6*y0*atom%ae(nt)*gei(kr)
                                                SEE_A(i,j) = SEE_A(i,j) - x0
                                                
                                                endif
                                            endif
                                        enddo
                                    endif
                                endif
                            endif
                        endif
                    enddo
                    
!----------------------------------------------------------------------
!----- We repeat the loop for rho^K_Q(j,j') with j>j'
!----- and for rho^K_Q(j,j') with j=j' and Q<0
!----------------------------------------------------------------------                 

                    do j = 1, nrhos
                    if (ntab(j) == atom%nterml(nt)) then
                        if (j2tab(j) == jp2tab(j) .and. qtab(j) == 0) then
                            else
                                
                                js2 = jp2tab(j)
                            jt2 = j2tab(j)
                            if (js2 /= j2 .and. jt2 /= jp2) then
                                else
                                kp = ktab(j)
                                kp2 = kp*2
                                qp = -qtab(j)
                                if (q == qp) then
                                    qp2 = qp*2
                                    irp = irtab(j)
                                    if (ir == irp) then
                                        do kr = 0, 2, 1
                                            kr2 = kr*2
                                            if (kr <= ll2) then
                                                if (kr > k+kp .or. kr < iabs(k-kp)) then
                                                    else
                                                    x1 = dfloat(lu2+1)
                                                    x2 = dsqrt(dfloat(3*(k2+1)*(kp2+1)*(kr2+1)))
                                                    x3 = 1.d0
                                                    if (mod(2+lu2-is2+j2+qp2,4) /= 0) then
                                                            x3 = -1.d0
                                                        endif
                                                    x4 = w6js(ll2,ll2,kr2,2,2,lu2)
                                                    x5 = w3js(k2,kp2,kr2,q2,-qp2,0)
                                                    x6 = x1*x2*x3*x4*x5
                                                    y0 = 0.d0
                                                        
!--------------------------------------------------------------------
!----- This is calculated only if j=js
!--------------------------------------------------------------------
                                                    if (j2 == js2) then
                                                        if (kr2 > jp2+jt2 .or. kr2 < iabs(jp2-jt2)) then
                                                            else
                                                            y1 = dsqrt(dfloat((jp2+1)*(jt2+1)))
                                                            y2 = w6js(ll2,ll2,kr2,jt2,jp2,is2)
                                                            y3 = w6js(k2,kp2,kr2,jt2,jp2,j2)
                                                            y0 = y1*y2*y3
                                                            endif
                                                        endif
                                                                                        
!---------------------------------------------------------------------
!----- This is calculated only if j=jt
!---------------------------------------------------------------------
                                                    if (jp2 == jt2) then
                                                        if (kr2 > j2+js2 .or. kr2 < iabs(j2-js2)) then
                                                            else
                                                            y1 = dsqrt(dfloat((j2+1)*(js2+1)))
                                                            y2 = 1.d0
                                                            if (mod(js2-jp2+k2+kp2+kr2,4) /= 0) then
                                                                    y2 = -1.d0
                                                                endif
                                                            y3 = w6js(ll2,ll2,kr2,js2,j2,is2)
                                                            y4 = w6js(k2,kp2,kr2,js2,j2,jp2)
                                                            y0 = y0+y1*y2*y3*y4
                                                            endif
                                                        endif

                                                    x0 = 0.5d0*x6*y0*atom%ae(nt)*gei(kr)
                                                        sign = 1.d0
                                                        if (mod(js2-jt2-qp2,4) /= 0) then
                                                            sign = -1.d0
                                                        endif
                                                        if (irp == 1) sign0 = sign
                                                        if (irp == 2) sign0 = -sign
                                                    SEE_A(i,j) = SEE_A(i,j) - sign0*x0
                                                    endif
                                                endif
                                            enddo
                                        endif
                                    endif
                                endif
                            endif
                        endif
                    enddo
                    
                endif
            enddo
            
!----------------------------------------------------------------------
!----- If stimulated emission is included...
!----------------------------------------------------------------------
            if (isti == 1) then
            
!----------------------------------------------------------------------
!----- Transfer rates due to stimulated emission
!----------------------------------------------------------------------

                do i = 1, nrhos
                
                    if (ntab(i) == atom%nterml(nt)) then
                    jl2 = j2tab(i)
                    jlp2 = jp2tab(i)
                    kl = ktab(i)
                      kl2 = kl*2
                      ql = qtab(i)
                      ql2 = ql*2
                      irl = irtab(i)
                      do j = 1, nrhos
                          if (ntab(j) == atom%ntermu(nt)) then
                              ju2 = j2tab(j)
                              if (iabs(jl2-ju2) > 2 .or. jl2+ju2 == 0) then
                                else
                                  jup2 = jp2tab(j)
                                  if (iabs(jlp2-jup2) > 2 .or. jlp2+jup2 == 0) then
                                    else
                                      qu = qtab(j)
                                      qu2 = qu*2
                                      if (ql2 == qu2) then
                                          iru = irtab(j)
                                          if (irl == iru) then
                                              ku = ktab(j)
                                              ku2 = ku*2
                                              do kr = 0, 2, 1
                                                  kr2 = kr*2
                                                  if (kr > kl+ku .or. kr < iabs(kl-ku)) then
                                                    else
                                                      x1 = dfloat(lu2+1)
                                                      x2 = dsqrt(dfloat(3*(jl2+1)*(jlp2+1)*(ju2+1)*(jup2+1)*(kl2+1)*(ku2+1)*(kr2+1)))
                                                      x3 = 1.d0
                                                      if (mod(kr2+ku2+qu2+jup2-ju2,4) /= 0) then
                                                            x3 = -1.d0
                                                        endif
                                                      x4 = w9js(jl2,ju2,2,jlp2,jup2,2,kl2,ku2,kr2)
                                                      x5 = w6js(lu2,ll2,2,jl2,ju2,is2)
                                                      x6 = w6js(lu2,ll2,2,jlp2,jup2,is2)
                                                      x7 = w3js(kl2,ku2,kr2,-ql2,qu2,0)
                                                      x0 = x1*x2*x3*x4*x5*x6*x7*atom%ae(nt)*gei(kr)
                                                      SEE_A(i,j) = SEE_A(i,j) + x0
                                                    endif
                                                enddo
                                            endif
                                        endif
                                    endif
                                endif
                            endif
                        enddo
                        
!----------------------------------------------------------------------
!----- We repeat the loop for rho^K_Q(j,j') with j>j'
!----- and for rho^K_Q(j,j') with j=j' and Q<0
!----------------------------------------------------------------------                         
                    do j = 1, nrhos
                        if (ntab(j) == atom%ntermu(nt)) then
                             if (j2tab(j) ==  jp2tab(j) .and. qtab(j) == 0) then
                                else
                                 ju2 = jp2tab(j)
                                 if (iabs(jl2-ju2)> 2 .or. jl2+ju2 == 0) then
                                    else
                                     jup2 = j2tab(j)
                                     if (iabs(jlp2-jup2) > 2 .or. jlp2+jup2 == 0) then
                                        else
                                         qu = -qtab(j)
                                         qu2 = qu*2
                                         if (ql2 == qu2) then
                                             iru = irtab(j)
                                             if (irl == iru) then
                                                 ku = ktab(j)
                                                 ku2 = ku*2
                                                 do kr = 0, 2, 1
                                                     kr2 = kr*2
                                                     if (kr > kl+ku .or. kr < iabs(kl-ku)) then
                                                        else
                                                         x1 = dfloat(lu2+1)
                                                         x2 = dsqrt(dfloat(3*(jl2+1)*(jlp2+1)*(ju2+1)*(jup2+1)*(kl2+1)*(ku2+1)*(kr2+1)))
                                                         x3 = 1.d0
                                                         if (mod(kr2+ku2+qu2+jup2-ju2,4) /= 0) then
                                                                x3 = -1.d0
                                                            endif
                                                          x4 = w9js(jl2,ju2,2,jlp2,jup2,2,kl2,ku2,kr2)
                                                         x5 = w6js(lu2,ll2,2,jl2,ju2,is2)
                                                         x6 = w6js(lu2,ll2,2,jlp2,jup2,is2)
                                                         x7 = w3js(kl2,ku2,kr2,-ql2,qu2,0)
                                                         x0 = x1*x2*x3*x4*x5*x6*x7*atom%ae(nt)*gei(kr)
                                                         sign = 1.d0
                                                         if (mod(ju2-jup2-qu2,4) /= 0) sign = -1.d0
                                                         if (iru == 1) sign0 = sign
                                                         if (iru == 2) sign0 = -sign
                                                         SEE_A(i,j) = SEE_A(i,j) + sign0*x0
                                                        endif
                                                    enddo
                                                endif
                                            endif
                                        endif
                                    endif
                                endif
                            endif
                        enddo
                    
                    endif
                enddo
                
!----------------------------------------------------------------------
!----- Relaxation rates due to stimulated emission
!----------------------------------------------------------------------
            do i = 1, nrhos
                if (ntab(i) == atom%ntermu(nt)) then
                    j2 = j2tab(i)
                    jp2 = jp2tab(i)
                    k = ktab(i)
                    k2 = k*2
                    q = qtab(i)
                    q2 = q*2
                    ir = irtab(i)

                        do j = 1, nrhos
                        if (ntab(j) == atom%ntermu(nt)) then
                            js2 = j2tab(j)
                            jt2 = jp2tab(j)
                            if (js2 /= j2 .and. jt2 /= jp2) then
                                else
                                kp = ktab(j)
                                kp2 = kp*2
                                qp = qtab(j)
                                if (q == qp) then
                                    qp2 = qp*2
                                    irp = irtab(j)
                                    if (ir == irp) then
                                        do kr = 0, 2, 1
                                            kr2 = kr*2
                                            if (kr <= lu2) then
                                                if (kr > k+kp .or. kr < iabs(k-kp)) then
                                                    else
                                                    x1 = dfloat(lu2+1)
                                                    x2 = dsqrt(dfloat(3*(k2+1)*(kp2+1)*(kr2+1)))
                                                    x3 = 1.d0
                                                    if (mod(2+ll2-is2+j2+kr2+qp2,4) /= 0) then
                                                            x3 = -1.d0
                                                        endif
                                                    x4 = w6js(lu2,lu2,kr2,2,2,ll2)
                                                    x5 = w3js(k2,kp2,kr2,q2,-qp2,0)
                                                    x6 = x1*x2*x3*x4*x5
                                                    y0 = 0.d0

!--------------------------------------------------------------------
!----- This is calculated only if j=js
!--------------------------------------------------------------------
                                                    if (j2 == js2) then
                                                        if (kr2 > jp2+jt2 .or. kr2 < iabs(jp2-jt2)) then
                                                            else
                                                            y1 = dsqrt(dfloat((jp2+1)*(jt2+1)))
                                                            y2 = w6js(lu2,lu2,kr2,jt2,jp2,is2)
                                                            y3 = w6js(k2,kp2,kr2,jt2,jp2,j2)
                                                            y0 = y1*y2*y3
                                                            endif
                                                        endif

!---------------------------------------------------------------------
!----- This is calculated only if j=jt
!---------------------------------------------------------------------
                                                    if (jp2 == jt2) then
                                                        if (kr2 > j2+js2 .or. kr2 < iabs(j2-js2)) then
                                                            else
                                                            y1 = dsqrt(dfloat((j2+1)*(js2+1)))
                                                            y2 = 1.d0
                                                            if (mod(js2-jp2+k2+kp2+kr2,4) /= 0) then
                                                                    y2 = -1.d0
                                                                endif
                                                            y3 = w6js(lu2,lu2,kr2,js2,j2,is2)
                                                            y4 = w6js(k2,kp2,kr2,js2,j2,jp2)
                                                            y0 = y0+y1*y2*y3*y4
                                                            endif
                                                        endif

                                                    x0 = 0.5d0*x6*y0*atom%ae(nt)*gei(kr)
                                                    SEE_A(i,j) = SEE_A(i,j) - x0

                                                    endif
                                                endif
                                            enddo
                                        endif
                                    endif
                                endif
                            endif
                        enddo

!----------------------------------------------------------------------
!----- We repeat the loop for rho^K_Q(j,j') with j>j'
!----- and for rho^K_Q(j,j') with j=j' and Q<0
!----------------------------------------------------------------------                 

                        do j = 1, nrhos
                        if (ntab(j) == atom%ntermu(nt)) then
                            if (j2tab(j) == jp2tab(j) .and. qtab(j) == 0) then
                                else

                                    js2 = jp2tab(j)
                                jt2 = j2tab(j)
                                if (js2 /= j2 .and. jt2 /= jp2) then
                                    else
                                    kp = ktab(j)
                                    kp2 = kp*2
                                    qp = -qtab(j)
                                    if (q == qp) then
                                        qp2 = qp*2
                                        irp = irtab(j)
                                        if (ir == irp) then
                                            do kr = 0, 2, 1
                                                kr2 = kr*2
                                                if (kr <= lu2) then
                                                    if (kr > k+kp .or. kr < iabs(k-kp)) then
                                                        else
                                                        x1 = dfloat(lu2+1)
                                                        x2 = dsqrt(dfloat(3*(k2+1)*(kp2+1)*(kr2+1)))
                                                        x3 = 1.d0
                                                        if (mod(2+ll2-is2+j2+kr2+qp2,4) /= 0) then
                                                                x3 = -1.d0
                                                            endif
                                                        x4 = w6js(lu2,lu2,kr2,2,2,ll2)
                                                        x5 = w3js(k2,kp2,kr2,q2,-qp2,0)
                                                        x6 = x1*x2*x3*x4*x5
                                                        y0 = 0.d0

!--------------------------------------------------------------------
!----- This is calculated only if j=js
!--------------------------------------------------------------------
                                                        if (j2 == js2) then
                                                            if (kr2 > jp2+jt2 .or. kr2 < iabs(jp2-jt2)) then
                                                                else
                                                                y1 = dsqrt(dfloat((jp2+1)*(jt2+1)))
                                                                y2 = w6js(lu2,lu2,kr2,jt2,jp2,is2)
                                                                y3 = w6js(k2,kp2,kr2,jt2,jp2,j2)
                                                                y0 = y1*y2*y3
                                                                endif
                                                            endif

!---------------------------------------------------------------------
!----- This is calculated only if j=jt
!---------------------------------------------------------------------
                                                        if (jp2 == jt2) then
                                                            if (kr2 > j2+js2 .or. kr2 < iabs(j2-js2)) then
                                                                else
                                                                y1 = dsqrt(dfloat((j2+1)*(js2+1)))
                                                                y2 = 1.d0
                                                                if (mod(js2-jp2+k2+kp2+kr2,4) /= 0) then
                                                                        y2 = -1.d0
                                                                    endif
                                                                y3 = w6js(lu2,lu2,kr2,js2,j2,is2)
                                                                y4 = w6js(k2,kp2,kr2,js2,j2,jp2)
                                                                y0 = y0+y1*y2*y3*y4
                                                                endif
                                                            endif

                                                        x0 = 0.5d0*x6*y0*atom%ae(nt)*gei(kr)
                                                            sign = 1.d0
                                                            if (mod(js2-jt2-qp2,4) /= 0) then
                                                                sign = -1.d0
                                                            endif
                                                            if (irp == 1) sign0 = sign
                                                            if (irp == 2) sign0 = -sign
                                                        SEE_A(i,j) = SEE_A(i,j) - sign0*x0
                                                        endif
                                                    endif
                                                enddo
                                            endif
                                        endif
                                    endif
                                endif
                            endif
                        enddo

                    endif
                enddo               
 
            
            endif   ! Stimulated emission if
                        
        enddo   ! Loop over transitions
        
        close(12)
        
!----------------------------------------------------------------------
!----- We add the imaginary term produced by the energy difference between 
!----- the two levels of the coherence
!----------------------------------------------------------------------
      do i = 1, nrhos
        if (j2tab(i) /= jp2tab(i)) then
            if (irtab(i) == 1) then
                    SEE_A(i,i+1) = SEE_A(i,i+1) + 2.d0*PI*rnutab(i)
                endif
            if (irtab(i) == 2) then
                    SEE_A(i,i-1) = SEE_A(i,i-1) - 2.d0*PI*rnutab(i)
                endif
            endif
        enddo

!-----------------------------------------------------------------------
!----- If idep=1, we add the term due to depolarizing collisions (only for the fundamental level)
!-----------------------------------------------------------------------
        if (idep == 1) then
            do i = 1, nrhos
                if (ntab(i) == 1 .and. ktab(i) /= 0) then
                    SEE_A(i,i) = SEE_A(i,i) - in_params%delta_collision
                endif
            enddo           
        endif
        
!-----------------------------------------------------------------------
!----- If imag=1, we add the term produced by the magnetic field. We first build
!----- the magnetic matrix SEE_mag_A and then add this contribution to the total matrix
!-----------------------------------------------------------------------        
      if (imag == 1) then
        
!-----------------------------------------------------------------------
!----- We get the magnetic parameters
!-----------------------------------------------------------------------

! Which component we are computing
            if (component == 1) then
            thb = in_params%thetabd * PI / 180.d0       
            chb = in_params%chibd * PI / 180.d0
        else
            thb = in_params%thetabd2 * PI / 180.d0          
            chb = in_params%chibd2 * PI / 180.d0
        endif
        
! Random azimuth solution. Since the density matrix in the magnetic field
! reference frame is independent on azimuth, in the random azimuth solution
! I think one can do the solution for any azimuth. The reason is that, even 
! if the SEE are solved for the density matrix in the vertical reference frame,
! they are later transformed to the magnetic field reference frame before 
! calculating the emission/absorption coefficients
        if (in_params%chibd == 999.d0) chb = 40.d0 * PI / 180.d0
        
        call rotmatk (1,0.d0,-thb,-chb,dr,di)

! Which component we are computing
        if (component == 1) then
            flarmor = in_params%bgauss * (PE / (4.d0*PI*PME*PC))
        else
            flarmor = in_params%bgauss2 * (PE / (4.d0*PI*PME*PC))
        endif
        bcoeff = 2.d0 * PI * flarmor

!-----------------------------------------------------------------------
!----- We build the SEE_mag_A matrix
!-----------------------------------------------------------------------
          do i = 1, nrhos
              nterm = ntab(i)
            l2 = lsto2(nterm)
              j2 = j2tab(i)
              jp2 = jp2tab(i)
              k = ktab(i)
              k2 = k*2
              q = qtab(i)
              q2 = q*2
              ir = irtab(i)
              do j = 1, nrhos
                  if (ntab(j) == nterm) then
                      js2 = j2tab(j)
                      jt2 = jp2tab(j)
                      kp = ktab(j)
                      if (iabs(k-kp) > 1 .or. k+kp == 0) then
                        else
                          kp2 = kp*2
                          qp = qtab(j)
                          if (iabs(q-qp) <= 1) then
                              qp2 = qp*2
                              irp = irtab(j)
                              x1 = 1.d0
                              if (mod(j2+jp2-qp2,4) /= 0) x1 = -1.d0
                              x2 = w3js(k2,kp2,2,-q2,qp2,q2-qp2)
                              call bigpar(l2,is2,k2,kp2,j2,jp2,js2,jt2,x3)
                              if (ir == 1 .and. irp == 1) x4 = di(0,qp-q)
                              if (ir == 1 .and. irp == 2) x4 = dr(0,qp-q)
                              if (ir == 2 .and. irp == 1) x4 = -dr(0,qp-q)
                              if (ir == 2 .and. irp == 2) x4 = di(0,qp-q)
                              x5 = dsqrt(dfloat((k2+1)*(kp2+1)))
                              SEE_mag_A(i,j) = SEE_mag_A(i,j) + bcoeff*x1*x2*x3*x4*x5
                            endif
                        endif
                    endif
                enddo

!----------------------------------------------------------------------
!----- We repeat the loop for rho^K_Q(j,j') with j>j'
!----- and for rho^K_Q(j,j') with j=j' and Q<0
!---------------------------------------------------------------------- 
            do j = 1, nrhos
                if (ntab(j) == nterm) then
                    if (j2tab(j) == jp2tab(j) .and. qtab(j) == 0) then
                        else
                        js2 = jp2tab(j)
                         jt2 = j2tab(j)
                         kp = ktab(j)
                         if (iabs(k-kp) > 1 .or. k+kp == 0) then
                            else
                             kp2 = kp*2
                             qp = -qtab(j)
                             if (iabs(q-qp) <= 1) then
                                 qp2 = qp*2
                                 irp = irtab(j)
                                 x1 = 1.d0
                                 if (mod(j2+jp2-qp2,4) /= 0) x1 = -1.d0
                                  x2 = w3js(k2,kp2,2,-q2,qp2,q2-qp2)
                                  call bigpar(l2,is2,k2,kp2,j2,jp2,js2,jt2,x3)
                                  if (ir == 1.and.irp == 1) x4 = di(0,qp-q)
                                  if (ir == 1.and.irp == 2) x4 = dr(0,qp-q)
                                if (ir == 2.and.irp == 1) x4 = -dr(0,qp-q)
                                 if (ir == 2.and.irp == 2) x4 = di(0,qp-q)
                                 x5 = dsqrt(dfloat((k2+1)*(kp2+1)))
                                 sign = 1.d0
                                 if (mod(js2-jt2-qp2,4) /= 0) sign = -1.d0
                                 if (irp == 1) sign0 = sign
                                 if (irp == 2) sign0 = -sign
                                 SEE_mag_A(i,j) = SEE_mag_A(i,j) + bcoeff*x1*x2*x3*x4*x5*sign0
                                endif
                            endif
                        endif
                    endif
                enddo
            enddo
            
        endif  ! Include magnetic field
        
! Now we add the magnetic contribution to the system matrix
        SEE_A = SEE_A + SEE_mag_A
        
! We substitute the first row with the trace equation (multiplied by 10^7)
        SEE_A(1,:) = 0.d0
        do j = 1, nrhos
            if (ktab(j) == 0) then
                SEE_A(1,j) = 1.d7 * dsqrt(dfloat(j2tab(j)+1))
            endif
        enddo

! Now the right hand side of the equation
        SEE_b(1) = 1.d7

! And we now solve the system

        select case(linear_solver)
        case(0)
            if (verbose_mode == 1) then
                print *, 'Solving the system of equations using a LU decomposition'
                print *, 'Size of the system : ', nrhos,'x',nrhos
            endif
            allocate(indx(nrhos))

            call ludcmp(SEE_A,indx,d)
            
            if (error_code == 1) then
                deallocate(indx)
                return
            else                
                call lubksb(SEE_A,indx,SEE_b)
            endif
            
            deallocate(indx)
        case(1)         
            if (verbose_mode == 1) then
                print *, 'Solving the system of equations using an Incomplete LU BiConjugate Gradient'
                print *, 'Size of the system : ', nrhos,'x',nrhos
            endif
            print *, 'Not implemented'
!           call slapsolver(SEE_A,SEE_b)
        end select

    
! Write the solution
!       open(unit=7,file='tanti.res',status='unknown')
!       do i = 1, nrhos 
!           write(7,FMT='(7(1x,i5),1x,e15.7)') i,ntab(i),j2tab(i),jp2tab(i),ktab(i),qtab(i),&
!               irtab(i),SEE_b(i)
!       enddo
!       close(7)
        
!       deallocate(aesto)
!       deallocate(ntlsto)
!       deallocate(ntusto)
        
!       deallocate(SEE_A)
!       deallocate(SEE_mag_A)
!       deallocate(SEE_b)
        
    end subroutine fill_SEE
end module SEE

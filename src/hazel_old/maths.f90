module maths
use vars
! use Fast_Fourier
use singleton
!use nrutil, ONLY : assert_eq, imaxloc, outerprod, swap
implicit none
contains

!----------------------------------------------------------------
! This function calculates the 3-j symbol
! J_i and M_i have to be twice the actual value of J and M
!----------------------------------------------------------------
   function w3js(j1,j2,j3,m1,m2,m3)
   integer :: m1, m2, m3, j1, j2, j3
    integer :: ia, ib, ic, id, ie, im, ig, ih, z, zmin, zmax, jsum
    real(kind=8) :: w3js, cc, denom, cc1, cc2

!       w3js = w3js_regge(j1/2,j2/2,j3/2,m1/2,m2/2,m3/2)
!       return
      w3js = 0.d0
      if (m1+m2+m3 /= 0) return
      ia = j1 + j2
      if (j3 > ia) return
      ib = j1 - j2
      if (j3 < abs(ib)) return
        if (abs(m1) > j1) return
        if (abs(m2) > j2) return
        if (abs(m3) > j3) return
        
      jsum = j3 + ia
      ic = j1 - m1
      id = j2 - m2
      
        ie = j3 - j2 + m1
        im = j3 - j1 - m2
        zmin = max0(0,-ie,-im)
        ig = ia - j3
        ih = j2 + m2
        zmax = min0(ig,ih,ic)
        cc = 0.d0

        do z = zmin, zmax, 2
            denom = fact(z/2)*fact((ig-z)/2)*fact((ic-z)/2)*fact((ih-z)/2)*&
                fact((ie+z)/2)*fact((im+z)/2)
        if (mod(z,4) /= 0) denom = -denom
            cc = cc + 1.d0/denom
        enddo

        cc1 = fact(ig/2)*fact((j3+ib)/2)*fact((j3-ib)/2)/fact((jsum+2)/2)
      cc2 = fact((j1+m1)/2)*fact(ic/2)*fact(ih/2)*fact(id/2)*fact((j3-m3)/2)*fact((j3+m3)/2)
      cc = cc * sqrt(1.d0*cc1*cc2)
        if (mod(ib-m3,4) /= 0) cc = -cc
        w3js = cc
        if (abs(w3js) < 1.d-8) w3js = 0.d0      
1000        return
   end function w3js
    
    
!----------------------------------------------------------------
! This function calculates the 6-j symbol
! J_i and M_i have to be twice the actual value of J and M
!----------------------------------------------------------------
    function w6js(j1,j2,j3,l1,l2,l3)
    integer :: j1,j2,j3,l1,l2,l3
    integer :: ia, ib, ic, id, ie, iif, ig, ih, sum1, sum2, sum3, sum4
    integer :: w, wmin, wmax, ii, ij, ik
   real(kind=8) :: w6js, omega, theta1, theta2, theta3, theta4, theta, denom

      w6js = 0.d0
        ia = j1 + j2
        if (ia < j3) return
        ib = j1 - j2
        if (abs(ib) > j3) return
        ic = j1 + l2
        if (ic < l3) return
        id = j1 - l2
        if (abs(id) > l3) return
        ie = l1 + j2
        if (ie < l3) return
        iif = l1 - j2
        if (abs(iif) > l3) return
        ig = l1 + l2
        if (ig < j3) return
        ih = l1 - l2
        if (abs(ih) > j3) return
      sum1=ia + j3
      sum2=ic + l3
      sum3=ie + l3
      sum4=ig + j3
        wmin = max0(sum1, sum2, sum3, sum4)
        ii = ia + ig
        ij = j2 + j3 + l2 + l3
        ik = j3 + j1 + l3 + l1
        wmax = min0(ii,ij,ik)
        omega = 0.d0
        do w = wmin, wmax, 2
          denom = fact((w-sum1)/2)*fact((w-sum2)/2)*fact((w-sum3)/2)&
                *fact((w-sum4)/2)*fact((ii-w)/2)*fact((ij-w)/2)&
                *fact((ik-w)/2)
            if (mod(w,4) /= 0) denom = -denom
            omega = omega + fact(w/2+1) / denom
        enddo     
      theta1 = fact((ia-j3)/2)*fact((j3+ib)/2)*fact((j3-ib)/2)/fact(sum1/2+1)
      theta2 = fact((ic-l3)/2)*fact((l3+id)/2)*fact((l3-id)/2)/fact(sum2/2+1)
      theta3 = fact((ie-l3)/2)*fact((l3+iif)/2)*fact((l3-iiF)/2)/fact(sum3/2+1)
      theta4 = fact((ig-j3)/2)*fact((j3+ih)/2)*fact((j3-ih)/2)/fact(sum4/2+1)
      theta = theta1 * theta2 * theta3 * theta4
        w6js = omega * sqrt(theta)
      if (abs(w6js) < 1.d-8) w6js = 0.d0
1000        return
   end function w6js

!----------------------------------------------------------------
! This function calculates the 9-j symbol
! J_i and M_i have to be twice the actual value of J and M
!----------------------------------------------------------------
   function w9js(j1,j2,j3,j4,j5,j6,j7,j8,j9)
    integer :: j1,j2,j3,j4,j5,j6,j7,j8,j9
    integer :: i, kmin, kmax, k
    real(kind=8) :: x, s, x1, x2, x3, w9js

      kmin = abs(j1-j9)
        kmax = j1 + j9
        i = abs(j4-j8)
        if (i > kmin) kmin = i
        i = j4 + j8
        if (i < kmax) kmax = i
        i = abs(j2-j6)
        if (i > kmin) kmin = i
        i = j2 + j6
        if (i < kmax) kmax = i
        x = 0.d0
        do k = kmin, kmax, 2
            s = 1.d0
            if (mod(k,2) /= 0) s = -1.d0
            x1 = w6js(j1,j9,k,j8,j4,j7)
        x2 = w6js(j2,j6,k,j4,j8,j5)
        x3 = w6js(j1,j9,k,j6,j2,j3)
        x = x + s*x1*x2*x3*dfloat(k+1)
        enddo
      w9js = x
      return
   end function w9js

!----------------------------------------------------------------
! Reorder the Regge square to make it similar to eq. (2.11) in Rasch & Yu (SIAM J. Sci. Comput. Vol 25, n 4, 1416)
!----------------------------------------------------------------
   function reorder_regge(regge, factor)
    integer :: reorder_regge(3,3), regge(3,3), temp(3,3), loc_max(2), col_max, row_max, loc_min(2), col_min, row_min, t
    integer :: i, j, tot, minus1, factor

        tot = sum(regge(1,:))
        minus1 = (-1)**tot
        factor = 1
        
        temp = regge
                
        loc_max = maxloc(temp)
        row_max = loc_max(2)
        col_max = loc_max(1)
        
        loc_min = minloc(temp)
        row_min = loc_min(2)
        col_min = loc_min(1)
        
! Put the maximum and minimum in the same row
        if (col_min == col_max) then
            temp = transpose(regge)
            
            row_max = loc_max(1)
            col_max = loc_max(2)

            row_min = loc_min(1)
            col_min = loc_min(2)

        endif
        
        reorder_regge = temp
        if (row_min /= 1) then
            reorder_regge(:,1) = temp(:,row_min)
            reorder_regge(:,row_min) = temp(:,1)
            factor = factor * minus1
        endif
        temp = reorder_regge
        
        loc_min = minloc(temp)
        row_min = loc_min(2)
        col_min = loc_min(1)
                        
! Put smallest one into R_11
        reorder_regge = temp
        if (col_min /= 1) then
            reorder_regge(1,:) = temp(col_min,:)
            reorder_regge(col_min,:) = temp(1,:)
            factor = factor * minus1                
            temp = reorder_regge
        endif
                
! Put largest one into R_21
        loc_max = maxloc(temp)
        row_max = loc_max(2)
        col_max = loc_max(1)
        
        if (col_max /= 2) then
            reorder_regge = temp
            reorder_regge(2,:) = temp(col_max,:)
            reorder_regge(col_max,:) = temp(2,:)
            factor = factor * minus1
        endif
        
        
   end function reorder_regge

!----------------------------------------------------------------
! Calculate the 3j symbol by indexing a 1D array with all the 3j symbols Rasch & Yu (SIAM J. Sci. Comput. Vol 25, n 4, 1416)
!----------------------------------------------------------------   
    function w3js_regge(j1,j2,j3,m1,m2,m3)
    integer :: j1, j2, j3, m1, m2, m3, regge(3,3), regge_order(3,3), L, X, T, B, S, c1, c2, c3, c4, c5, c, factor
    real(kind=8) :: w3js_regge
        
        w3js_regge = 0.d0
        if (abs(m1) > j1 .or. abs(m2) > j2 .or. abs(m3) > j3 .or. (m1+m2+m3) /= 0) return
                
        regge(1,1) = -j1+j2+j3
        regge(2,1) = j1-j2+j3
        regge(3,1) = j1+j2-j3
    
        regge(1,2) = j1-m1
        regge(2,2) = j2-m2
        regge(3,2) = j3-m3
    
        regge(1,3) = j1+m1
        regge(2,3) = j2+m2
        regge(3,3) = j3+m3
        
        regge_order = reorder_regge(regge,factor)
        
        L = regge_order(2,1)
        X = regge_order(1,2)
        T = regge_order(3,3)
        B = regge_order(2,2)
        S = regge_order(1,1)
        
        c1 = L*(24+L*(50+L*(35+L*(10+L))))/120
        c2 = X*(6+X*(11+X*(6+X)))/24
        c3 = T*(2+T*(3+T))/6
        c4 = B*(B+1)/2
        c5 = S
        c = c1+c2+c3+c4+c5+1

        w3js_regge = factor*threej(c)
    end function w3js_regge 
    
!----------------------------------------------------------------
! Calculate the list of 3j symbols and store them in a 1D array
!----------------------------------------------------------------   
    subroutine init_regge(Lmax)
    integer :: S, B, T, X, L, Lmax, cmax, c1, c2, c3, c4, c5, c, k, j1, j2, j3, m1, m2, m3
        
        cmax = Lmax*(274+Lmax*(225+Lmax*(85+Lmax*(15+Lmax))))/120+1

        allocate(threej(0:cmax))

        k = 0
        do L = 0, Lmax
            do X = 0, L
                do T = 0, X
                    do B = 0, T
                        do S = 0, B
                            j1 = (B+L-T+X)
                            j2 = (B+S-T+X)
                            j3 = (L+S)
                            m1 = (B+L-T-X)
                            m2 = (-B+S-T+X)
                            m3 = (T+T-L-S)

                            c1 = L*(24+L*(50+L*(35+L*(10+L))))/120
                            c2 = X*(6+X*(11+X*(6+X)))/24
                            c3 = T*(2+T*(3+T))/6
                            c4 = B*(B+1)/2
                            c5 = S
                            c = c1+c2+c3+c4+c5+1

                            threej(c) = w3js(j1,j2,j3,m1,m2,m3)

                            k = k + 1
                        enddo
                    enddo
                enddo
            enddo
        enddo

        open(unit=12,file='regge_3j.dat',action='write',form='unformatted')
        write(12) threej
        close(12)
        
    end subroutine init_regge

!----------------------------------------------------------------
! This function calculates the rotation matrix of order k
! The results are returned in dr(iq,iqp) and di(iq,iqp) representing
! the real and imaginary part of the rotation matrix respectevely.
! The angles alpha, beta and gamma are given in radians
!----------------------------------------------------------------
   subroutine rotmatk(k,alpha,beta,gamma,dr,di)
    integer :: k
    real(kind=8) :: alpha, beta, gamma
   real(kind=8) :: dr(-1:1,-1:1), di(-1:1,-1:1)
    integer :: iq, iqp, it, it1, it2, it3, itmax, itmin, ie1, ie2
    real(kind=8) :: a1, a2, ca, sa, cg, sg, bo2, cb2, sb2, redm, s, den

        dr = 0.d0
        di = 0.d0
        do iq = -k, k
            do iqp = -k, k
                a1 = alpha * dfloat(iq)
                a2 = gamma * dfloat(iqp)
                ca = dcos(a1)
                sa = dsin(a1)
                cg = dcos(a2)
                sg = dsin(a2)
                bo2 = beta / 2.d0
                cb2 = dcos(bo2)
                sb2 = dsin(bo2)
                it1 = k + iq
                it2 = k - iqp
                itmax = it1
                if (it2 < itmax) itmax = it2
                itmin = 0
                it3 = iq - iqp
                if (it3 > 0) itmin = it3
                redm = 0.d0
                if (beta /= 0.d0) then
                    do it = itmin, itmax
                        s = 1.d0
                        if (mod(it,2) /= 0) s = -1.d0
                        ie1 = 2*(k-it)+iq-iqp
                        ie2 = 2*it+iqp-iq
                        den = fact(it1-it) * fact(it2-it) * fact(it) * fact(it-it3)
                        redm = redm + s*(cb2**ie1)*(sb2**ie2) / den
                    enddo
                    redm = redm*sqrt(fact(k+iq)*fact(k-iq)*fact(k+iqp)*fact(k-iqp))
                else
                    if (iq == iqp) redm = 1.d0
                endif
                dr(iq,iqp) = redm*(ca*cg-sa*sg)
                di(iq,iqp) = -redm*(ca*sg+sa*cg)
            enddo               
        enddo

    end subroutine rotmatk


!----------------------------------------------------------------
! This function calculates the reduced rotation matrix of order 0,1,2,...kmax
! The results are returned in complex format dr(K,iq,iqp) 
! The angle beta is given in radians
! The maximum value of kmax is 6
!----------------------------------------------------------------       
   subroutine reduced_matrix(kmax,beta,reduced)
    integer :: kmax, K, M, N, t, min_t, max_t
    real(kind=8) :: beta, cos_beta, sin_beta, suma, f1, f2
    real(kind=8) :: reduced(0:kmax,-kmax:kmax,-kmax:kmax)

        cos_beta = cos(beta/2.d0)
      sin_beta = sin(beta/2.d0)

        do K = 0, kmax
            do M = -K, K
                do N = -K, K
                    max_t = min(K+M,K-N)
                    min_t = max(0,M-N)
                    suma = 0.d0
                    f1 = sqrt(fact(K+M)*fact(K-M)*fact(K+N)*fact(K-N))
                    do t = min_t, max_t                         
                        f2 = fact(K+M-t)*fact(K-N-t)*fact(t)*fact(t+N-M)
                        suma = suma + (-1.d0)**t*cos_beta**(2*K+M-N-2*t)*sin_beta**(2*t-M+N) / f2
                    enddo                       
                    reduced(k,M,N) = f1 * suma
                enddo
            enddo               
        enddo


    end subroutine reduced_matrix

!----------------------------------------------------------------
! This function calculates the rotation matrix of order 0,1,2,...kmax
! The results are returned in complex format dr(K,iq,iqp) 
! The angles alpha, beta and gamma are given in radians
! The maximum value of kmax is 6
!----------------------------------------------------------------       
   subroutine rotmatsu(kmax,alpha,beta,gamma,dr)
   integer :: kmax, K, M, N
    real(kind=8) :: alpha, beta, gamma, ca, sa, cg, sg
    complex(kind=8) :: dr(0:kmax,-kmax:kmax,-kmax:kmax), ii
    real(kind=8) :: reduced(0:kmax,-kmax:kmax,-kmax:kmax)

        ii = cmplx(0.d0,1.d0)
        call reduced_matrix(kmax,beta,reduced)
        
        do K = 0, kmax
            do M = -K, K
                do N = -K, K
                    dr(K,M,N) = reduced(K,M,N) * exp(-ii*(alpha*M+gamma*N))
                enddo
            enddo
        enddo

   end subroutine 

!----------------------------------------------------------------
! This function returns the Kronecker delta
!----------------------------------------------------------------       
    function delta(A,B)
    integer :: A, B, delta
        delta = 0.d0
        if (A == B) delta = 1.d0
    end function delta      

!----------------------------------------------------------------
! This function calculates the t^K_P(i) symbol
!----------------------------------------------------------------       
   function tkp(K,P)
    integer :: K, P
    complex(kind=8) :: tkp(0:3), ii

        tkp = 0.d0
        ii = cmplx(0.d0,1.d0)
        if (K == 0) then
            tkp(0) = delta(P,0)
            return
        endif

        if (K == 1) then
            tkp(3) = sqrt(3.d0/2.d0) * delta(P,0)
            return
        endif

        if (K == 2) then
            tkp(0) = sqrt(1.d0/2.d0) * delta(P,0)
            tkp(1) = -sqrt(3.d0)/2.d0 * (delta(P,-2) + delta(P,2))
            tkp(2) = ii * sqrt(3.d0)/2.d0 * (delta(P,-2) - delta(P,2))
            return
        endif                       
    end function tkp

!----------------------------------------------------------------
! This function calculates the tensor T^K_P(i,Omega) taking into account a double rotation
!----------------------------------------------------------------       
   subroutine tkq2(alpha,beta,gamma,alphap,betap,gammap,tr)
    real(kind=8) :: alpha, beta, gamma, alphap, betap, gammap
    complex(kind=8) :: tr(0:2,-2:2,0:3), tkpi(0:3), dr1(0:2,-2:2,-2:2), dr2(0:2,-2:2,-2:2), suma
    integer :: i, K, Q, P, Qp

        call rotmatsu(2,alpha,beta,gamma,dr1)
        call rotmatsu(2,alphap,betap,gammap,dr2)
        
        tr = 0.d0

        do i = 0, 3
            do K = 0, 2
                do Q = -K, K
                    suma = 0.d0
                    do P = -K, K
                        tkpi = tkp(K,P)
                        do Qp = -K, K
                            suma = suma + tkpi(i) * dr1(K,P,Qp) * dr2(K,Qp,Q)
                        enddo
                    enddo
                    tr(K,Q,i) = suma
                enddo
            enddo
        enddo           

    end subroutine tkq2

!----------------------------------------------------------------
! Calculates part of eq. 7.41 from Landi Degl'Innocenti & Landolfi (2004)
!----------------------------------------------------------------                   
   subroutine bigpar(l2,is2,k2,kp2,j2,jp2,js2,jt2,x)
    integer :: l2, is2, k2, kp2, j2, jp2, js2, jt2
    real(kind=8) :: x, x1, x2, x3, g

      x = 0.d0
      if (jp2 /= jt2 .and. j2 /= js2) return
      if (jp2 == jt2) then
        if (iabs(j2-js2) > 2 .or. j2+js2 == 0) then
            else
            call gamma(l2,is2,j2,js2,g)
            x1 = g
            x2 = w6js(k2,kp2,2,js2,j2,jp2)
            x = x1*x2
            endif
        endif

      if (j2 /= js2) return
      if (iabs(jt2-jp2) > 2 .or. jt2+jp2 == 0) return
      call gamma(l2,is2,jt2,jp2,g)
      x1 = g
      x2 = 1.d0
      if (mod(k2-kp2,4) /= 0) x2 = -1.d0
      x3 = w6js(k2,kp2,2,jt2,jp2,j2)
      x = x + x1*x2*x3
    end subroutine bigpar
      
!----------------------------------------------------------------
! Calculates eq. 7.42 from Landi Degl'Innocenti & Landolfi (2004)
!----------------------------------------------------------------           
   subroutine gamma(l2,is2,j2,jp2,g)
    integer :: l2, is2, j2, jp2
    real(kind=8) :: g
    real(kind=8) :: x1, x2, x3

      g = 0.d0
      if (j2 == jp2) then
            g = g + dsqrt(dfloat(j2*(j2+2)*(j2+1))/4.d0)
        endif
      x1 = 1.d0
      if (mod(2+l2+is2+j2,4) /= 0) x1 = -1.d0
      x2 = dsqrt(dfloat((j2+1)*(jp2+1)*is2*(is2+2)*(is2+1))/4.d0)
      x3 = w6js(j2,jp2,2,is2,is2,l2)
      g = g + x1*x2*x3
    end subroutine gamma

!----------------------------------------------------------------
! Calculates the factorial of the integers up to 301 and save it into the fact array
!----------------------------------------------------------------
   subroutine factrl
    integer :: i
      fact(0) = 1.d0
      do i=1,101
            fact(I) = fact(I-1) * dble(I)
        enddo
   END subroutine factrl
        
! ---------------------------------------------------------
! LU decomposition of a matrix
!  INPUT:
!       - a is the matrix to decompose
!       
!  OUTPUT:
!       - a is the LU decomposition of a
!       - indx is a vector that records the row permutation effected by the partial pivoting
!       - d takes values +1/-1 depending on whether the number of row interchanges was odd or even
! ---------------------------------------------------------
    subroutine ludcmp(a,indx,d,error)
    integer, INTENT(INOUT) :: indx(:)
    real*8, INTENT(INOUT) :: a(:,:), d
    real*8, parameter :: TINY = 1.d-20
    integer :: i, imax, j, k, n
    real*8 :: aamax, dum, sum, vv(size(a,1))
    integer :: values_start(8), values_end(8), error
    real(kind=8) :: start_time, end_time
    
!       call date_and_time(values=values_start)
        
        d = 1.d0
        n = size(a,1)
        error = 0
        
        do i = 1, n
            aamax = 0.d0    
            aamax = maxval(dabs(a(i,:)))
            if (aamax == 0.d0) then
                print *, 'Singular matrix in LU decomposition'
                error = 1
                return
            endif
            vv(i) = 1.d0 / aamax
        enddo
        
!$OMP PARALLEL DO PRIVATE(I,J,K,sum,aamax,imax,d,dum) SHARED (A,indx)
        do j = 1, n
            do i = 1, j-1
                sum = a(i,j)
                do k = 1, i-1
                    sum = sum - a(i,k) * a(k,j)
                enddo
                a(i,j) = sum
            enddo
            aamax = 0.d0
            do i = j, n
                sum = a(i,j)
                do k = 1, j-1
                    sum = sum - a(i,k) * a(k,j)
                enddo
                a(i,j) = sum
                dum = vv(i) * dabs(sum)
                if (dum >= aamax) then
                    imax = i
                    aamax = dum
                endif               
            enddo
            if (j /= imax) then
                do k = 1, n
                    dum = a(imax,k)
                    a(imax,k) = a(j,k)
                    a(j,k) = dum
                enddo
                d = -d
                vv(imax) = vv(j)
            endif
            indx(j) = imax
            if (a(j,j) == 0.d0) a(j,j) = TINY
            if (j /= n) then
                dum = 1.d0 / a(j,j)
                do i = j+1, n
                    a(i,j) = a(i,j) * dum
                enddo
            endif
        enddo
!$OMP END PARALLEL DO

!       call date_and_time(values=values_end)
!       start_time = values_start(5) * 3600 + values_start(6) * 60 + values_start(7) + 0.001 * values_start(8)
!       end_time = values_end(5) * 3600 + values_end(6) * 60 + values_end(7) + 0.001 * values_end(8)
!       print *, 'Time : ', end_time-start_time
    
    end subroutine ludcmp
    
! ---------------------------------------------------------
! Solves the set of equations AX=b where A is the LU decomposition of a matrix
!  INPUT:
!       - a is the LU decomposition of the system matrix
!       - b is the right hand side vector of the system
!       - indx is the vector returned by ludcmp
!  OUTPUT:
!       - b is the solution of the system
! ---------------------------------------------------------
    subroutine lubksb(a,indx,b)
    real*8, INTENT(IN) :: a(:,:)
    real*8, INTENT(INOUT) :: b(:)
    integer, INTENT(IN) :: indx(:)
    integer :: i, ii, n, j, ll
    real*8 :: sum
        n = size(a,1)
        ii = 0
        do i = 1, n
            ll = indx(i)
            sum = b(ll)
            b(ll) = b(i)
            if (ii /= 0) then
                do j = ii, i-1
                    sum = sum - a(i,j) * b(j)
                enddo
            else if (sum /= 0.d0) then
                ii = i
            endif
            b(i) = sum
        enddo
        do i = n, 1, -1
            sum = b(i)
            do j = i+1, n
                sum = sum - a(i,j) * b(j)
            enddo
            b(i) = sum / a(i,i)
        enddo
    end subroutine lubksb       
    
! ---------------------------------------------------------
! This subroutine solves a linear system of equations using the SLAP routines
! It uses an Incomplete LU BiConjugate Gradient Sparse Ax=b solver
! ---------------------------------------------------------       
!   subroutine slapsolver(a,b)
!   real(kind=8), INTENT(INOUT) :: a(:,:), b(:)
!   integer :: i, j, k, j2
!   real(kind=8), allocatable :: a_sparse(:), x(:), rwork(:)
!   integer, allocatable :: ia(:), ja(:), iwork(:)
!   integer :: n, nelt, isym, itol, itmax, iter, ierr, iunit, lenw, leniw, nsave
!   real(kind=8) :: tol, err
!   integer, allocatable :: indx(:)
!   
!       n = size(b)
!       allocate(x(n))
!       k = 0
!       do i = 1, n
!           do j = 1, n
!               if (a(j,i) /= 0.d0) k = k + 1
!           enddo
!       enddo
!       nelt = k
!       if (verbose_mode == 1) then
!           print *, 'Number of non-zero elements : ', nelt
!       endif
!       allocate(a_sparse(nelt))
!       allocate(ia(nelt))
!       allocate(ja(nelt))      
!       
!       k = 1
!       do i = 1, n
!           do j = 1, n
!               if (a(j,i) /= 0.d0) then
!                   a_sparse(k) = a(j,i)
!                   ia(k) = j
!                   ja(k) = i
!                   k = k + 1
!               endif
!           enddo
!       enddo
!       
!       x = 0.d0
!       
!       isym = 0
!       itol = 1
!       tol = 1.d-10
!       itmax = 10000
!       iunit = 0
!       lenw = n*n
!       leniw = n*n
!       nsave = 10
!       allocate(iwork(leniw))
!       allocate(rwork(lenw))
!       
!       if (verbose_mode == 1) then 
!           print *, 'Solving the system...'
!       endif
!       call ds2y( n, nelt, ia, ja, a_sparse, isym )
!       call dslubc(n, b, x, nelt, ia, ja, a_sparse, isym, itol, tol, &
!           itmax, iter, err, ierr, iunit, rwork, lenw, iwork, leniw )
!       if (verbose_mode == 1) then 
!           print *, 'Number of iterations : ', iter, 'of ', itmax
!           print *, 'Error : ', err, ierr
!       endif
!       if (iter >= itmax+1 .or. ierr /= 0) then
!           if (verbose_mode == 1) then 
!               print *, 'Couldn''t reach convergence. Trying other initialization...'
!           endif
!           write(17,*) 'Couldn''t reach convergence. Trying other initialization...'
!           x = 1.d0
!           isym = 0
!           itol = 0
!           tol = 1.d-10
!           itmax = 10000
!           iunit = 0
!           lenw = n*n
!           leniw = n*n
!           nsave = 10
!           call dslucs(n, b, x, nelt, ia, ja, a_sparse, isym, itol, tol, &
!           itmax, iter, err, ierr, iunit, rwork, lenw, iwork, leniw )
!           if (verbose_mode == 1) then 
!               print *, 'Number of iterations : ', iter, 'of ', itmax
!               print *, 'Error : ', err, ierr
!           endif
!       endif
!       
!       if (iter >= itmax+1 .or. ierr /= 0) then
!           if (verbose_mode == 1) then 
!               print *, 'Couldn''t reach convergence. Trying GMRES...'
!           endif
!           write(17,*) 'Couldn''t reach convergence. Trying GMRES...'
!           x = 0.d0
!           isym = 0
!           itol = 0
!           tol = 1.d-10
!           itmax = 10000
!           iunit = 0
!           lenw = n*n
!           leniw = n*n
!           nsave = 10
!           call dsdgmr(n, b, x, nelt, ia, ja, a_sparse, isym, nsave, itol, tol, &
!           itmax, iter, err, ierr, iunit, rwork, lenw, iwork, leniw )
!           if (verbose_mode == 1) then 
!               print *, 'Number of iterations : ', iter, 'of ', itmax
!               print *, 'Error : ', err, ierr
!           endif
!       endif
!       
!       if (iter >= itmax+1 .or. ierr /= 0) then
!           if (verbose_mode == 1) then 
!               print *, 'Couldn''t reach convergence. Trying other initialization...'
!           endif
!           write(17,*) 'Couldn''t reach convergence. Trying other initialization...'
!           x = 1.d0
!           isym = 0
!           itol = 0
!           tol = 1.d-10
!           itmax = 10000
!           iunit = 0
!           lenw = n*n
!           leniw = n*n
!           nsave = 10
!           call dsdgmr(n, b, x, nelt, ia, ja, a_sparse, isym, nsave, itol, tol, &
!           itmax, iter, err, ierr, iunit, rwork, lenw, iwork, leniw )
!           if (verbose_mode == 1) then 
!               print *, 'Number of iterations : ', iter, 'of ', itmax
!               print *, 'Error : ', err, ierr
!           endif
!       endif
!       
!       b = x
!       
!       deallocate(x)
!       deallocate(a_sparse)
!       deallocate(ia)
!       deallocate(ja)
!       deallocate(iwork)
!       deallocate(rwork)
!       
!       if (iter >= itmax+1 .or. ierr /= 0) then
!           if (verbose_mode == 1) then 
!               print *, 'Couldn''t reach convergence. Trying LU decomposition...'
!           endif
!           write(17,*) 'Couldn''t reach convergence. Trying LU decomposition...'
!           allocate(indx(n))
!           call ludcmp(a,indx,tol)
!           call lubksb(a,indx,b)
!           deallocate(indx)
!       endif                                           
!       
!   end subroutine slapsolver

! ---------------------------------------------------------
! Read lines that don't mind in a file
! INPUT:
!     - unit : unit to whom the file is associated
!     - n : number of lines to read
! ---------------------------------------------------------
   subroutine lb(unit, n)
   integer, INTENT(IN) :: unit, n
   integer :: i
          do i = 1, n
                     read(unit,*)
          enddo
   end subroutine lb
    
! ---------------------------------------------------------
! Returns the sizes of the eigenvectors of the Paschen Back effect
! ---------------------------------------------------------
    subroutine paschen_size(l2,is2,j2max,njlargest)
    integer :: l2, is2
    integer :: m2, j2max, j2min, j2a, j2b, nj, njlargest
    
        j2min = abs(l2-is2)
        j2max = l2 + is2
                
! First find the largest value of nj
        njlargest = 0
        do m2 = -j2max, j2max, 2
            j2a = abs(m2)
            if (j2a < j2min) j2a = j2min
            j2b = j2max
            nj = 1 + (j2b-j2a) / 2
            if (nj > njlargest) njlargest = nj
        enddo
    end subroutine paschen_size
            
! ---------------------------------------------------------
! Calculates the eigenvalues and eigenvectors of the Paschen-Back matrix
! INPUT:
!     l2 = 2*L
!     is2 = 2*S
!     e0 = a vector containing the energy in cm^-1 of the level J (Jmin <= J <= Jmax) organized in terms of the J2 value
!          ex. e0(3) is the energy of the J=3/2 level. (Some elements remain undefined)
!     b = magnetic field in gauss
! OUTPUT:
!     njlev = number of levels for each m2
!     e = a vector containing the eigenvalues of the matrix. It is organized as e(m2,jsmall), where m2 is 2*M and
!         jsmall is a quantum number which takes values between 1 and the number of J compatible with M
!     c = eigenvectors, organized as c(m2,jsmall,J2)
! --------------------------------------------------------- 
    subroutine paschen(l2,is2,e0,B,njlev,e,c,j2m,njl)
    integer :: l2, is2, j2m, njl
    integer :: njlev(-j2m:j2m)
    real(kind=8) :: e(-j2m:j2m,njl), c(-j2m:j2m,njl,0:j2m)
    real(kind=8) :: b, e0(0:jlimit2), bb
    real(kind=8), allocatable :: d(:), ud(:), z(:,:)
    integer :: m2, j2max, j2min, j2a, j2b, nj, j2vero, j, jp, jp2vero, njlargest, i, jsmall, jaux, j2
    real(kind=8) :: x1, x2, x3, x4
    
        bb = PE / (4.d0*PI*PME*PC**2) * B
        j2min = abs(l2-is2)
        j2max = l2 + is2
                        
        do m2 = -j2max, j2max, 2

! Build the matrix which corresponds to the value of M
            j2a = abs(m2)
            if (j2a < j2min) j2a = j2min
            j2b = j2max
            nj = 1 + (j2b-j2a) / 2
            njlev(m2) = nj
            
            allocate(d(nj))
            d = 0.d0
            allocate(ud(nj))
            ud = 0.d0
            allocate(z(nj,nj))
                        
! First the diagonal part
            do j = 1, nj
                j2vero = j2a+2*(j-1)
            d(j) = e0(j2vero) 
            x1 = 1.d0
            if (mod(2*j2vero+l2+is2+m2,4) /= 0) x1 = -1.d0
            x2 = dfloat(j2vero+1)*dsqrt(dfloat(is2*(is2+2)*(is2+1))/4.d0)
            x3 = w3js(j2vero,j2vero,2,-m2,m2,0)
            x4 = w6js(j2vero,j2vero,2,is2,is2,l2)
            d(j) = d(j) + bb * (dfloat(m2)/2.d0+x1*x2*x3*x4)
            enddo
                
! If we use the Pashen-Back regime
            if (use_paschen_back == 1) then
! Then the upper diagonal (only if the matrix is at least 2x2)
                if (nj /= 1) then
                    do jp = 2, nj
                    jp2vero = j2a+2*(jp-1)
                    j2vero = jp2vero-2
                    x1 = 1.d0
                    if (mod(j2vero+jp2vero+l2+is2+m2,4) /= 0) x1 = -1.d0
                    x2 = dsqrt(dfloat((j2vero+1)*(jp2vero+1)*is2*(is2+2)*(is2+1))/4.d0)
                    x3 = w3js(j2vero,jp2vero,2,-m2,m2,0)
                    x4 = w6js(j2vero,jp2vero,2,is2,is2,l2)
                    ud(jp) = bb*x1*x2*x3*x4
                    enddo
                endif
            endif
            
! Initializes the matrix z with the eigenvectors
            z = 0.d0
            do i = 1, nj
                z(i,i) = 1.d0
            enddo
            
! Diagonalize the matrix
            call tqli(d,ud,nj,nj,z)
            
            do jsmall = 1, nj
                e(m2,jsmall) = d(jsmall)
                do jaux = 1, nj
                    j2 = j2a + 2*(jaux-1)
                    c(m2,jsmall,j2) = z(jaux,jsmall)
                enddo
            enddo
            
            deallocate(d)
            deallocate(ud)
            deallocate(z)
            
        enddo
        
    end subroutine paschen
    
! ---------------------------------------------------------
! Diagonalizes a tridiagonal matrix
! ---------------------------------------------------------
    subroutine tqli(d,e,n,np,z)
    integer :: n, np, i, l, m, iter, k
    real(kind=8) :: d(np),e(np),z(np,np), g, r, s, c, p, f, b, dd
      
        if (n.gt.1) then
      do 11 i=2,n
      e(i-1)=e(i)
 11   continue
      e(n)=0.d0
      do 15 l=1,n
      iter=0
 1    do 12 m=l,n-1
      dd=dabs(d(m))+dabs(d(m+1))
      if(dabs(e(m))+dd.eq.dd) go to 2
 12   continue
      m=n
 2    if(m.ne.l) then
      if(iter.eq.30) then
        print *, 'too many iterations'
        stop
      endif
      iter=iter+1
      g=(d(l+1)-d(l))/(2.d0*e(l))
      r=dsqrt(g**2+1.d0)
      g=d(m)-d(l)+e(l)/(g+dsign(r,g))
      s=1.d0
      c=1.d0
      p=0.d0
      do 14 i=m-1,l,-1
      f=s*e(i)
      b=c*e(i)
      if(dabs(f).ge.dabs(g)) then
      c=g/f
      r=dsqrt(c**2+1.d0)
      e(i+1)=f*r
      s=1.d0/r
      c=c*s
      else
      s=f/g
      r=dsqrt(s**2+1.d0)
      e(i+1)=g*r
      c=1.d0/r
      s=s*c
      endif
      g=d(i+1)-p
      r=(d(i)-g)*s+2.d0*c*b
      p=s*r
      d(i+1)=g+p
      g=c*r-b
      do 13 k=1,n
      f=z(k,i+1)
      z(k,i+1)=s*z(k,i)+c*f
      z(k,i)=c*z(k,i)-s*f
 13   continue
 14   continue
      d(l)=d(l)-p
      e(l)=g
      e(m)=0.d0
      go to 1
      endif
 15   continue
      endif
      return
      end subroutine tqli
        
! ---------------------------------------------------------
! Sets the imaginary part of a complex to a given number
! ---------------------------------------------------------
        subroutine set_imaginary(A,b)
        complex(kind=8) :: A
        real(kind=8) :: b
            A = cmplx(real(A),b)
        end subroutine set_imaginary
        
! ---------------------------------------------------------
! Sets the imaginary part of a complex to a given number
! ---------------------------------------------------------
        subroutine set_real(A,b)
        complex(kind=8) :: A
        real(kind=8) :: b
            A = cmplx(b,aimag(A))
        end subroutine set_real
        


!--------------------------------------------------------------
! Dawson's integral
!--------------------------------------------------------------
    function dawson(dv)
    real(kind=8) :: dv(:)
    real(kind=8) :: dawson(size(dv)), c(6), x
    real(kind=8), parameter :: H=0.4, A1=2./3., A2=0.4, A3=2./7.
    integer :: i, n, j, n0
    real(kind=8) :: xx, xp, e1, e2, d1, d2, sum, x2
                
        n = size(dv)
        do i = 1, 6
            c(i)=exp(-((2.*float(i)-1.)*H)**2)
        enddo

        do i = 1, n
            if (abs(dv(i)) < 0.2) then
                x = dv(i)
                x2 = dv(i)**2
                dawson(i) = x*(1.-A1*x2*(1.-A2*x2*(1.-A3*x2)))
            else
                x = dv(i)
                xx = abs(dv(i))
                n0 = 2*nint(0.5*xx/H)
                xp = xx-float(n0)*H
                e1 = exp(2.*xp*H)
                e2 = e1**2
                d1 = float(n0+1)
                d2 = d1-2.
                sum = 0.
                do j = 1,6
                  sum = sum + c(j)*(e1/d1+1./(d2*e1))
                  d1 = d1+2.
                  d2 = d2-2.
                  e1 = e2*e1
                enddo
                dawson(i) = 0.5641895835*sign(exp(-xp**2),x)*sum
            endif
        enddo

    end function dawson

!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles using a different damping for each frequency
! point
! See Humlicek (1982) JQSRT 27, 437
!--------------------------------------------------------------
    function profile(da, dv)
    real(kind=8) :: da
    real(kind=8) :: dv(:)
    complex(kind=8) :: profile(size(dv))
    complex(kind=8) :: w4, z, t, u, v4
    real(kind=8) :: s
    integer :: i, n

                
        n = size(dv)

! For small a, evaluate the power expansion up to second order (Landi degl'Innocenti & Landolfi 2004, p. 168)
        if (da < 1.d-3) then
            profile = cmplx(exp(-dv**2) + 2.d0*da/sqrt(PI)*(2*dv*dawson(dv)-1.d0), 2.d0*dawson(dv) / sqrt(PI) - 2.d0*da*dv*exp(-dv**2))
            return
        endif
        do i = 1, n
! There is a problem for very small dampings. We force it to be 0 using a threshold
            z = cmplx(dv(i), da)
            t = cmplx(da, -dv(i))
            s = dabs(dv(i)) + da
            u = t*t


            if (s >= 15.d0) then
                w4 = t * 0.5641896d0 / (0.5d0+u)
            elseif (s >= 5.5) then
                w4 = t*(1.410474d0+u*0.5641896d0)/(0.75d0+u*(3.d0+u))
            elseif (da >= 0.195d0*dabs(dv(i))-0.176d0) then
                w4 = (16.4955d0+t*(20.20933d0+t*(11.96482d0+t*(3.778987d0+t*0.5642236d0)))) / &
                    (16.4955d0+t*(38.82363d0+t*(39.27121d0+t*(21.69274d0+t*(6.699398d0+t)))))
            else 
                w4 = t*(36183.31d0-u*(3321.9905d0-u*(1540.787d0-u*(219.0313d0-u*(35.76683d0-&
                    u*(1.320522d0-u*0.56419d0))))))
                v4 = (32066.6d0-u*(24322.84d0-u*(9022.228d0-u*(2186.181d0-u*(364.2191d0-&
                    u*(61.57037d0-u*(1.841439d0-u)))))))
                w4 = exp(u) - w4/v4
            endif
            profile(i) = w4
        enddo

    end function profile
    
! ---------------------------------------------------------
! Given xa(:) and ya(:) which tabulate a function, returns the interpolation using
! linear interpolation of vector x(:) in y(:)
! ---------------------------------------------------------     
    subroutine lin_interpol(xa,ya,x,y)
    real*8, INTENT(IN) :: xa(:), ya(:), x(:)
    real*8, INTENT(INOUT) :: y(:)
    integer :: i, n
    integer :: locat(1), loc

        n = size(x)

        do i = 1, n
            locat = minloc(dabs(xa-x(i)))
            loc = locat(1)
            if (loc == 1) then
                y(i) = ya(loc)
            else 
                y(i) = (ya(loc)-ya(loc-1))/(xa(loc)-xa(loc-1)) * (x(i)-xa(loc-1)) + ya(loc-1)
            endif
        enddo

    end subroutine lin_interpol 

! ---------------------------------------------------------
! Given etaI, etaQ, etaU, etaV, rhoQ, rhoU and rhoV, fill the absorption matrix
! ---------------------------------------------------------     
    subroutine fill_absorption_matrix(matrix,etaI,etaQ,etaU,etaV,rhoQ,rhoU,rhoV)
    real(kind=8) :: matrix(4,4), etaI, etaQ, etaU, etaV, rhoQ, rhoU, rhoV

        matrix(1,1) = etaI
        matrix(2,2) = etaI
        matrix(3,3) = etaI
        matrix(4,4) = etaI
        matrix(1,2) = etaQ
        matrix(2,1) = etaQ          
        matrix(1,3) = etaU
        matrix(3,1) = etaU
        matrix(1,4) = etaV
        matrix(4,1) = etaV           
        matrix(2,3) = rhoV
        matrix(3,2) = -rhoV
        matrix(2,4) = -rhoU
        matrix(4,2) = rhoU                       
        matrix(3,4) = rhoQ           
        matrix(4,3) = -rhoQ  
        
    end subroutine fill_absorption_matrix

!--------------------------------------------------------------
! Inversion of a 4x4 matrix
!--------------------------------------------------------------
    subroutine invert(a)
    real*8 :: a(4,4)
    real*8 :: b(4,4), det, maxim, fabsmax
    ! First some tests of singularity
        b = dabs(a)
        maxim = maxval(b)
        fabsmax = 1.d0 / maxim
        if (maxim == 0.d0) then
            print *, 'Singularity in the inversion'
!           stop
        endif

        a = a * fabsmax

    b(1,1) = a(2,2) * a(3,3) * a(4,4) + a(2,3) * a(3,4) * a(4,2)&
        + a(2,4) * a(3,2) * a(4,3) - a(2,2) * a(3,4) * a(4,3)&
        - a(2,3) * a(3,2) * a(4,4) - a(2,4) * a(3,3) * a(4,2)
    b(2,1) = a(2,3) * a(3,1) * a(4,4) + a(2,4) * a(3,3) * a(4,1)&
        + a(2,1) * a(3,4) * a(4,3) - a(2,3) * a(3,4) * a(4,1)&
        - a(2,4) * a(3,1) * a(4,3) - a(2,1) * a(3,3) * a(4,4)
    b(3,1) = a(2,4) * a(3,1) * a(4,2) + a(2,1) * a(3,2) * a(4,4)&
        + a(2,2) * a(3,4) * a(4,1) - a(2,4) * a(3,2) * a(4,1)&
        - a(2,1) * a(3,4) * a(4,2) - a(2,2) * a(3,1) * a(4,4)
    b(4,1) = a(2,1) * a(3,3) * a(4,2) + a(2,2) * a(3,1) * a(4,3)&
        + a(2,3) * a(3,2) * a(4,1) - a(2,1) * a(3,2) * a(4,3)&
        - a(2,2) * a(3,3) * a(4,1) - a(2,3) * a(3,1) * a(4,2)
    b(1,2) = a(3,2) * a(4,4) * a(1,3) + a(3,3) * a(4,2) * a(1,4)&
        + a(3,4) * a(4,3) * a(1,2) - a(3,2) * a(4,3) * a(1,4)&
        - a(3,3) * a(4,4) * a(1,2) - a(3,4) * a(4,2) * a(1,3)
    b(2,2) = a(3,3) * a(4,4) * a(1,1) + a(3,4) * a(4,1) * a(1,3)&
        + a(3,1) * a(4,3) * a(1,4) - a(3,3) * a(4,1) * a(1,4)&
        - a(3,4) * a(4,3) * a(1,1) - a(3,1) * a(4,4) * a(1,3)
    b(3,2) = a(3,4) * a(4,2) * a(1,1) + a(3,1) * a(4,4) * a(1,2)&
        + a(3,2) * a(4,1) * a(1,4) - a(3,4) * a(4,1) * a(1,2)&
        - a(3,1) * a(4,2) * a(1,4) - a(3,2) * a(4,4) * a(1,1)
    b(4,2) = a(3,1) * a(4,2) * a(1,3) + a(3,2) * a(4,3) * a(1,1)&
        + a(3,3) * a(4,1) * a(1,2) - a(3,1) * a(4,3) * a(1,2)&
        - a(3,2) * a(4,1) * a(1,3) - a(3,3) * a(4,2) * a(1,1)
    b(1,3) = a(4,2) * a(1,3) * a(2,4) + a(4,3) * a(1,4) * a(2,2)&
        + a(4,4) * a(1,2) * a(2,3) - a(4,2) * a(1,4) * a(2,3)&
        - a(4,3) * a(1,2) * a(2,4) - a(4,4) * a(1,3) * a(2,2)
    b(2,3) = a(4,3) * a(1,1) * a(2,4) + a(4,4) * a(1,3) * a(2,1)&
        + a(4,1) * a(1,4) * a(2,3) - a(4,3) * a(1,4) * a(2,1)&
        - a(4,4) * a(1,1) * a(2,3) - a(4,1) * a(1,3) * a(2,4)
    b(3,3) = a(4,4) * a(1,1) * a(2,2) + a(4,1) * a(1,2) * a(2,4)&
        + a(4,2) * a(1,4) * a(2,1) - a(4,4) * a(1,2) * a(2,1)&
        - a(4,1) * a(1,4) * a(2,2) - a(4,2) * a(1,1) * a(2,4)
    b(4,3) = a(4,1) * a(1,3) * a(2,2) + a(4,2) * a(1,1) * a(2,3)&
        + a(4,3) * a(1,2) * a(2,1) - a(4,1) * a(1,2) * a(2,3)&
        - a(4,2) * a(1,3) * a(2,1) - a(4,3) * a(1,1) * a(2,2)
    b(1,4) = a(1,2) * a(2,4) * a(3,3) + a(1,3) * a(2,2) * a(3,4)&
        + a(1,4) * a(2,3) * a(3,2) - a(1,2) * a(2,3) * a(3,4)&
        - a(1,3) * a(2,4) * a(3,2) - a(1,4) * a(2,2) * a(3,3)
    b(2,4) = a(1,3) * a(2,4) * a(3,1) + a(1,4) * a(2,1) * a(3,3)&
        + a(1,1) * a(2,3) * a(3,4) - a(1,3) * a(2,1) * a(3,4)&
        - a(1,4) * a(2,3) * a(3,1) - a(1,1) * a(2,4) * a(3,3)
    b(3,4) = a(1,4) * a(2,2) * a(3,1) + a(1,1) * a(2,4) * a(3,2)&
        + a(1,2) * a(2,1) * a(3,4) - a(1,4) * a(2,1) * a(3,2)&
        - a(1,1) * a(2,2) * a(3,4) - a(1,2) * a(2,4) * a(3,1)
    b(4,4) = a(1,1) * a(2,2) * a(3,3) + a(1,2) * a(2,3) * a(3,1)&
        + a(1,3) * a(2,1) * a(3,2) - a(1,1) * a(2,3) * a(3,2)&
        - a(1,2) * a(2,1) * a(3,3) - a(1,3) * a(2,2) * a(3,1)

        det = a(1,1) * b(1,1) + a(1,2) * b(2,1) + a(1,3) * b(3,1) + a(1,4) * b(4,1)

        a = b * (fabsmax / det)

    end subroutine invert

!--------------------------------------------------------------
! Returns in b the evolution operator if the propagation matrix is constant
! given by exp(-A*t). It uses Eqs. A5.20 of Landi Degl'Innocenti & Landolfi 2004
!--------------------------------------------------------------
    subroutine evol_operator(a,t,b)
    real(kind=8) :: a(4,4), t
    real(kind=8) :: b(4,4), det, maxim, fabsmax, w(3), y(3), w_y_scalar, sigma
    real(kind=8) :: Lambda_plus, Lambda_minus, tt, Theta, w_mod_sq, y_mod_sq, factor1, factor2, factor3, factor4, w0, temp1, temp2
    real(kind=8), dimension(4,4) :: M1, M2, M3, M4
    integer :: i

        w0 = a(1,1)

        w(1) = a(1,2)
        w(2) = a(1,3)
        w(3) = a(1,4)
        w_mod_sq = sum(w*w)

        y(1) = a(3,4)
        y(2) = a(4,2)
        y(3) = a(2,3)
        y_mod_sq = sum(y*y)

        w_y_scalar = sum(w*y)
        sigma = 1.d0
        if (w_y_scalar < 0.d0) sigma = -1.d0

! Eigenvalues
        temp1 = sqrt((w_mod_sq - y_mod_sq)**2 / 4.d0 + w_y_scalar**2)
        temp2 = (w_mod_sq - y_mod_sq) / 2.d0
        if (temp1+temp2 < 0.d0) then
            Lambda_plus = 0.d0
        else
            Lambda_plus = sqrt(temp1 + temp2)
        endif
        if (temp1-temp2 < 0.d0) then
            Lambda_minus = 0.d0
        else
            Lambda_minus = sqrt(temp1 - temp2)
        endif
        
        Theta = 2.d0 * temp1
        tt = (w_mod_sq + y_mod_sq) / 2.d0

        M1 = 0.d0
        M2 = 0.d0
        M3 = 0.d0
        M4 = 0.d0

! M1
        do i = 1, 4
            M1(i,i) = 1.d0
        enddo

! M2
        M2(1,2) = Lambda_minus * w(1) - sigma * Lambda_plus * y(1)
        M2(1,3) = Lambda_minus * w(2) - sigma * Lambda_plus * y(2)
        M2(1,4) = Lambda_minus * w(3) - sigma * Lambda_plus * y(3)
        M2(2,1) = M2(1,2)
        M2(2,3) = sigma * Lambda_plus * w(3) + Lambda_minus * y(3)
        M2(2,4) = -sigma * Lambda_plus * w(2) - Lambda_minus * y(2)
        M2(3,1) = M2(1,3)
        M2(3,2) = -M2(2,3)
        M2(3,4) = sigma * Lambda_plus * w(1) + Lambda_minus * y(1)
        M2(4,1) = M2(1,4)
        M2(4,2) = -M2(2,4)
        M2(4,3) = -M2(3,4)

! M3
        M3(1,2) = Lambda_plus * w(1) + sigma * Lambda_minus * y(1)
        M3(1,3) = Lambda_plus * w(2) + sigma * Lambda_minus * y(2)
        M3(1,4) = Lambda_plus * w(3) + sigma * Lambda_minus * y(3)
        M3(2,1) = M3(1,2)
        M3(2,3) = -sigma * Lambda_minus * w(3) + Lambda_plus * y(3)
        M3(2,4) = sigma * Lambda_minus * w(2) - Lambda_plus * y(2)
        M3(3,1) = M3(1,3)
        M3(3,2) = -M3(2,3)
        M3(3,4) = -sigma * Lambda_minus * w(1) + Lambda_plus * y(1)
        M3(4,1) = M3(1,4)
        M3(4,2) = -M3(2,4)
        M3(4,3) = -M3(3,4)

! M4
        M4(1,1) = tt
        M4(1,2) = w(3)*y(2) - w(2)*y(3)
        M4(1,3) = w(1)*y(3) - w(3)*y(1)
        M4(1,4) = w(2)*y(1) - w(1)*y(2)
        M4(2,1) = -M4(1,2)
        M4(2,2) = w(1)**2 + y(1)**2 - tt
        M4(2,3) = w(1)*w(2) + y(1)*y(2)
        M4(2,4) = w(3)*w(1) + y(3)*y(1)
        M4(3,1) = -M4(1,3)
        M4(3,2) = M4(2,3)
        M4(3,3) = w(2)**2 + y(2)**2 - tt
        M4(3,4) = w(2)*w(3) + y(2)*y(3)
        M4(4,1) = -M4(1,4)
        M4(4,2) = M4(2,4)
        M4(4,3) = M4(3,4)
        M4(4,4) = w(3)**2 + y(3)**2 - tt

        M2 = M2 / Theta
        M3 = M3 / Theta
        M4 = 2.d0 * M4 / Theta

        factor1 = 0.5d0*(cosh(t*Lambda_plus) + cos(t*Lambda_minus))
        factor2 = sin(t*Lambda_minus)
        factor3 = sinh(t*Lambda_plus)
        factor4 = 0.5d0*(cosh(t*Lambda_plus) - cos(t*Lambda_minus))

        b = factor1*M1 - factor2*M2 - factor3*M3 + factor4*M4
        b = exp(-t*w0) * b
        
    end subroutine evol_operator

! ---------------------------------------------------------
! Return the input array but shifted by a sub-pixel amount
! Note that it works only for any length of series which has no prime factor greater than 23
! ---------------------------------------------------------
        function fft_shift(x, sh)
        real(kind=8) :: x(:), fft_shift(size(x)), sh
        complex(kind=8) :: k(size(x)), ff(size(x))
        integer :: nx, i, n21
        
            nx = size(x)

            do i = 1, nx
            k(i) = i-1
         enddo
         n21 = nx/2+1
         do i = n21+1, nx
            k(i) = n21-nx+(i-1-n21)
         enddo

         k = 2.d0*PI*k/(1.d0*nx)
         k = k * cmplx(0.d0,1.d0)

         ff = x

         ff = fft(ff)

         fft_shift = real(fft(ff * exp(k*sh), inv=.TRUE.))
        
        end function fft_shift

        function erf(x)
        real(kind=8) :: erf, x, dumerfc, t, z

            z = abs(x)
            t = 1.0 / ( 1.0 + 0.5 * z )

            dumerfc =       t * exp(-z * z - 1.26551223 + t *       &
             ( 1.00002368 + t * ( 0.37409196 + t *      &
             ( 0.09678418 + t * (-0.18628806 + t *      &
             ( 0.27886807 + t * (-1.13520398 + t *      &
             ( 1.48851587 + t * (-0.82215223 + t * 0.17087277 )))))))))

            if ( x.lt.0.0 ) dumerfc = 2.0 - dumerfc
     
            erf = 1.0 - dumerfc

        end function erf
        
! ---------------------------------------------------------
! Solve a one-dimensional equality using the secant method
! ---------------------------------------------------------
    function confidenceLevel(x, nu, P)
    real(kind=8) :: confidenceLevel, x, nu, P
        confidenceLevel = gammq(0.5*nu,0.5*x) + P - 1.d0
        return 
    end function confidenceLevel

! ---------------------------------------------------------
! Solve a one-dimensional equality using the secant method
! ---------------------------------------------------------
    function secantConfidenceLevel(nu,P)
    real(kind=8) :: func, x1, x2, xacc, fl, f, secantConfidenceLevel, xl, swap, dx, nu,P
    integer :: j
    integer, parameter :: maxit=30
        xacc = 1e-3
        x1 = nu-0.1
        x2 = nu+14.d0
        fl=confidenceLevel(x1,nu,P)
        f=confidenceLevel(x2,nu,P)      
      if(abs(fl).lt.abs(f))then
        secantConfidenceLevel=x1
        xl=x2
        swap=fl
        fl=f
        f=swap
      else
        xl=x1
        secantConfidenceLevel=x2
      endif
      do j=1,maxit
        dx=(xl-secantConfidenceLevel)*f/(f-fl)
        xl=secantConfidenceLevel
        fl=f
        secantConfidenceLevel=secantConfidenceLevel+dx
        f=confidenceLevel(secantConfidenceLevel,nu,P)
        if(abs(dx).lt.xacc.or.f.eq.0.)return
      enddo
      print *, 'rtsec exceed maximum iterations'
      end function secantConfidenceLevel

! ---------------------------------------------------------
! Return the log of the gamma function
! ---------------------------------------------------------
      function gammln(xx)
      real(kind=8) :: cof(6),stp,half,one,fpf,x,tmp,ser,gammln, dx, xx
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,-1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data half,one,fpf/0.5d0,1.0d0,5.5d0/
      integer :: j
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*log(tmp)-tmp
      ser=one
      do j=1,6
        x=x+one
        ser=ser+cof(j)/x
      enddo
      gammln=tmp+log(stp*ser)
      return
      end function gammln

! ---------------------------------------------------------
! 
! ---------------------------------------------------------
     subroutine gcf(gammcf,a,x,gln)
     real(kind=8) :: gammcf, a, x, gln, gold, a0, a1, b0, b1, fac, an, ana, anf, g
     real(kind=8), parameter :: eps=3e-7
     integer, parameter :: itmax=200
     integer :: n
     
      gln=gammln(a)
      gold=0.
      a0=1.
      a1=x
      b0=0.
      b1=1.
      fac=1.
      do n=1,itmax
        an=float(n)
        ana=an-a
        a0=(a1+a0*ana)*fac
        b0=(b1+b0*ana)*fac
        anf=an*fac
        a1=x*a0+anf*a1
        b1=x*b0+anf*b1
        if(a1.ne.0.)then
          fac=1./a1
          g=b1*fac
          if(abs((g-gold)/g).lt.eps) then
            gammcf=exp(-x+a*log(x)-gln)*g
            return
          endif
          gold=g
        endif
      enddo
        print *, 'a too large, itmax too small'
      end subroutine gcf
      
      subroutine gser(gamser,a,x,gln)
      real(kind=8) :: gamser, a, x, gln, ap, sum, del
      real(kind=8), parameter :: eps=3e-7
      integer, parameter :: itmax=200
      integer :: n
          
      gln=gammln(a)
      if(x.le.0.)then        
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do n=1,itmax
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*eps) then
            gamser=sum*exp(-x+a*log(x)-gln)
            return
        endif
      enddo
      print *, 'a too large, itmax too small'
      end subroutine gser
      
    function gammq(a,x)
    real(kind=8) :: gammq, a, x, gln, gamser, gammcf
      if(x.lt.0..or.a.le.0.) then
        print *, x, a
        stop
    endif
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammq=1.-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      end function gammq


    
end module maths

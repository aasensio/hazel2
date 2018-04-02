c_______________________________________________________________
c StokesFR
c
c los perfiles de stokes estan ordenados en longitud de onda y 
c     luego parametro Stokes
c stok=i(l1,l2,...),q(l1,..),....v(l1....)
c
c las FR salen ordenadas en longitud de onda,perfil y tau
c rt(tau1:i(l1,l2,...),q(l1,..),....v(l1....);tau2:i......) 
c Las FR son a perturbaciones absolutas y no van convolucionadas
c con la macroturbulencia ni la PSF
c las FR estan evaluadas en todos los puntos en tau
c
c Basilio 26/02/2016
c_______________________________________________________________

        include 'PARAMETER'       !incluye kt,kn,kl,kld
        parameter (kld4=4*kld,kldt=kld*kn,kldt4=4*kldt) 

        real*4 stok(kld4)
        real*4 rt(kldt4),rp(kldt4),rh(kldt4),rv(kldt4)
        real*4 rg(kldt4),rf(kldt4),rmic(kldt4)
        real*4 rmac(kld4)
               
        call leyendo !carga en commons toda la informacion necesaria
              
        call StokesFRsub(stok,rt,rp,rh,rv,rg,rf,rmic,rmac)
                
        call escribe_stokes(stok)        

        end

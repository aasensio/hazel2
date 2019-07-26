c subroutine to read all the atomic parameters using leelineasii 
c and storing it in _all variables
  
        subroutine lee_all_lines                              
     
        include 'PARAMETER'         !includs kt,kn,kl,kld
       
        integer nxx,ixx,iln
        integer ntl,nlin(kl),npas(kl),nble(kl)
        real*4  dlamda(kld)
        character*2 atom,atom_all(kl)
        integer istage,istage_all(kl)
        real*4 wlengt,wlengt_all(kl)
        real*4 zeff,zeff_all(kl)
        real*4 energy, energy_all(kl)
        real*4 loggf, loggf_all(kl)
        integer mult(2), mult_all(2,kl)
        character*1 design(2), design_all(2,kl)
        real*4 tam(2), tam_all(2,kl)
        real*4 alfa, alfa_all(kl)
        real*4 sigma, sigma_all(kl)
        
        common/Malla/ntl,nlin,npas,nble,dlamda             !se carga en lee_malla.f
        common/Lineas_all/atom_all,istage_all,wlengt_all,zeff_all,  
     &           energy_all,loggf_all,mult_all,design_all,tam_all,  
     &           alfa_all,sigma_all
                       
c reading atomic parameters using leelineasii
        ixx=0
        do iln=1,ntl 
           do ible=1,nble(iln)
             ixx=ixx+1
             nxx=nlin(ixx) 
             if(nxx.eq.0)then
                nxx=nlin(ixx-1)
                call leelineasii(nxx,atom,istage,wlengt,zeff,energy,
     &                         loggf,mult,design,tam,alfa,sigma)
                loggf=-20.
                wlengt=5000.
             else
                call leelineasii(nxx,atom,istage,wlengt,zeff,energy,
     &                         loggf,mult,design,tam,alfa,sigma)
             end if 
             
             atom_all(ixx)=atom
             istage_all(ixx)=istage
             wlengt_all(ixx)=wlengt
             zeff_all(ixx)=zeff
             energy_all(ixx)=energy
             loggf_all(ixx)=loggf
             mult_all(1,ixx)=mult(1)
             mult_all(2,ixx)=mult(2)
             design_all(1,ixx)=design(1)
             design_all(2,ixx)=design(2)
             tam_all(1,ixx)=tam(1)
             tam_all(2,ixx)=tam(2)
             alfa_all(ixx)=alfa
             sigma_all(ixx)=sigma
             
           end do
        end do   
        
        return
        end
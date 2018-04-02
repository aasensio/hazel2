      FUNCTION SAHA(THETA,CHI,U1,U2,PE) 
C     SAHA-EGGERT EQUATION  
      SAHA=U2*10.**(9.0805126-THETA*CHI)/(U1*PE*THETA**2.5) 
c      print*,'saha en saha =          ',saha
c      print*,theta,chi,u1,u2,pe
      RETURN
      END   


      FUNCTION SAHAdouble(THETA,CHI,U1,U2,PE) 
      real*8 THETA,PE
      real*4 CHI,U1,U2

      SAHAdouble=U2*10.d0**(9.0805126d0-THETA*CHI)/(U1*PE*THETA**2.5d0) 

      RETURN
      END   

      FUNCTION SAHAdoubleinv(THETA,CHI,U1,U2,PE) 
      real*8 THETA,PE
      real*4 CHI,U1,U2

      SAHAdoubleinv=8.3078262e-10*10.d0**(THETA*CHI)*THETA**2.5d0*PE*U1/U2

      RETURN
      END   

      FUNCTION dSAHA(THETA,CHI,dU1,dU2)
c	calcula la detivada respecto a t, del logaritmo de la
C     SAHA-EGGERT EQUATION
c	la derivada del logaritmo de saha respecto a pe es -1/pe
      dSAHA=du2-du1+(theta/5040.)*(2.5+chi*THETA*alog(10.))
      RETURN
      END   

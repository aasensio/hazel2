!+-----------------------------------------------------------------------+
!| Program       : Direct.f (subfile DIRserial.f)                        |
!| Last modified : 04-12-2001                                            |
!| Written by    : Joerg Gablonsky                                       |
!| SUBROUTINEs, which differ depENDing on the serial or parallel version.|
!+-----------------------------------------------------------------------+

!+-----------------------------------------------------------------------+
!| SUBROUTINE for sampling.                                              |
!+-----------------------------------------------------------------------+
      SUBROUTINE DIRSamplef(c,ArrayI,delta,sample,new,length,&
                dwrit,logfile,f,free,maxI,point,fcn,x,l,fmin,&
                minpos,u,n,maxfunc,maxdeep,oops,fmax,&
                 IFeasiblef,IInfesiblef,&
                iidata, iisize, ddata, idsize, cdata, icsize) 
      IMPLICIT NONE

!+-----------------------------------------------------------------------+
!| JG 07/16/01 fcn must be declared external.                            |
!+-----------------------------------------------------------------------+
      EXTERNAL fcn

      INTEGER n,maxfunc,maxdeep,oops
      INTEGER maxI,ArrayI(n),sample,new
      INTEGER length(maxfunc,n),free,point(maxfunc),i
!+-----------------------------------------------------------------------+
!| JG 07/16/01 Removed fcn.                                              |
!+-----------------------------------------------------------------------+
      DOUBLE PRECISION c(maxfunc,n),delta,x(n)
      DOUBLE PRECISION l(n),u(n),f(maxfunc,2)
      DOUBLE PRECISION fmin
      INTEGER pos,j,dwrit,logfile,minpos
      INTEGER helppoint,kret
!+-----------------------------------------------------------------------+
!| JG 01/22/01 Added variable to keep track of the maximum value found.  |
!|             Added variable to keep track IF feasible point was found. |
!+-----------------------------------------------------------------------+
      DOUBLE PRECISION fmax
      INTEGER  IFeasiblef,IInfesiblef
!+-----------------------------------------------------------------------+
!| Variables to pass user defined data to the function to be optimized.  |
!+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      DOUBLE PRECISION ddata(idsize)
      CHARACTER*40 cdata(icsize)
      
!+-----------------------------------------------------------------------+
!| Set the pointer to the first function to be evaluated,                |
!| store this position also in helppoint.                                |
!+-----------------------------------------------------------------------+
      pos = new
      helppoint = pos
!+-----------------------------------------------------------------------+
!| Iterate over all points, where the function should be                 |
!| evaluated.                                                            |
!+-----------------------------------------------------------------------+
      DO 40,j=1,maxI + maxI
!+-----------------------------------------------------------------------+
!| Copy the position into the helparrayy x.                              |
!+-----------------------------------------------------------------------+
         DO 60,i=1,n
           x(i) = c(pos,i)
60       CONTINUE
!+-----------------------------------------------------------------------+
!| Call the function.                                                    |
!+-----------------------------------------------------------------------+
         CALL DIRinfcn(fcn,x,l,u,n,f(pos,1),kret,&
                      iidata, iisize, ddata, idsize, cdata, icsize)
!+-----------------------------------------------------------------------+
!| Remember IF an infeasible point has been found.                       |
!+-----------------------------------------------------------------------+
         IInfesiblef = max(IInfesiblef,kret)
         IF (kret .eq. 0) then
!+-----------------------------------------------------------------------+
!| IF the function evaluation was O.K., set the flag in                  |
!| f(pos,2). Also mark that a feasible point has been found.             |
!+-----------------------------------------------------------------------+
           f(pos,2) = 0.D0
            IFeasiblef = 0
!+-----------------------------------------------------------------------+
!| JG 01/22/01 Added variable to keep track of the maximum value found.  |
!+-----------------------------------------------------------------------+
           fmax = max(f(pos,1),fmax)
         END if
         IF (kret .ge. 1) then
!+-----------------------------------------------------------------------+
!|  IF the function could not be evaluated at the given point,            |
!| set flag to mark this (f(pos,2) and store the maximum                 |
!| box-sidelength in f(pos,1).                                           |
!+-----------------------------------------------------------------------+
           f(pos,2) = 2.D0
           f(pos,1) = fmax
         END if
!+-----------------------------------------------------------------------+
!|  IF the function could not be evaluated due to a failure in            |
!| the setup, mark this.                                                 |
!+-----------------------------------------------------------------------+
         IF (kret .eq. -1) then
           f(pos,2) = -1.D0
         END if
!+-----------------------------------------------------------------------+
!| Set the position to the next point, at which the function             |
!| should e evaluated.                                                   |
!+-----------------------------------------------------------------------+
         pos = point(pos)
40    CONTINUE
      pos = helppoint
!+-----------------------------------------------------------------------+
!| Iterate over all evaluated points and see, IF the minimal             |
!| value of the function has changed.  IF this has happEND,               |
!| store the minimal value and its position in the array.                |
!| Attention: Only valied values are checked!!                           |
!+-----------------------------------------------------------------------+
      DO 50,j=1,maxI + maxI
        IF ((f(pos,1) .LT. fmin) .and. (f(pos,2) .eq. 0)) THEN
          fmin = f(pos,1) 
          minpos = pos
        END IF
        pos = point(pos)
50    CONTINUE
      END

!+-----------------------------------------------------------------------+
!| Problem-specific Initialisation                                       |
!+-----------------------------------------------------------------------+
      SUBROUTINE DIRInitSpecific(x,n, iidata, iisize, ddata, idsize, cdata, icsize)
      IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION x(n)
!+-----------------------------------------------------------------------+
!| Variables to pass user defined data to the function to be optimized.  |
!+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      DOUBLE PRECISION ddata(idsize)
      CHARACTER*40 cdata(icsize)
!+-----------------------------------------------------------------------+
!| Problem - specific variables !                                        |
!+-----------------------------------------------------------------------+

!+-----------------------------------------------------------------------+
!| END of problem - specific variables !                                 |
!+-----------------------------------------------------------------------+

!+-----------------------------------------------------------------------+
!| Start of problem-specific initialisation                              |
!+-----------------------------------------------------------------------+

!+-----------------------------------------------------------------------+
!| END of problem-specific initialisation                                |
!+-----------------------------------------------------------------------+
      END

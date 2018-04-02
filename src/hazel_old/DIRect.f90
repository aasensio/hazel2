!+-----------------------------------------------------------------------+
!| Program       : Direct.f                                              |
!| Last modified : 07-16-2001                                            |
!| Written by    : Joerg Gablonsky (jmgablon@unity.ncsu.edu)             |
!|                 North Carolina State University                       |
!|                 Dept. of Mathematics                                  |
!| DIRECT is a method to solve problems of the form:                     |
!|              min f: Q --> R,                                          |
!| where f is the function to be minimized and Q is an n-dimensional     |
!| hyperrectangle given by the the following equation:                   |
!|                                                                       |
!|       Q={ x : l(i) <= x(i) <= u(i), i = 1,...,n }.                    |
!| Note: This version of DIRECT can also handle hidden constraints. By   |
!|       this we mean that the function may not be defined over the whole|
!|       hyperrectangle Q, but only over a subset, which is not given    |
!|       analytically.                                                   |
!|                                                                       |
!| We now give a brief outline of the algorithm:                         |
!|                                                                       |
!|   The algorithm starts with mapping the hyperrectangle Q to the       |
!|   n-dimensional unit hypercube. DIRECT then samples the function at   |
!|   the center of this hypercube and at 2n more points, 2 in each       |
!|   coordinate direction. Uisng these function values, DIRECT then      |
!|   divides the domain into hyperrectangles, each having exactly one of |
!|   the sampling points as its center. In each iteration, DIRECT chooses|
!|   some of the existing hyperrectangles to be further divided.         |
!|   We provide two different strategies of how to decide which          |
!|   hyperrectangles DIRECT divides and several different convergence    |
!|   criteria.                                                           |
!|                                                                       |
!|   DIRECT was designed to solve problems where the function f is       |
!|   Lipschitz continues. However, DIRECT has proven to be effective on  |
!|   more complex problems than these.                                   |
!+-----------------------------------------------------------------------+

      SUBROUTINE Direct(fcn, x, n, eps, maxf, maxT, fmin, l, u,&
                       algmethod, Ierror, logfile, &
                       fglobal, fglper, volper, sigmaper,&
                       iidata, iisize, ddata, idsize, cdata, icsize,&
                       verbose)

!+-----------------------------------------------------------------------+
!|    SUBROUTINE Direct                                                  |
!| On entry                                                              |
!|     fcn -- The argument containing the name of the user-supplied      |
!|            SUBROUTINE that returns values for the function to be      |
!|            minimized.                                                 |
!|       n -- The dimension of the problem.                              |
!|     eps -- Exceeding value. If eps > 0, we use the same epsilon for   |
!|            all iterations. If eps < 0, we use the update formula from |
!|            Jones:                                                     |
!|               eps = max(1.D-4*abs(fmin),epsfix),                      |
!|            where epsfix = abs(eps), the absolute value of eps which is|
!|            passed to the function.                                    |
!|    maxf -- The maximum number of function evaluations.                |
!|    maxT -- The maximum number of iterations.                          |
!|            Direct stops when either the maximum number of iterations  |
!|            is reached or more than maxf function-evalutions were made.|
!|       l -- The lower bounds of the hyperbox.                          |
!|       u -- The upper bounds of the hyperbox.                          |
!|algmethod-- Choose the method, that is either use the original method  |
!|            as described by Jones et.al. (0) or use our modification(1)|
!| logfile -- File-Handle for the logfile. DIRECT expects this file to be|
!|            opened and closed by the user outside of DIRECT. We moved  |
!|            this to the outside so the user can add extra informations |
!|            to this file before and after the call to DIRECT.          |
!| fglobal -- Function value of the global optimum. If this value is not |
!|            known (that is, we solve a real problem, not a testproblem)|
!|            set this value to -1.D100 and fglper (see below) to 0.D0.  |
!|  fglper -- Terminate the optimization when the percent error          |
!|                100(f_min - fglobal)/max(1,abs(fglobal)) < fglper.     |
!|  volper -- Terminate the optimization when the volume of the          |
!|            hyperrectangle S with f(c(S)) = fmin is less then volper   |
!|            percent of the volume of the original hyperrectangle.      |
!|sigmaper -- Terminate the optimization when the measure of the         |
!|            hyperrectangle S with f(c(S)) = fmin is less then sigmaper.|
!|                                                                       |
!| User data that is passed through without being changed:               |
!|  iidata -- Integer Array of user data. This array of size iisize is   |
!|            passed to the function to be optimized and can be used to  |
!|            transfer data to this function. The contents are not       |
!|            changed by DIRECT.                                         |
!|  iisize -- Size of array iidata.                                      |
!|   ddata -- DOUBLE PRECISION array of user data. See iidata.           |
!|  idsize -- Size of array ddata.                                       |
!|   cdata -- Character array. See iidata.                               |
!|  icsize -- Size of array ddata.                                       |
!|                                                                       |
!| On return                                                             |
!|                                                                       |
!|       x -- The final point obtained in the optimization process.      |
!|            X should be a good approximation to the global minimum     |
!|            for the function within the hyper-box.                     |
!|                                                                       |
!|    fmin -- The value of the function at x.                            |
!|  Ierror -- Error flag. If Ierror is lower 0, an error has occured. The|
!|            values of Ierror mean                                      |
!|            Fatal errors :                                             |
!|             -1   u(i) <= l(i) for some i.                             |
!|             -2   maxf is too large.                                   |
!|             -3   Initialization in DIRpreprc failed.                  |
!|             -4   Error in DIRSamplepoints, that is there was an error |
!|                  in the creation of the sample points.                |
!|             -5   Error in DIRSamplef, that is an error occured while  |
!|                  the function was sampled.                            |
!|             -6   Error in DIRDoubleInsert, that is an error occured   |
!|                  DIRECT tried to add all hyperrectangles with the same|
!|                  size and function value at the center. Either        |
!|                  increase maxdiv or use our modification (Jones = 1). |
!|            Termination values :                                       |
!|              1   Number of function evaluations done is larger then   |
!|                  maxf.                                                |
!|              2   Number of iterations is equal to maxT.               |
!|              3   The best function value found is within fglper of    |
!|                  the (known) global optimum, that is                  |
!|                   100(fmin - fglobal/max(1,|fglobal|))  < fglper.     |
!|                  Note that this termination signal only occurs when   |
!|                  the global optimal value is known, that is, a test   |
!|                  function is optimized.                               |
!|              4   The volume of the hyperrectangle with fmin at its    |
!|                  center is less than volper percent of the volume of  |
!|                  the original hyperrectangle.                         |
!|              5   The measure of the hyperrectangle with fmin at its   |
!|                  center is less than sigmaper.                        |
!|                                                                       |
!| SUBROUTINEs used :                                                    |
!|                                                                       |
!| DIRheader, DIRInitSpecific, DIRInitList, DIRpreprc, DIRInit, DIRChoose|
!| DIRDoubleInsert, DIRGet_I, DIRSamplepoints, DIRSamplef, DIRDivide     |
!| DIRInsertList, DIRreplaceInf, DIRWritehistbox, DIRsummary, Findareas  |
!|                                                                       |
!| Functions used :                                                      |
!|                                                                       |
!| DIRgetMaxdeep, DIRgetlevel                                            |
!+-----------------------------------------------------------------------+

      IMPLICIT NONE
!+-----------------------------------------------------------------------+
!| Parameters                                                            |
!+-----------------------------------------------------------------------+

!+-----------------------------------------------------------------------+
!| The maximum of function evaluations allowed.                          |
!| The maximum dept of the algorithm.                                    |
!| The maximum number of divisions allowed.                              |
!| The maximal dimension of the problem.                                 |
!+-----------------------------------------------------------------------+
      INTEGER maxfunc, maxdeep, maxdiv, MaxDim, mdeep
      PARAMETER (Maxfunc = 2000)
      PARAMETER (maxdeep = 600)
      PARAMETER (maxdiv = 3000)
      PARAMETER (MaxDim = 64)

!+-----------------------------------------------------------------------+
!| Global Variables.                                                     |
!+-----------------------------------------------------------------------+
      INTEGER JONES
      COMMON /directcontrol/ JONES

!+-----------------------------------------------------------------------+
!| EXTERNAL Variables.                                                   |
!+-----------------------------------------------------------------------+
      EXTERNAL fcn
      INTEGER n, maxf, maxT, algmethod, Ierror, logfile, dwrit, verbose
      DOUBLE PRECISION  x(n),fmin,eps,l(n),u(n)
      DOUBLE PRECISION fglobal, fglper, volper, sigmaper

!+-----------------------------------------------------------------------+
!| User Variables.                                                       |
!| These can be used to pass user defined data to the function to be     |
!| optimized.                                                            |
!+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      DOUBLE PRECISION ddata(idsize)
      Character*40 cdata(icsize)

!+-----------------------------------------------------------------------+
!| Place to define, if needed, some application-specific variables.      |
!| Note: You should try to use the arrays defined above for this.        |
!+-----------------------------------------------------------------------+

!+-----------------------------------------------------------------------+
!| End of application - specific variables !                             |
!+-----------------------------------------------------------------------+


!+-----------------------------------------------------------------------+
!| Internal variables :                                                  |
!|       f -- values of functions.                                       |
!|divfactor-- Factor used for termination with known global minimum.     |
!|  anchor -- anchors of lists with deepness i, -1 is anchor for list of |
!|            NaN - values.                                              |
!|       S -- List of potentially optimal points.                        |
!|   point -- lists                                                      |
!|    free -- first free position                                        |
!|       c -- midpoints of arrays                                        |
!|  thirds -- Precalculated values of 1/3^i.                             |
!|  levels -- Length of intervals.                                       |
!|  length -- Length of intervall (index)                                |
!|       t -- actual iteration                                           |
!|       j -- loop-variable                                              |
!| actdeep -- the actual minimal interval-length index                   |
!|  Minpos -- position of the actual minimum                             |
!|    file -- The filehandle for a datafile.                             |
!|  maxpos -- The number of intervalls, which are truncated.             |
!|    help -- A help variable.                                           |
!| numfunc -- The actual number of function evaluations.                 |
!|   file2 -- The filehandle for an other datafile.                      |
!|  ArrayI -- Array with the indexes of the sides with maximum length.   |
!|    maxi -- Number of directions with maximal side length.             |
!|    oops -- Flag which shows if anything went wrong in the             |
!|            initialisation.                                            |
!|   cheat -- Obsolete. If equal 1, we don't allow Ktilde > kmax.        |
!|  writed -- If writed=1, store final division to plot with Matlab.     |
!|   List2 -- List of indicies of intervalls, which are to be truncated. |
!|       i -- Another loop-variable.                                     |
!|actmaxdeep-- The actual maximum (minimum) of possible Interval length. |
!|  oldpos -- The old index of the minimum. Used to print only, if there |
!|            is a new minimum found.                                    |
!|  tstart -- The start of the outer loop.                               |
!|   start -- The postion of the starting point in the inner loop.       |
!| Newtosample -- The total number of points to sample in the inner loop.|
!|       w -- Array used to divide the intervalls                        |
!|    kmax -- Obsolete. If cheat = 1, Ktilde was not allowed to be larger|
!|            than kmax. If Ktilde > kmax, we set ktilde = kmax.         |
!|   delta -- The distance to new points from center of old hyperrec.    |
!|    pos1 -- Help variable used as an index.                            |
!| version -- Store the version number of DIRECT.                        |
!| oldmaxf -- Store the original function budget.                        |
!|increase -- Flag used to keep track if function budget was increased   |
!|            because no feasible point was found.                       |
!| freeold -- Keep track which index was free before. Used with          |
!|            SUBROUTINE DIRReplaceInf.                                  |
!|actdeep_div-- Keep track of the current depths for divisions.          |
!|    oldl -- Array used to store the original bounds of the domain.     |
!|    oldu -- Array used to store the original bounds of the domain.     |
!|  epsfix -- If eps < 0, we use Jones update formula. epsfix stores the |
!|            absolute value of epsilon.                                 |
!|iepschange-- flag iepschange to store if epsilon stays fixed or is     |
!|             changed.                                                  |
!|DIRgetMaxdeep-- Function to calculate the level of a hyperrectangle.   |
!|DIRgetlevel-- Function to calculate the level and stage of a hyperrec. |
!|    fmax -- Keep track of the maximum value of the function found.     |
!|Ifeasiblef-- Keep track if a feasible point has  been found so far.    |
!|             Ifeasiblef = 0 means a feasible point has been found,     |
!|             Ifeasiblef = 1 no feasible point has been found.          |
!|   dwrit -- Controls the level of output. So far not used a lot, set to|
!|            0. If its value is set to 2, more output, exspecially from |
!|            Subroutine DIRChoose.                                      |
!+-----------------------------------------------------------------------+
      DOUBLE PRECISION  f(maxfunc,2), divfactor
      INTEGER anchor(-1:maxdeep),S(maxdiv,2)
      INTEGER point(maxfunc), free
      DOUBLE PRECISION  c(maxfunc,MaxDim)
      DOUBLE PRECISION  thirds(0:maxdeep),levels(0:maxdeep)
      INTEGER length(maxfunc,MaxDim),t,j,actdeep
      INTEGER Minpos,file,maxpos,help,numfunc,file2
      INTEGER ArrayI(MaxDim),maxi,oops,cheat,writed
      INTEGER List2(MaxDim,2),i,actmaxdeep,oldpos
      INTEGER tstart,start,Newtosample
      DOUBLE PRECISION  w(MaxDim),kmax, delta
      INTEGER pos1
!+-----------------------------------------------------------------------+
!| JG 09/25/00 Version counter.                                          |
!+-----------------------------------------------------------------------+
      INTEGER version
      INTEGER oldmaxf,increase, freeold
!+-----------------------------------------------------------------------+
!| JG 09/24/00 Add another actdeep to keep track of the current depths   |
!|             for divisions.                                            |
!+-----------------------------------------------------------------------+
      INTEGER actdeep_div
      DOUBLE PRECISION  oldl(MaxDim), oldu(MaxDim)
!+-----------------------------------------------------------------------+
!|JG 01/13/01 Added epsfix for epsilon update. If eps < 0, we use Jones  |
!|            update formula. epsfix stores the absolute value of epsilon|
!|            then. Also added flag iepschange to store if epsilon stays |
!|            fixed or is changed.                                       |
!+-----------------------------------------------------------------------+
      DOUBLE PRECISION epsfix
      INTEGER iepschange

      INTEGER DIRGetMaxdeep, DIRgetlevel
!+-----------------------------------------------------------------------+
!| JG 01/22/01 fmax is used to keep track of the maximum value found.    |
!+-----------------------------------------------------------------------+
      DOUBLE PRECISION fmax
!+-----------------------------------------------------------------------+
!| JG 01/22/01 Ifeasiblef is used to keep track if a feasible point has  |
!|             been found so far. Ifeasiblef = 0 means a feasible point  |
!|             has been found, Ifeasiblef = 1 if not.                    |
!| JG 03/09/01 IInfeasible is used to keep track if an infeasible point  |
!|             has been found. IInfeasible > 0 means a infeasible point  |
!|             has been found, IInfeasible = 0 if not.                   |
!+-----------------------------------------------------------------------+
      INTEGER Ifeasiblef, IInfesiblef

!+-----------------------------------------------------------------------+
!+-----------------------------------------------------------------------+
!|                            Start of code.                             |
!+-----------------------------------------------------------------------+
!+-----------------------------------------------------------------------+
      writed = 0
      dwrit  = 0
      JONES  = algmethod
!+-----------------------------------------------------------------------+
!| Save the upper and lower bounds.                                      |
!+-----------------------------------------------------------------------+
      DO 150,i=1,n
        oldu(i) = u(i)
        oldl(i) = l(i)
150   CONTINUE

!+-----------------------------------------------------------------------+
!| Set version.                                                          |
!+-----------------------------------------------------------------------+
      version = 204
!+-----------------------------------------------------------------------+
!| Set parameters.                                                       |
!|    If cheat > 0, we do not allow \tilde{K} to be larger than kmax, and|
!|    set \tilde{K} to set value if necessary. Not used anymore.         |
!+-----------------------------------------------------------------------+
      cheat = 0
      kmax  = 1.D10
      mdeep = maxdeep
!+-----------------------------------------------------------------------+
!| Write the header of the logfile.                                      |
!+-----------------------------------------------------------------------+
      CALL DIRheader(logfile, version, x, n, eps, maxf, maxT, l, u,&
                    algmethod, maxfunc, maxdeep, fglobal, fglper,&
                    Ierror, epsfix, iepschange, volper, sigmaper,&
                    iidata, iisize, ddata, idsize, cdata,&
                    icsize, verbose)
!+-----------------------------------------------------------------------+
!| If an error has occured while writing the header (we do some checking |
!| of variables there), return to the main program.                      |
!+-----------------------------------------------------------------------+
      IF (Ierror .lt. 0) then
         RETURN
      END IF
!+-----------------------------------------------------------------------+
!| If the known global minimum is equal 0, we cannot divide by it.       |
!| Therefore we set it to 1. If not, we set the divisionfactor to the    |
!| absolute value of the global minimum.                                 |
!+-----------------------------------------------------------------------+
      IF (fglobal .EQ. 0.D0) then
         divfactor = 1.D0
      ELSE
         divfactor = abs(fglobal)
      END IF

!+-----------------------------------------------------------------------+
!| Start of application-specific initialisation.                         |
!+-----------------------------------------------------------------------+
      CALL DIRInitSpecific(x,n,&
        iidata, iisize, ddata, idsize, cdata, icsize)
!+-----------------------------------------------------------------------+
!| End of application-specific initialisation.                           |
!+-----------------------------------------------------------------------+

!+-----------------------------------------------------------------------+
!| Save the budget given by the user. The variable maxf will be changed  |
!| if in the beginning no feasible points are found.                     |
!+-----------------------------------------------------------------------+
      oldmaxf = maxf
      increase = 0
!+-----------------------------------------------------------------------+
!| Initialiase the lists.                                                |
!+-----------------------------------------------------------------------+
      CALL DIRInitList(anchor,free,point,f,maxfunc,maxdeep)
!+-----------------------------------------------------------------------+
!| Call the routine to initialise the mapping of x from the n-dimensional|
!| unit cube to the hypercube given by u and l. If an error occured,     |
!| give out a error message and return to the main program with the error|
!| flag set.                                                             |
!| JG 07/16/01 Changed call to remove unused data.                       |
!+-----------------------------------------------------------------------+
      CALL DIRpreprc(u,l,n,l,u,oops)
      IF (oops .GT. 0) THEN
        if (verbose == 1) Write(*,10005)
        if (verbose == 1) Write(logfile,10005)
        IError = -3
        Return
      END IF
      tstart = 2
!+-----------------------------------------------------------------------+
!| Initialise the algorithm DIRECT.                                      |
!+-----------------------------------------------------------------------+
!+-----------------------------------------------------------------------+
!| Added variable to keep track of the maximum value found.              |
!+-----------------------------------------------------------------------+
      CALL DIRInit(f,fcn,c,length,actdeep,point,anchor,free,&
        dwrit,logfile,ArrayI,maxI,List2,w,x,l,u,fmin,minpos,&
        thirds,levels,maxfunc,maxdeep,n,MaxDim,fmax,Ifeasiblef,&
        IInfesiblef, Ierror,&
        iidata, iisize, ddata, idsize, cdata, icsize)
!+-----------------------------------------------------------------------+
!| Added error checking.                                                 |
!+-----------------------------------------------------------------------+
      IF (Ierror .lt. 0) then
         IF (Ierror .eq. -4) THEN
            if (verbose == 1) Write(*,10006)
            if (verbose == 1) Write(logfile,10006)
            return
         END IF
         IF (Ierror .eq. -5) THEN
            if (verbose == 1) Write(*,10007)
            if (verbose == 1) Write(logfile,10007)
            return
         END IF
      END IF
      
      numfunc = 1 + maxI + maxI
      actmaxdeep = 1
      oldpos = 0
      tstart = 2
!+-----------------------------------------------------------------------+
!| If no feasible point has been found, give out the iteration, the      |
!| number of function evaluations and a warning. Otherwise, give out     |
!| the iteration, the number of function evaluations done and fmin.      |
!+-----------------------------------------------------------------------+
      IF (Ifeasiblef .gt. 0) then
        if (verbose == 1) write(*,10012) tstart-1,numfunc
        if (verbose == 1) write(logfile,10012) t,numfunc
      ELSE
        if (verbose == 1) Write(*,10002) numfunc, fmin, fmax
        if (verbose == 1) Write(logfile,10003) tstart-1,numfunc,fmin
      END IF
!+-----------------------------------------------------------------------+
!+-----------------------------------------------------------------------+
!| Main loop!                                                            |
!+-----------------------------------------------------------------------+
!+-----------------------------------------------------------------------+
      DO 10,t=tstart,MaxT
!+-----------------------------------------------------------------------+
!| Choose the sample points. The indices of the sample points are stored |
!| in the list S.                                                        |
!+-----------------------------------------------------------------------+
        actdeep = actmaxdeep
        CALL DIRChoose(anchor,S,maxdeep,f,fmin,eps,levels,maxpos,&
          length,maxfunc,maxdeep,maxdiv,n,logfile,dwrit,cheat,kmax,&
          Ifeasiblef)
!+-----------------------------------------------------------------------+
!| Add other hyperrectangles to S, which have the same level and the same|
!| function value at the center as the ones found above (that are stored |
!| in S). This is only done if we use the original DIRECT algorithm.     |
!| JG 07/16/01 Added Errorflag.                                          |
!+-----------------------------------------------------------------------+
        IF (algmethod .EQ. 0) THEN
          CALL DIRDoubleInsert(anchor, S, maxpos, point, f,&
                maxdeep, maxfunc, maxdiv, Ierror)
          IF (Ierror .eq. -6) THEN
              if (verbose == 1) Write(*,10020)
              if (verbose == 1) Write(*,10021)
              if (verbose == 1) Write(*,10022)
              if (verbose == 1) Write(*,10023)
              if (verbose == 1) Write(logfile,10020)
              if (verbose == 1) Write(logfile,10021)
              if (verbose == 1) Write(logfile,10022)
              if (verbose == 1) Write(logfile,10023)
              return
          END IF
        ENDIF
        oldpos = minpos
!+-----------------------------------------------------------------------+
!| Initialise the number of sample points in this outer loop.            |
!+-----------------------------------------------------------------------+
        Newtosample = 0
        DO 20,j=1,maxpos
           actdeep = S(j,2)
!+-----------------------------------------------------------------------+
!| If the actual index is a point to sample, do it.                      |
!+-----------------------------------------------------------------------+
           IF (S(j,1) .GT. 0) THEN
!+-----------------------------------------------------------------------+
!| JG 09/24/00 Calculate the value delta used for sampling points.       |
!+-----------------------------------------------------------------------+
              actdeep_div = DIRGetmaxdeep(S(j,1),length,maxfunc,n)
              delta = thirds(actdeep_div+1)
              actdeep = S(j,2)
!+-----------------------------------------------------------------------+
!| If the current dept of division is only one under the maximal allowed |
!| dept, stop the computation.                                           |
!+-----------------------------------------------------------------------+
              IF (actdeep+1 .GE. mdeep) THEN
                 if (verbose == 1) Write(*,10004)
                 if (verbose == 1) write(logfile,10004)
                 Ierror = -6
                 GOTO 100
              END IF
              actmaxdeep = max(actdeep,actmaxdeep)
              help = S(j,1)
              IF (.NOT. (anchor(actdeep) .EQ. help)) THEN
                pos1 = anchor(actdeep)
                DO WHILE (.NOT. (point(pos1) .EQ. help))
                  pos1 = point(pos1)
                END DO
                point(pos1) = point(help)
              ELSE
                anchor(actdeep) = point(help)
              END IF
              IF (actdeep .lt. 0) then
                actdeep = f(help,1)
              END IF
!+-----------------------------------------------------------------------+
!| Get the Directions in which to decrease the intervall-length.         |
!+-----------------------------------------------------------------------+
              CALL DIRGet_I(length,help,ArrayI,maxI,n,maxfunc)
!+-----------------------------------------------------------------------+
!| Sample the function. To do this, we first calculate the points where  |
!| we need to sample the function. After checking for errors, we then do |
!| the actual evaluation of the function, again followed by checking for |
!| errors.                                                               |
!+-----------------------------------------------------------------------+
              CALL DIRSamplepoints(c,ArrayI,delta,help,&
                  start,length,dwrit,logfile,f,free,maxI,point,&
                  fcn,x,l,fmin,minpos,u,n,maxfunc,maxdeep,oops)
              IF (oops .GT. 0) THEN
                if (verbose == 1) Write(*,10006)
                if (verbose == 1) Write(logfile,10006)
                IError = -4
                return
              END IF
              Newtosample = newtosample + maxI
!+-----------------------------------------------------------------------+
!| JG 01/22/01 Added variable to keep track of the maximum value found.  |
!+-----------------------------------------------------------------------+
              CALL DIRSamplef(c,ArrayI,delta,help,start,&
                 length,dwrit,logfile,f,free,maxI,point,fcn,x,l,&
                 fmin,minpos,u,n,maxfunc,maxdeep,oops,fmax,&
                 Ifeasiblef,IInfesiblef,&
                 iidata, iisize, ddata, idsize, cdata, icsize)
              IF (oops .GT. 0) THEN
                if (verbose == 1) Write(*,10007)
                if (verbose == 1) Write(logfile,10007)
                IError = -5
                return
              END IF
!+-----------------------------------------------------------------------+
!| Divide the intervalls.                                                |
!+-----------------------------------------------------------------------+
              CALL DIRDivide(start,actdeep_div,length,point,&
                  ArrayI,help,List2,w,maxI,f,maxfunc,maxdeep,n)
!+-----------------------------------------------------------------------+
!| Insert the new intervalls into the list (sorted).                     |
!+-----------------------------------------------------------------------+
              CALL DIRInsertList(start,anchor,point,f,maxI,length,&
                         maxfunc,maxdeep,n,help)
!+-----------------------------------------------------------------------+
!| Increase the number of function evaluations.                          |
!+-----------------------------------------------------------------------+
              numfunc = numfunc + maxI + maxI
           END IF
!+-----------------------------------------------------------------------+
!| End of main loop.                                                     |
!+-----------------------------------------------------------------------+
20      CONTINUE
!+-----------------------------------------------------------------------+
!| If there is a new minimum, show the actual iteration, the number of   |
!| function evaluations, the minimum value of f (so far) and the position|
!| in the array.                                                         |
!+-----------------------------------------------------------------------+
        IF (oldpos .LT. minpos) THEN
          if (verbose == 1) Write(*,10002) numfunc,fmin, fmax
          if (verbose == 1) Write(logfile,10003) t,numfunc,fmin
        END IF
!+-----------------------------------------------------------------------+
!| If no feasible point has been found, give out the iteration, the      |
!| number of function evaluations and a warning.                         |
!+-----------------------------------------------------------------------+
        IF (Ifeasiblef .gt. 0) then
          if (verbose == 1) write(*,10012) t,numfunc
          if (verbose == 1) write(logfile,10012) t,numfunc
        END IF
!+-----------------------------------------------------------------------+
!+-----------------------------------------------------------------------+
!|                       Termination Checks                              |
!+-----------------------------------------------------------------------+

!+-----------------------------------------------------------------------+
!| JG 01/22/01 Calculate the index for the hyperrectangle at which       |
!|             fmin is assumed. We then calculate the volume of this     |
!|             hyperrectangle and store it in delta. This delta can be   |
!|             used to stop DIRECT once the volume is below a certain    |
!|             percentage of the original volume. Since the original     |
!|             is 1 (scaled), we can stop once delta is below a certain  |
!|             percentage, given by volper.                              |
!+-----------------------------------------------------------------------+
        Ierror = Jones
        Jones = 0
        actdeep_div = DIRGetlevel(minpos,length,maxfunc,n)
        Jones = Ierror
!+-----------------------------------------------------------------------+
!| JG 07/16/01 Use precalculated values to calculate volume.             |
!+-----------------------------------------------------------------------+
        delta = thirds(actdeep_div)*100
        IF (delta .LE. volper) THEN
           Ierror = 4
           if (verbose == 1) Write(*,10011) delta, volper
           if (verbose == 1) Write(logfile,10011) delta, volper
           GOTO 100
        END IF
!+-----------------------------------------------------------------------+
!| JG 01/23/01 Calculate the measure for the hyperrectangle at which     |
!|             fmin is assumed. If this measure is smaller then sigmaper,|
!|             we stop DIRECT.                                           |
!+-----------------------------------------------------------------------+
        actdeep_div = DIRGetlevel(minpos,length,maxfunc,n)
        delta = levels(actdeep_div)
        IF (delta .LE. sigmaper) THEN
           Ierror = 5
           if (verbose == 1) Write(*,10013) delta, sigmaper
           if (verbose == 1) Write(logfile,10013) delta, sigmaper
           GOTO 100
        END IF
!+-----------------------------------------------------------------------+
!| If the best found function value is within fglper of the (known)      |
!| global minimum value, terminate. This only makes sense if this optimal|
!| value is known, that is, in test problems.                            |
!+-----------------------------------------------------------------------+
        IF ((100*(fmin - fglobal)/divfactor) .LE. fglper) THEN
           Ierror = 3
           if (verbose == 1) Write(*,10010)
           if (verbose == 1) Write(logfile,10010)
           GOTO 100
        END IF
!+-----------------------------------------------------------------------+
!| Find out if there are infeasible points which are near feasible ones. |
!| If this is the case, replace the function value at the center of the  |
!| hyper rectangle by the lowest function value of a nearby function.    |
!| If no infeasible points exist (IInfesiblef = 0), skip this.           |
!+-----------------------------------------------------------------------+
        IF (IInfesiblef .gt. 0) THEN
           CALL DIRreplaceInf(free,freeold,f,c,thirds,length,anchor,&
            point,u,l,maxfunc,maxdeep,maxdim,n,logfile, fmax)
        ENDIF
        freeold = free
!+-----------------------------------------------------------------------+
!| If iepschange = 1, we use the epsilon change formula from Jones.      |
!+-----------------------------------------------------------------------+
        IF (iepschange .eq. 1) then
           eps = max(1.D-4*abs(fmin),epsfix)
        END IF
!+-----------------------------------------------------------------------+
!| If no feasible point has been found yet, set the maximum number of    |
!| function evaluations to the number of evaluations already done plus   |
!| the budget given by the user.                                         |
!| If the budget has already be increased, increase it again. If a       |
!| feasible point has been found, remark that and reset flag. No further |
!| increase is needed.                                                   |
!+-----------------------------------------------------------------------+
        IF (increase .eq. 1) then        
           maxf = numfunc + oldmaxf
           IF (Ifeasiblef .eq. 0) then
           if (verbose == 1) write(logfile,10031) maxf
             increase = 0
           END IF
        END IF
!+-----------------------------------------------------------------------+
!| Check if the number of function evaluations done is larger than the   |
!| allocated budget. If this is the case, check if a feasible point was  |
!| found. If this is a case, terminate. If no feasible point was found,  |
!| increase the budget and set flag increase.                            |
!+-----------------------------------------------------------------------+
        IF (numfunc .GT. maxf) THEN
           IF (Ifeasiblef .eq. 0) then
              Ierror = 1
              if (verbose == 1) Write(*,10008)
              if (verbose == 1) Write(logfile,10008)
              GOTO 100
           ELSE
              increase = 1
              if (verbose == 1) write(logfile,10030) numfunc
              maxf = numfunc+ oldmaxf
           END IF
        END IF
10    CONTINUE
!+-----------------------------------------------------------------------+
!+-----------------------------------------------------------------------+
!| End of main loop.                                                     |
!+-----------------------------------------------------------------------+
!+-----------------------------------------------------------------------+

!+-----------------------------------------------------------------------+
!| The algorithm stopped after maxT iterations.                          |
!+-----------------------------------------------------------------------+
      Ierror = 2
      if (verbose == 1) Write(*,10009)
      if (verbose == 1) Write(logfile,10009)

100   CONTINUE
!+-----------------------------------------------------------------------+
!| Store the position of the minimum in x.                               |
!+-----------------------------------------------------------------------+
      DO 50,i=1,n
         x(i) = c(Minpos,i)*l(i)+l(i)*u(i)
         u(i) = oldu(i)
         l(i) = oldl(i)
50    CONTINUE
!+-----------------------------------------------------------------------+
!| Store the number of function evaluations in maxf.                     |
!+-----------------------------------------------------------------------+
      maxf = numfunc

!+-----------------------------------------------------------------------+
!| If needed, save the final division in a file for use with Matlab.     |
!+-----------------------------------------------------------------------+
      writed = 0
      IF (writed .EQ. 1) THEN
        file  = 12
        file2 =  0
        CALL DIRWritehistbox(point,f,thirds,c,anchor,actdeep,file,&
                    l,u,file2,maxfunc,maxdeep,n,MaxDim,length)
      END IF

!+-----------------------------------------------------------------------+
!| Give out a summary of the run.                                        |
!+-----------------------------------------------------------------------+
      CALL DIRsummary(logfile,x,l,u,n,fmin,fglobal, numfunc, Ierror,&
      verbose)

!+-----------------------------------------------------------------------+
!| Format statements.                                                    |
!+-----------------------------------------------------------------------+
10002  FORMAT(i5," & ",f18.10," & ",f18.10," \\\\ ")
10003  FORMAT(i5,'       ', i5,'     ',f18.10)
10004  FORMAT('WARNING : Maximum number of levels reached. Increase maxdeep.')
10005  FORMAT('WARNING : Initialisation in DIRpreprc failed.')
10006  FORMAT('WARNING : Error occured in routine DIRsamplepoints.')
10007  FORMAT('WARNING : Error occured in routine DIRsamplef.')
10008  FORMAT('DIRECT stopped: numfunc >= maxf.')
10009  FORMAT('DIRECT stopped: maxT iterations.')
10010  FORMAT('DIRECT stopped: fmin within fglper of global minimum.')
10011  FORMAT('DIRECT stopped: Volume of S_min is ',d9.2, '% < ',d9.2,'% of the original volume.')
10012  FORMAT('No feasible point found in ',I4,' iterations ', 'and ',I5,' function evaluations.')
10013  FORMAT('DIRECT stopped: Measure of S_min = ',d9.2,' < ' ,d9.2,'.')
10020  FORMAT('WARNING : Capacity of array S in DIRDoubleInsert reached. Increase maxdiv.')
10021  FORMAT('This means that there are a lot of hyperrectangles')
10022  FORMAT('with the same function value at the center. We')
10023  FORMAT('suggest to use our modification instead (Jones = 1)')

10030 FORMAT('DIRECT could not find a feasible ', 'point after ', I5, ' function evaluations. ', 'DIRECT continues until a feasible point is found.')
10031 FORMAT('DIRECT found a feasible point. ', 'The adjusted budget is now set to ',I5,'.')
      END
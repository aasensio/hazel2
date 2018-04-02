      function planck(t,x)
      implicit real*4 (d)
      a=1.43880/x
      c=(1.1910627e-5/x)/x**4
      entry black(t)
      alpha=a/t
      if(alpha-2.e-4.lt.0.) goto 2
      planck=c/(exp(alpha)-1.)
   1  black=planck
      return
   2  planck=c/alpha
      goto 1
      entry dplnck(dt,dx)
      da=1.43880d0/dx
      dc=(1.1910627d-5/dx)/dx**4
      entry dblack(dt)  
      dalpha=da/dt  
      if(dalpha-1.d-8.lt.0.d0) goto 4
      dplnck=dc/(exp(dalpha)-1.0) 
   3  dblack=dplnck 
      return
   4  dplnck=dc/dalpha  
      goto 3
      entry tplck(p,x)  
      tplck=(1.43880/x)/alog(1.+(1.1910627e-5/x)/(p*x**4))  
      return
      entry dtplck(dp,dx)   
      dtplck=(1.43880d0/dx)/dlog(1.d0+(1.1910627d-5/dx)/(dp*dx**4)) 
      return
      end   

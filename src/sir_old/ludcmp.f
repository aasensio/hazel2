      subroutine ludcmp(a,n,np,indx,d)

	implicit real*4 (a-h,o-z)

      include 'PARAMETER'  !por nmax==kn
c      parameter (nmax=100,tiny=1.0d-20)
c      dimension a(np,np),indx(n),vv(nmax)

      parameter (tiny=1.0d-20)
      dimension a(np,np),indx(n),vv(kn)

      d=1.d0

      do 12 i=1,n

        aamax=0.d0

        do 11 j=1,n

          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))

11      continue

C        if (aamax.eq.0.d0) pause 'singular matrix.'
        if (aamax.eq.0.d0) aamax=1.e-6


        vv(i)=1.d0/aamax

12    continue

      do 19 j=1,n

        if (j.gt.1) then

          do 14 i=1,j-1

            sum=a(i,j)

            if (i.gt.1)then

              do 13 k=1,i-1

                sum=sum-a(i,k)*a(k,j)

13            continue

              a(i,j)=sum

            endif

14        continue

        endif

        aamax=0.d0

        do 16 i=j,n

          sum=a(i,j)

          if (j.gt.1)then

            do 15 k=1,j-1

              sum=sum-a(i,k)*a(k,j)

15          continue

            a(i,j)=sum

          endif

          dum=vv(i)*abs(sum)

          if (dum.ge.aamax) then

            imax=i

            aamax=dum

          endif

16      continue

        if (j.ne.imax)then

          do 17 k=1,n

            dum=a(imax,k)

            a(imax,k)=a(j,k)

            a(j,k)=dum

17        continue

          d=-d

          vv(imax)=vv(j)

        endif

        indx(j)=imax

        if(j.ne.n)then

          if(a(j,j).eq.0.d0)a(j,j)=tiny

          dum=1./a(j,j)

          do 18 i=j+1,n

            a(i,j)=a(i,j)*dum

18        continue

        endif

19    continue

      if(a(n,n).eq.0.d0)a(n,n)=tiny

      return

      end


c     output: m is the number of R^2 below the threshold d
      SUBROUTINE oddsratio(x, y, u, w, a, al, au, q, d, m, wind)
     
      integer x, y, n11, n10, n01, n00, m, wind((x-1)*(y-1),2)
      double precision a, al, au, q, u(x-1, y-1), w(x-1, y-1)
      double precision s, v, d, tl, tu 

      i = 1
      do 10 n01=1, x-1
       do 20 n11=1, y-1
       n10 = y - n11
      n00 = x - n01
      v = dble(n11*n00)/dble(n10*n01)
      u(n01,n11) = (v - a)**2
      s = dexp(q*dsqrt(1.0D0/n11 + 1.0D0/n10 + 1.0D0/n01 + 1.0D0/n00))
      tl = v / s
      tu = v * s
      w(n01,n11) = (tl - al)**2
      w(n01, n11) = w(n01, n11) + u(n01, n11)
      if(w(n01, n11) .LE. d)then
      wind(i, 1) = n01
      wind(i, 2) = n11
      i = i + 1
      endif
  20  continue    
  10  continue
      return
      end
      
      SUBROUTINE sortw(x, w, ind)
     
      integer x, ind
      double precision w(x)

      call SORTRX(x, w, ind)
      return
      end


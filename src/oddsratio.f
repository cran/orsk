C input: x, y, a, al, au, q, d, t
C output
c     w: squared error
c     wind: index corresponding to squared error

      SUBROUTINE oddsratio(x, y, a, al, au, q, d, t, w, wind)
      
      implicit none
      integer i, x, y, n11, n10, n01, n00, wind((x-1)*(y-1),2), t
      double precision d, a, al, au, q, u, w(x-1, y-1)
      double precision s, v

      i = 1
      do 10 n01=1, x-1
            n00 = x - n01
         do 20 n11=1, y-1
            n10 = y - n11
            v = dlog(dble(n11*n00)/dble(n10*n01))
            u = (v - dlog(a))**2
            s = q*dsqrt(1.0D0/n11 + 1.0D0/n10 + 1.0D0/n01 + 1.0D0/n00)
               w(n01, n11) = (v - s - dlog(al))**2
            if(t.EQ.1) then
               w(n01, n11) = (v - s - dlog(al))**2
            else
               if(t.EQ.2) then
                  w(n01, n11) = (v + s - dlog(au))**2
               else 
                  if(t.EQ.3) then        
                     w(n01, n11)=(v-s-dlog(al))**2+(v+s-dlog(au))**2
                  endif
               endif
            endif
            w(n01, n11) = w(n01, n11) + u
            if(w(n01, n11) .LE. d)then
               wind(i, 1) = n01
               wind(i, 2) = n11
               i = i + 1
            endif
 20      continue    
 10   continue
      return
      end
      

      implicit real*8 (a-h,o-z)
      character*80, str
      dimension ff(999)
      call getarg (1, str)
      open (7, file=str)
      call getarg (2, str)
      n0 = 0
20    read (7,*,end=30) nn
      n0 = n0 + nn
      goto 20
30    rewind 7
      open (8, file=str)
      if (iargc().ge.3) then
        call getarg (3, str)
        read (str, *) nf
      else
        nf = 3
      endif
      write (6,*) n0, nf
 10   read (7, *, end= 100) np, fx, fy, fz, sx, sy, sz
      np1 = np -1
      dx = (sx - fx)/np1
      dy = (sy - fy)/np1
      dz = (sz - fz)/np1
      do i = 1, np
        read (8, *) (ff(j),j=1,nf)
        if (i.eq.1.or.i.eq.np) then
          write (6, 500) fx, fy, fz, (ff(j)*1.e-12, j=1,nf)
        else
          write (6, 600) fx, fy, fz, (ff(j)*1.e-12, j=1,nf)
        endif
        fx = fx + dx
        fy = fy + dy
        fz = fz + dz
      enddo
      write (6,*)
      write (6,*)
      goto 10
100   continue
500   format (' 1 ', 3f10.6, 999f9.4)
600   format (' 0 ', 3f10.6, 999f9.4)
      end

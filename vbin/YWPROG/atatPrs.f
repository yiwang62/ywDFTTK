      program atat_of_phonon_elastic_constant
c      USE NUMERICAL_LIBRARIES
c
c -------- energy in Ryd, temperature in K, coordinate in a.u.
c
      implicit real*8(a-h,o-z)

      real*8, allocatable, target::v(:), f(:),
     &  XDATA(:), B(:), SSPOLY(:) 
      real*8 STAT(10)
      integer, allocatable, target::iperm(:)
      character*80 str, mode
      parameter (NSMOOTH = 10001)
      dimension vs(NSMOOTH), es(NSMOOTH)
c
      narg = iargc()
      call getarg(1, mode)
      if (mode.eq."-BM") then
        ibm = 1
        call getarg(2, str)
        read (str, *) vol0
        nv = (narg - 2) / 2
        i = 3
      else if (mode.eq."-LOG") then
        ibm = 0
        call getarg(2, str)
        read (str, *) vol0
        vol0 = log(vol0)
        nv = (narg - 2) / 2
        i = 3
      else
        ibm = 0
        read (mode, *) vol0
        nv = (narg - 1) / 2
        i = 2
      endif

      allocate (v(nv), f(nv),iperm(nv), XDATA(nv))

      j = 0
      do while (i.le.narg)
        j = j + 1
        iperm(j) = j
        call getarg(i, str)
        read (str, *) v(j)
        call getarg(i+1, str)
        read (str, *) XDATA(j)
        i = i + 2
        if (mode.eq."-LOG") then
          v(j) = log(v(j))
          XDATA(j) = log(XDATA(j))
        endif
      enddo

      CALL DSVRGP (nv, v, v, iperm)
      do i = 1, nv
        f(i) = XDATA(iperm(i))
      enddo

      if (ibm.eq.1) goto 1001

      dv = (v(nv) - v(1))/dble(NSMOOTH-1)
      do i = 1, NSMOOTH
        vs(i) = v(1) + dble(i-1)*dv
      enddo
      call DCSIEZ (nv, v, f, NSMOOTH, vs, es)
      pressure = -DQDDER(1, vol0,NSMOOTH,vs,es,.true.)
      goto 100

1001  allocate (B(nv), SSPOLY(nv))
      NDEG = min(3,nv-1)
      do i = 1,nv
        XDATA(i) = v(i)**(-2.d0/3.d0)
      enddo
      CALL DRCURV (nv, XDATA, f, NDEG, B, SSPOLY, STAT)
      x0 = vol0**(-2.d0/3.d0)
      x1 = (-2.d0/3.d0)*x0/vol0
      x2 = (-2.d0/3.d0)*(-5.d0/3.d0)*x0/vol0/vol0
      t1 = 0.d0
      t2 = 0.d0
      do i = 1, NDEG
        t1 = t1 + dble(i)*B(i+1)*x0**(i-1)
        t2 = t2 + dble(i-1)*dble(i)*B(i+1)*x0**(i-2)
      enddo
      d1 = t1*x1
      d2 = t2*x1*x1 + t1*x2
      pressure = -d1
      Bcld = vol0*d2
 100  write (6, '(f20.16)') pressure
c      write (0, *) "x0", x0, x1, x2, t1, t2, d1, d2
      end

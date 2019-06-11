      program atat_of_phonon
c      USE NUMERICAL_LIBRARIES
c
c -------- energy in Ryd, temperature in K, coordinate in a.u.
c
      implicit real*8(a-h,o-z)
      common/stktrp/KITRP

      real*8, allocatable, target:: 
     &  t(:),  vmin(:), fmin(:), smin(:), hmin(:), bmin(:),
     &  fion_min(:), sion_min(:), hion_min(:), 
     &  v(:), f(:, :), s(:, :),
     &  fion(:, :), sion(:, :),
     &  cvion(:, :), hion(:, :), 
     &  felec(:, :), selec(:, :), cvelec(:, :),
     &  fatat(:), satat(:), hatat(:),
     &  runv(:), runf(:), fdum(:), sdum(:), cdum(:), 
     &  veclat(:,:), latva(:,:)

      dimension tmpv(9)

      character*1024 wkdir, fn, mode*16, doskey*16
c
      UtoJ = dconst('ElectronVolt')*dconst('Avogadro')
      UtoGPa = dconst('ElectronVolt')*1.0d21
      Pressure = 0.d0
      NRUNX = 10001
      KITRP = 7
      mode = "vdos_e"
      doskey = "int2"

      narg = iargc()
      i = 1
      do while (i.le.narg)
	  call getarg(i, wkdir)
	  if (wkdir(1:4).eq."vdos".or.
     &      wkdir(1:4).eq."fvib".or.
     &      wkdir(1:4).eq."fitf") then
	    mode = wkdir
	  else if (wkdir(1:4).eq."-int") then
	    doskey = wkdir(2:5)
	  else if (wkdir(1:2).eq."-Pr") then
            i = i + 1
	    call getarg(i, wkdir)
            read (wkdir, *) Pressure
            Pressure = Pressure/UtoGPa
	  endif
          i = i + 1
      enddo

      read (5, *) natom, nt, nv
      allocate (
     &  t(nt),  vmin(nt), fmin(nt), smin(nt), hmin(nt), bmin(nt),
     &  fion_min(nt), sion_min(nt), hion_min(nt), 
     &  v(nv), f(nv, nt), s(nv, nt),
     &  fion(nv, nt), sion(nv, nt), 
     &  cvion(nv, nt), hion(nv, nt), 
     &  felec(nv, nt), selec(nv, nt), cvelec(nv, nt), 
     &  fatat(nt), satat(nt), hatat(nt), 
     &  veclat(nv, 9), latva(nt,9),
     &  runv(NRUNX), runf(NRUNX), 
     &  fdum(nt), sdum(nt), cdum(nt), STAT=istat)
      if (istat.ne.0) then 
        print *, "Cannot allocate main dynamic memory"
        stop
      endif 

      do i = 1, nt
        read (5,*) t(i), fatat(i), satat(i)
        fatat(i) = fatat(i)/dble(natom)
        satat(i) = satat(i)/dble(natom)
      enddo

      icalat = 1 
      nv = 0
 10   nv = nv+1
      read (5, *, end=100) temp, ee0, wkdir
      v(nv) = temp
      lwk = istrlen(wkdir)

      if (mode(1:4).eq."vdos") then
        fn = wkdir(1:lwk)//"/vdos.out"
	  if (doskey(1:4).eq."int0") then
            call sum_d_free_energy(nt, t, fn, fdum, sdum, cdum)
	  else if (doskey(1:4).eq."int1") then
            call int_d_free_energy
     &        (nt, t, fn, fdum, sdum, cdum, natom, 0)
	  else if (doskey(1:4).eq."int2") then
            call int_d_free_energy
     &        (nt, t, fn, fdum, sdum, cdum, natom, 1)
	  endif
        do i = 1, nt
          fion(nv,i) = (ee0 + fdum(i))/dble(natom)
          sion(nv,i) = (sdum(i))/dble(natom)
          hion(nv,i) = fion(nv,i) + t(i)*sion(nv,i)
          cvion(nv,i) = (cdum(i))/dble(natom)
        enddo
      else if (mode(1:4).eq."fitf") then
        fn = wkdir(1:lwk)//"/fitfc.log"
        call f_free_energy(nt, t, fn, fdum)
        do i = 1, nt
          fion(nv,i) = ee0/dble(natom)+fdum(i)
        enddo
        do i = 1, nt
	    sion(nv,i) = -ymDDER(t(i),nt,t,fdum, 1)
            hion(nv,i) = fion(nv,i) + t(i)*sion(nv,i)
	    cdum(i) = hion(nv,i)
        enddo
        do i = 1, nt
	    cvion(nv,i) = ymDDER(t(i),nt,t,cdum, 1)
        enddo
      else if (mode(1:4).eq."fvib") then
        fn = wkdir(1:lwk)//"/fvib"
        open (7, file=fn)
        read (7,*) (fdum(i), i = 1,nt)
        do i = 1, nt
          fion(nv,i) = (ee0 + fdum(i))/dble(natom)
        enddo
        close (7)
        do i = 1, nt
	    sion(nv,i) = -ymDDER(t(i),nt,t,fdum, 1)/dble(natom)
            hion(nv,i) = fion(nv,i) + t(i)*sion(nv,i)
	    cdum(i) = hion(nv,i)
        enddo
        do i = 1, nt
	    cvion(nv,i) = ymDDER(t(i),nt,t,cdum, 1)
        enddo
      endif

      v(nv) = v(nv)/dble(natom)

      if (mode(5:6).eq."_e") then
        fn = wkdir(1:lwk)//"/fvib_ele"
        open (7, file=fn)
        do i = 1, nt
          read (7,*) tmp0, tmpf, tmps
          felec(nv,i) = tmpf/dble(natom)
          selec(nv,i) = tmps/dble(natom)
          cdum(i) = felec(nv,i) + t(i)*selec(nv,i)
        enddo
        do i = 1, nt
	    cvelec(nv,i) = ymDDER(t(i),nt,t,cdum, 1)
        enddo
        close (7)
      endif

      if (icalat.eq.1) then
        fn = wkdir(1:lwk)//"/POSCAR.static"
        open (7, file=fn, status="old", err=50)
        goto 60
 50     fn = wkdir(1:lwk)//"/POSCAR"
        open (7, file=fn, status="old", err=70)
 60     read (7,*)
        read (7,*) scale
        read (7,*) tmpv
        do i = 1, 9
          veclat(nv, i) = scale*tmpv(i)
        enddo
        close(7)
        goto 10
 70     icalat = 0
      endif
      goto 10
 100  nv = nv - 1

      do i = 1, nt
        do j = 1, nv
          if (mode(5:6).eq."_e") then
            f(j,i) = fion(j,i) + felec(j,i) + Pressure*v(j)
            s(j,i) = sion(j,i) + selec(j,i) + Pressure*v(j)
          else
            f(j,i) = fion(j,i) + Pressure*v(j)
            s(j,i) = sion(j,i) + Pressure*v(j)
          endif
        enddo
      enddo

      dv = (v(nv) - v(1))/dble(NRUNX-1) 
      do i = 1, NRUNX
	  runv(i) = v(1) + dble(i-1)*dv
      enddo
      do i = 1, nt
        if (nv.gt.1) then
	  call DCSIEZ (nv, v, f(1,i),NRUNX, runv, runf)
          vmin(i) = v_min(NRUNX, runv, runf, fmin(i))
        else
          vmin(i) = v(1)
          fmin(i) = f(1,i)
        endif
      enddo
      do i = 1, nt
        if (nv.gt.1) then
	    call DCSIEZ (nv, v, s(1,i), 1, vmin(i), smin(i))
	    call DCSIEZ (nv, v, fion(1,i), 1, vmin(i), fion_min(i))
	    call DCSIEZ (nv, v, sion(1,i), 1, vmin(i), sion_min(i))
        else
          smin(i) = s(1,i)
          fatat(i) = fatat(i) + ee0/dble(natom)
        endif
	  hmin(i) = fmin(i) + t(i)*smin(i)
	  hion_min(i) = fion_min(i) + t(i)*sion_min(i)
	  hatat(i) = fatat(i) + t(i)*satat(i)
      enddo
      write (6, 500)
 500  format("# T(K)     volome         F(eV)   S(J/K)     ",
     &       "    H(J/K)  a(-6/K) Cp(J/mol)   Cvion    Cpio",
     &       "n     Sion           Hion     T-D(K)       Cv",
     &       "  Bt(GPa)      dCp   Gamma")
      write (6, 510) Pressure*UtoGPa
 510  format("# Pressure = ", f12.4, "  GPa")
 	do i = 1, nt

	  vderv = DQDDER(1,t(i),nt,t,vmin,.true.)/vmin(i)/3.d0
	  Cp = DQDDER(1,t(i),nt,t,hmin,.true.)
          if (nv.gt.1) then
	    call DCSIEZ (nv, v, cvion(1,i), 1, vmin(i), Cv_ion)
	    call DCSIEZ (nv, v, sion(1,i), 1, vmin(i), S_ion)
	    call DCSIEZ (nv, v, cvelec(1,i), 1, vmin(i), Cv_ele)
	  else
	    Cv_ion = cvion(1,i)
	    S_ion = sion(1,i)
	    Cv_ele = cvelec(1,i)
	  endif
	  Cp_ion = ymDDER(t(i),nt,t,hion_min, 1)
	  Cv = Cv_ion + Cv_ele

          tdeb = debyet(Cv_ion, t(i))

	  if (nv.le.2) then
	    Bt = 0.d0
	  else
	    call DCSIEZ (nv, v, f(1,i),NRUNX, runv, runf)
c	    Bt = vmin(i)*DQDDER(2, vmin(i),NRUNX, runv, runf,.true.)
            CALL ERSET (3, 0, -1)
	    Bt = vmin(i)*ymBDER(nv, v, f(1,i), vmin(i), 4)
            CALL ERSET (3, 2,  2)
	  endif
          bmin(i) = Bt

	  dCp = (vderv*3.d0)**2*Bt*t(i)*vmin(i)
	  if (Cv.gt.0.d0) then
            ggam = vderv*3.d0*Bt*vmin(i)/Cv
	  else
	    ggam = 0.d0
	  endif
	  if (abs(ggam).ge.100.d0) ggam = 0.0d0

	  write (6, '(f6.1, (f11.6,f14.8,f9.4 f15.1, 
     &              2f9.4,
     &              2f9.4, f9.4, f15.1,
     &              f11.4, 
     &              f9.4, f9.4, f9.4, f8.4))')
     &        t(i), vmin(i), fmin(i), smin(i)*UtoJ, hmin(i)*UtoJ,
     &        vderv*1.e6, Cp*UtoJ,  
     &        Cv_ion*UtoJ, Cp_ion*UtoJ, S_ion*UtoJ, hion_min(i)*UtoJ, 
     &        tdeb,  
     &        Cv*UtoJ, Bt*UtoGPa, dCp*UtoJ, ggam
	enddo
      open (7, file=mode(1:len_trim(mode))//"_svib")
      open (8, file=mode(1:len_trim(mode))//"_fvib")
      write (8, "(i6,i6,99f14.8)") nv, nt, (v(j), j=1,nv)
      do i = 1,nt
        write (7, "(f6.1,99f11.4)") t(i), (s(j,i)*UtoJ, j=1,nv)
        write (8, "(f6.1,99f14.8)") t(i), 
     &    vmin(i), fmin(i), bmin(i), (f(j,i), j=1,nv)
      enddo
      close (7)
      close (8)
      if (icalat.eq.1) then
        do i = 1,nt
	  do j = 1, 9
	    call DCSIEZ (nv, v, veclat(1,j), 1, vmin(i), latva(i,j))
          enddo
c	  latva(i) = sqrt(tmpv(1)**2+tmpv(2)**2+tmpv(3)**2)
c	  latvb(i) = sqrt(tmpv(4)**2+tmpv(5)**2+tmpv(6)**2)
c	  latvc(i) = sqrt(tmpv(7)**2+tmpv(8)**2+tmpv(9)**2)
        enddo
        open (7, file=mode(1:len_trim(mode))//"_lat")
        write (7, '("# T(K)",9("  a(",i1,",",i1,")"," (-6/K)"))')
     &    ((i,j,j=1,3),i=1,3)
        do i = 1,nt
	  do j = 1, 9
	    tmp = DQDDER(1,t(i),nt,t,latva(1,j),.true.)
            if (abs(latva(i,j)).ge.0.00001d0) then
	       tmpv(j) = tmp/latva(i,j)
            else
	       tmpv(j) = 0.d0
            endif
          enddo
          write (7, "(f6.1,9(f8.4, f7.2))") t(i), 
     &      (latva(i,j), tmpv(j)*1.e6, j=1,9)

c	  aderv = DQDDER(1,t(i),nt,t,latva,.true.)/latva(i)
c	  bderv = DQDDER(1,t(i),nt,t,latvb,.true.)/latvb(i)
c	  cderv = DQDDER(1,t(i),nt,t,latvc,.true.)/latvc(i)
c          write (7, "(f6.1,3f10.6, 3f9.4)") t(i), 
c     &      latva(i),   latvb(i),   latvc(i),
c     &      aderv*1.e6, bderv*1.e6, cderv*1.e6
        enddo
        close (7)
      endif
      end
c
      subroutine sum_d_free_energy(nt, t, fndos, fdum, sdum, cdum)
c
      implicit real*8(a-h,o-z)
      character*1024 fndos
      dimension t(nt), fdum(nt), sdum(nt), cdum(nt)
      real*8, allocatable:: fdos(:)

      open (7, file=fndos)
      n = 0
      do while (.TRUE.)
        read (7, *, end=100) freq, dum
        n = n + 1
      enddo
100   allocate (fdos(n), STAT=istat)
      if (istat.ne.0) then 
        print *, "Cannot allocate dynamic memory for fdos"
        stop
      endif 

      rewind(7)
      w0 = 0.d0
      do i = 1, n
        read (7, *) freq, fdos(i)
        w0 = w0 + fdos(i)
      enddo
200   close(7)

      xfreq = freq/dble(n-1)
      do i = 1, n
        fdos(i) = fdos(i)*xfreq
      enddo
c      print *, "w0*dfreq = ", w0*xfreq
      dfreq = xfreq*dconst('Planck')
      do j = 1, nt
        temper = t(j)
	tkb = temper*dCONST('Boltzman')
	f0 = 0.d0
	s0 = 0.d0
	c0 = 0.d0
	  do i = 1, n
	      if (temper.gt.0.d0) t0 = dfreq*dble(i-1)
	      if (t0.gt.0.d0.and.temper.gt.0.d0) then
	    x = t0/tkb
	    xd2 = 0.5d0*x 
		f0 = f0 + fdos(i)*tkb*
     &            log(2.d0*sinh(xd2))
		s0 = s0 + 
     &            (x/(exp(x)-1.d0) - log(1.d0 
     &            - exp(-x)))*fdos(i)
		c0 = c0 + (xd2/sinh(xd2))**2*fdos(i)
	      else if (temper.eq.0.d0) then
	        f0 = f0 + 0.5d0*fdos(i)*dfreq*dble(i-1)
	      endif
	  enddo
        fdum(j) = f0/dconst('ElectronVolt')
        sdum(j) = s0*dCONST('Boltzman')/dconst('ElectronVolt')
        cdum(j) = c0*dCONST('Boltzman')/dconst('ElectronVolt')
      enddo
      deallocate(fdos)
      end
c
      subroutine int_d_free_energy(nt, t, fndos, 
     &           fdum, sdum, cdum, natom, imode)
c
      implicit real*8(a-h,o-z)
      character*1024 fndos
      dimension t(nt), fdum(nt), sdum(nt), cdum(nt)
      real*8, allocatable:: fdos(:), rdos(:),
     &  rx(:), ados(:), a(:), fx(:), sx(:), cx(:) 
      NRUNX = 1001

      open (7, file=fndos)
      n = 0
      do while (.TRUE.)
        read (7, *, end=100) freq, dum
        n = n + 1
      enddo
100   allocate (rdos(n), fdos(n), 
     &  rx(NRUNX), ados(NRUNX), a(NRUNX), 
     &  fx(NRUNX), sx(NRUNX), cx(NRUNX), 
     &  STAT=istat)
      if (istat.ne.0) then 
        print *, "Cannot allocate dynamic memory for fdos"
        stop
      endif 

      rewind(7)
      do i = 1, n
        read (7, *) rdos(i), fdos(i)
      enddo
200   close(7)
      if (rdos(2).lt.0.d0) write (0, '(a,a,a)') "******** ", 
     &    trim(fndos), " imaginary mode found!"
      df = (rdos(n)-rdos(1))/dble(n-1)
      do i = 1, n
        rdos(i) = (rdos(1) + dble(i-1)*df)
	  fdos(i) = fdos(i)
      enddo
      df = rdos(n)/dble(NRUNX-1)
      do i = 1, NRUNX
        rx(i) = dble(i-1)*df
      enddo
	call DCSIEZ (n, rdos, fdos, NRUNX, rx, ados)
	ados(1) = 0.d0

      call integr (ados, rx, NRUNX, a, 3)
      w0 = 3.d0*dble(natom)/a(NRUNX)
      write (0, '(a, f11.6, a, f11.6)') "******** Sum DOS = ",
     &    a(NRUNX), ",      Ideal DOS = ", 3.d0*dble(natom)
      if (imode.eq.1) then
        do i = 1, NRUNX
          ados(i) = ados(i)*w0
        enddo
      endif
c	call integr (ados, rx, NRUNX, a, 3)

      h = dconst('Planck')
      do j = 1, nt
        temper = t(j)
	  tkb = temper*dCONST('Boltzman')
	  do i = 1, NRUNX
	    hmu = h*rx(i)
	    fx(i) = 0.d0
	    sx(i) = 0.d0
	    cx(i) = 0.d0
	    if (hmu.gt.0.d0.and.temper.gt.0.d0) then
	      x = hmu/tkb
	      xd2 = 0.5d0*x 
		  fx(i) = tkb*log(2.d0*sinh(xd2))*ados(i)
		  sx(i) = (x/(exp(x)-1.d0) - log(1.d0 - exp(-x)))*ados(i)
		  cx(i) = (xd2/sinh(xd2))**2*ados(i)
	    else if (hmu.gt.0.d0.and.temper.eq.0.d0) then
	      fx(i) = 0.5d0*hmu*ados(i)
	    endif
	  enddo
	  call integr (fx, rx, NRUNX, a, 3)
        fdum(j) = a(NRUNX)/dconst('ElectronVolt')
	  call integr (sx, rx, NRUNX, a, 3)
        sdum(j) = a(NRUNX)*dCONST('Boltzman')/dconst('ElectronVolt')
	  call integr (cx, rx, NRUNX, a, 3)
        cdum(j) = a(NRUNX)*dCONST('Boltzman')/dconst('ElectronVolt')
      enddo
      deallocate(rdos, fdos, rx, ados, a, fx, sx, cx)
      end
c
      subroutine f_free_energy(nt, t, fnfrq, fdum)
c
      implicit real*8(a-h,o-z)
        character*1024 fnfrq, line
        dimension t(nt), fdum(nt)
      open (7, file=fnfrq)
      line = ""
      do while (line.ne."Phonon frequencies:")
        read (7, '(a)', end=100) line
      enddo

      do i = 1, nt
        fdum(i) = 0.d0
      enddo
      n = 0
      read (7, '(a)', end=100) line
      do while (line.ne."end")
        read (line, *) freq
        n = n + 1
          do i = 1, nt
              temper = t(i)
              tkb = temper*dCONST('Boltzman')
              if (temper.gt.0.d0) t0 = freq*dconst('Planck')
              if (t0.gt.0.d0.and.temper.gt.0.d0) then
                fdum(i) = fdum(i) + tkb*log(2.d0*sinh(0.5d0*t0/tkb))
              else if (temper.eq.0.d0) then
                fdum(i) = fdum(i) + 0.5d0*t0
              endif
          enddo
        read (7, '(a)', end=100) line
      enddo
100   close(7)
      do i = 1, nt
        fdum(i) = 3.d0*fdum(i)/dconst('ElectronVolt')/dble(n)
      enddo
      end
c     
      real*8 function ymBDER(nv, v, f, vol0, NORDER)
c     
      implicit real*8(a-h,o-z)
      dimension v(nv), f(nv),
     &  XDATA(nv), B(nv), SSPOLY(nv), STAT(10)
      NDEG = min(NORDER,nv-1)
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
      ymBDER = d2
c      write (0,*) "x0", ymBDER
      end
c
	real*8 function ymDDER(tx, nt, t, f, IDERIV)
c
      implicit real*8(a-h,o-z)
      dimension t(nt), f(nt)
      dimension runv(10001), runf(10001)
      if (IDERIV.eq.0) then
        call DCSIEZ (nt, t, f, 1, tx, value)
        ymDDER = value
        return
      endif
      NRUNX = 10001
      dv = (t(nt) - t(1))/dble(NRUNX -1)
      do i = 1, NRUNX
        runv(i) = t(1) + dble(i-1)*dv
      enddo
      call DCSIEZ (nt, t, f, NRUNX, runv, runf)
      ymDDER = DQDDER(IDERIV,tx,NRUNX,runv,runf,.true.)
      end
c
	real*8 function CSI_ymDDER(tx, nt, t, f, IDERIV, idummy)
c
      implicit real*8(a-h,o-z)
	dimension t(nt), f(nt)
	real*8, allocatable:: BSCOEF(:), XKNOT(:)
	KORDER = 3
      NKNOT = nt + KORDER
	allocate (BSCOEF(nt), XKNOT(NKNOT), STAT=istat)
      if (istat.ne.0) then 
        print *, "Cannot allocate dynamic memory in ymDDER"
        stop
      endif
C                                  Generate knot sequence
      CALL DBSNAK (nt, t, KORDER, XKNOT)
C                                  Interpolate
      CALL DBSINT (nt, t, f, KORDER, XKNOT, BSCOEF)
      CSI_ymDDER = DBSDER(IDERIV,tx,KORDER,XKNOT,nt,BSCOEF)
	deallocate (BSCOEF, XKNOT)
	end
c
	real*8 function old_ymDDER(tx, nt, t, f, NORDER, NSODX)
c
      implicit real*8(a-h,o-z)
	dimension t(nt), f(nt), a(0:NSODX,0:NSODX), b(0:NSODX), x(0:NSODX)
	dt = tx*0.15d0
	tb = tx - dt
	te = tx + dt
	tb = max(tb,t(1))
	te = min(te,t(nt))
	do i = 1, nt
	  if (t(i).gt.tb) goto 10
	enddo
	i = nt
 10   ib = i - 1
      ib = max (ib, 1)
	do i = 1, nt
	  if (t(i).ge.te) goto 20
	enddo
	i = nt
 20   ie = i
c      if ((ie-ib).le.4) ie = ib +5
      if ((ie-ib).le.4) then
        old_ymDDER = DQDDER(NORDER,tx,nt,t,f,.true.)
	  return
	endif
	do m = 0, NSODX
        b(m) = 0.d0	  
	  do i = ib, ie
	    b(m) = b(m) + f(i)*(t(i)-tx)**m
	  enddo
	  do n = 0, NSODX
	    a(n, m) = 0.d0
	    do i = ib, ie
	      a(n, m) = a(n,m) + (t(i)-tx)**n*(t(i)-tx)**m
	    enddo
	  enddo
	enddo
	CALL DLSLRG (NSODX+1, a, NSODX+1, b, 1, x)
	old_ymDDER = x(NORDER)
	do i = 1, NORDER
	  old_ymDDER = old_ymDDER/dble(i)
	enddo
	end
c
      real*8 function debyet(cv, temper)
c
      implicit real*8(a-h,o-z)
      common/stack/rcv
	external fdebye
	data zero /0.d0/
	if (temper.le.zero) then
	  debyet = zero
	  return
	endif
	rcv = cv/3.d0/dCONST('Boltzman')*dconst('ElectronVolt')

	A = 1.d0
	vA = fdebye(A/temper)
	do B = A, 100000.1d0, 10.d0
	  vB = fdebye(B/temper)
	  if (vA*vB.le.zero) goto 100
	  A = B
	  vA = vB
	enddo
	debyet = zero
	return

 100    ERRABS = 0.d0
	  ERRREL = 1.d-12
	  MAXFN = 999
	  A = A/temper
	  B = B/temper
        CALL DZBREN (fdebye, ERRABS, ERRREL, A, B, MAXFN)
	  debyet = B*temper
      end
c
c
	real*8 function fdebye(x)
c
      implicit real*8(a-h,o-z)
      common/stack/rcv
      external fdx
      A = 0.0d0
      ERRABS = 0.0d0
      ERRREL = 1.d-12
      CALL DQDAGS (fdx, A, x, ERRABS, ERRREL, r, ERREST)
	fdebye = 3.d0*(1.d0/x)**3*r - rcv
	return
	end
c
	real*8 function fdx(x)
c
      implicit real*8(a-h,o-z)
	fdx = exp(x)*x**4/(exp(x)-1.d0)**2 
	return
	end
c
	real*8 function rootvl(root)
c
      implicit real*8(a-h,o-z)
      common/rootst/p(7),x(7)
      common/stktrp/KITRP
      call interp(x,p,KITRP,root,vl,dum,.false.)
	rootvl = vl
	return
	end
c
      real*8 function v_min(n, runv, runf, fmin)
c
      implicit real*8(a-h,o-z)
      dimension runv(n), runf(n),dumf(7)
      common/stktrp/KITRP
      common/rootst/p(7),x(7)
	external rootvl
	data zero /0.d0/
      pmin = 1.e36
      do i = 1, n
        if (pmin.gt.runf(i)) then
          k = i
          pmin = runf(i)
        endif
      enddo
	k = max(4,k)
	k = min (k, n-3)
      do i=1,KITRP
        x(i) = runv(k+i-KITRP/2-1)
        p(i)     = DQDDER(1,x(i),n,runv,runf,.true.)
        dumf(i) = runf(k+i-KITRP/2-1)
      enddo

	A = x(1)
	vA = p(1)
	do i=1, KITRP
	  if(p(i)*vA.gt.0.d0) then
	    A = x(i)
	    vA = p(i)
	  endif
	enddo
	B = x(KITRP)
	vB = p(KITRP)
c	A = A - 1.e-2
c	B = B + 1.e-2
	do i=KITRP, 1, -1
	  if(p(i)*vB.gt.0.d0) then
	    B = x(i)
	    vB = p(i)
	  endif
	enddo

	if (1.eq.0.or.vA*vB.gt.0.d0) then
        call interp(p,x,KITRP,zero,xzero,dum,.false.)
	  if (xzero.lt.x(1).or.xzero.gt.x(KITRP)) xzero=x(KITRP/2+1)
	else
	  ERRABS = 0.d0
	  ERRREL = 1.d-8
	  MAXFN = 999
        CALL DZBREN (rootvl, ERRABS, ERRREL, A, B, MAXFN)
	  xzero = B
	endif
      v_min = xzero
      call DCSIEZ (KITRP, x, dumf, 1, xzero, fmin)
      end
c
      integer function istrlen(s)
c
      character*(*) s
      ls = len(s)      
      do i = 1, ls
        if (ichar(s(i:i)).le.32) goto 100
      enddo
      i = ls
 100  istrlen = i - 1
      return
      end
c
      subroutine interp(r,p,n,rs,ps,dps,deriv)
c
c***********************************************************************
c
c        interpolates function p on the 1-d mesh r by newtonian divided
c      differences.  see: apostol v.2, eq.8.60.
c        rs is the pt to be interpolated,ps is the interpolated value
c      of p at rs.
c        dif must be dimensioned .ge. n, # interp points.
c
c      called by: nrmliz,scf,tmat,vgen
c
c***********************************************************************
c
      implicit real*8(a-h,o-z)
      logical deriv,nodriv
      dimension r(1),p(1),dif(7),r1(7)
      data zero,one/0.0d0,1.0d0/
c
      nodriv=.not.deriv
      ps=p(1)
      dps=zero
      rprod=one
      nm1=n-1
c
c----- zeroth order differences ----------------------------------------
c
      do 10 i=1,n
      if (rs.eq.r(i)) go to 50
   10 dif(i)=p(i)
c
c----- cycle over difference orders ------------------------------------
c
      do 40 i=1,nm1
      nmi=n-i
      ind=i+1
c
c----- compute i th order differences & pack in lowest (n-i) of dif ----
c
      do 20 j=1,nmi
      jp1=j+1
      dif(j)=(dif(j)-dif(jp1))/(r(j)-r(ind))
   20 ind=ind+1
      rprod=rprod*(rs-r(i))
      radd=rprod*dif(1)
      if (nodriv) go to 40
      do 30 k=1,i
   30 dps=dps+radd/(rs-r(k))
   40 ps=ps+radd
      return
c
c----- this section executed if rs=r(i): (n-2)nd order interp ----------
c
   50 j=1
      do 70 i=1,n
      if (rs.eq.r(i)) go to 60
      dif(j)=p(i)
      r1(j)=r(i)
      j=j+1
      go to 70
   60 ps=p(i)
      if (nodriv) return
   70 continue
      nm2=n-2
      do 100 i=1,nm2
      nmi=nm1-i
      ind=i+1
      do 80 j=1,nmi
      jp1=j+1
      dif(j)=(dif(j)-dif(jp1))/(r1(j)-r1(ind))
   80 ind=ind+1
      rprod=rprod*(rs-r1(i))
      radd=rprod*dif(1)
      do 90 k=1,i
   90 dps=dps+radd/(rs-r1(k))
  100 continue
      return
      end
c
      subroutine integr(f,r,kmx,a,md)
c
c***********************************************************************
c
c        integrates function f on 1-d mesh r by quadratures, stores
c      partial integrals (integrals up to r(k)) in a.
c        see:  k s kunz, numerical analysis eq. 7.36 ff
c
c        modes:
c      md=1, integration starts at r=0 and a is not reversed.
c      md=2, integration starts at first mesh point and members
c             of array a are subtracted from a(kmax); a becomes that
c             fraction of the integral beyond r(k).
c      md=3, integration starts at first mesh point and a not reversed.
c
c***********************************************************************
c
      implicit real*8(a-h,o-z)
      dimension f(1),r(2),a(1)
      data s720,s646,s456,s346,s264,s251,s106,s74,s30,s25,s24,
     &          s19,s11,s10,s9,s7,s6,s5,s4,s2
     &     /720.d0,646.d0,456.d0,346.d0,264.d0,251.d0,106.d0,74.d0,
     &          30.d0,25.d0,24.d0,19.d0,11.d0,10.d0,9.d0,7.d0,6.d0,5.d0,
     &          4.d0,2.d0/
      data zero/0.0d0/
c
      h=r(2)-r(1)
      if (md.ne.1) go to 10
c
c----- begin integration at zero (function must =0 at 0) --------------
c
      k0=0
      f0=zero
      go to 20
c
c----- begin integration at first mesh point ---------------------------
c
   10 k0=1
      a(1)=zero
      f0=f(1)
   20 n=1
      km3 =k0+1
      km2 =k0+2
      km1 =k0+3
      kfst=k0+4
c
c----- quadrature formulas q41(0),q41(1),q41(2) for 1st 3 pts ----------
c
      a(km3)=h*(s251*f0+s646*f(km3)-s264*f(km2)+s106*f(km1)
     &        -s19*f(kfst))/s720
      a(km2)=a(km3)+h*(-s19*f0+s346*f(km3)+s456*f(km2)-s74*f(km1)
     &        +s11*f(kfst))/s720
      a(km1)=a(km2)+h*(s11*f0-s74*f(km3)+s456*f(km2)+s346*f(km1)
     &        -s19*f(kfst))/s720
c
c----- main part of integration: q31(2) --------------------------------
c
   30 klst=kmx
      do 40 k=kfst,klst
      a(k)=a(km1)+h*(s9*f(k)+s19*f(km1)-s5*f(km2)+f(km3))/s24
      km3=km2
      km2=km1
      km1=k
   40 continue
      if (mod(md,2).ne.0) return
c
c----- if md is even, reverse a ----------------------------------------
c
      do 60 k=1,kmx
   60 a(k)=a(kmx)-a(k)
      return
      end


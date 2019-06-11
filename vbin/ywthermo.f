      program atat_of_phonon
c      USE NUMERICAL_LIBRARIES
c
c -------- energy in Ryd, temperature in K, coordinate in a.u.
c
      implicit real*8(a-h,o-z)
      common/stktrp/KITRP
      common/rundis/NRUNX
      common/lowtem/toobig

      real*8, allocatable, target:: 
     &  t(:),  vmin(:), fmin(:), smin(:), hmin(:), bmin(:),
     &  Cp(:), Cv(:),  alpha(:),
     &  fion_min(:), sion_min(:), hion_min(:), 
     &  v(:), Ezero(:), f(:, :), s(:, :),
     &  runFion(:, :), fion(:, :), sion(:, :),
     &  Sspin(:), Espin(:),
     &  Fmag(:,:), Smag(:,:), Emag(:,:), Cmag(:,:),
     &  cvion(:, :), hion(:, :), 
     &  felec(:, :),selec(:, :),
     &  cvelec(:, :), gmu(:, :), exel(:, :),
     &  dVdT(:, :), dMudT(:, :), xgmu(:),
     &  elef(:),
     &  runv(:), runf(:), runp(:), fdum(:), sdum(:), cdum(:), 
     &  veclat(:,:), latva(:,:), f00(:), ffdum(:,:),
     &  vecCij(:,:) 

      dimension tmpv(9), revm(3,3), revm6(6), revlat(9), tmpCij(36),
     &  aseCij(36)

      character*1024 wkdir, contcar, cijfile,
     &  fn, mode*16, doskey*16, spinf*16
c
      call settim

      toobig = 1.d2
      UtoJ = dconst('ElectronVolt')*dconst('Avogadro')
      UtoGPa = dconst('ElectronVolt')*1.0d21
      Pressure = 0.d0
      NRUNX = 10001
      KITRP = 7
c     mode = "vdos_e"
      mode = "vdos"
      doskey = "int2"
      spinf = " "
      unitf = 1.e0;
      modeBt = 0
      modeFF = 0
      contcar = "Static.CON"
      cijfile = ""
      ihugoniot = 0
      izero = 0

      narg = iargc()
      i = 1
      do while (i.le.narg)
	  call getarg(i, wkdir)
	  if (wkdir(1:4).eq."vdos".or.
     &      wkdir(1:4).eq."fitf") then
	    mode = wkdir
	  else if (wkdir(1:4).eq."-Hug") then
	    ihugoniot = 1
	  else if (wkdir(1:3).eq."-NRUNX") then
            i = i + 1
	    call getarg(i, wkdir)
            read (wkdir, *) NRUNX
	  else if (wkdir(1:6).eq."-izero") then
            i = i + 1
	    call getarg(i, wkdir)
            read (wkdir, *) izero
	  else if (wkdir(1:4).eq."-int") then
	    doskey = wkdir(2:5)
	  else if (wkdir(1:4).eq."-THz") then
	    unitf = 1.e12
	  else if (wkdir(1:6).eq."-fitBt") then
	    modeBt = 1
	  else if (wkdir(1:5).eq."-fitF") then
	    modeFF = 1
	  else if (wkdir(1:4).eq."-mag") then
	    spinf = wkdir(2:5)
	  else if (wkdir(1:5).eq."-cont") then
            i = i + 1
	    call getarg(i, contcar)
	  else if (wkdir(1:5).eq."-cijf") then
            i = i + 1
	    call getarg(i, cijfile)
	  else if (wkdir(1:3).eq."-Pr") then
            i = i + 1
	    call getarg(i, wkdir)
            read (wkdir, *) Pressure
c      write (0, *) "Pre = ", Pressure
            Pressure = Pressure/UtoGPa
c      write (0, *) "new Pre = ", Pressure
	  endif
          i = i + 1
      enddo

      read (5, *) natom, nv, t0, t1, dt
      nt = (t1-t0)/dt + 1.5d0
      nv0 = nv

      allocate (
     &  t(nt),  vmin(nt), fmin(nt), smin(nt), hmin(nt), bmin(nt),
     &  Cp(nt), Cv(nt),  alpha(nt),
     &  fion_min(nt), sion_min(nt), hion_min(nt), 
     &  v(nv), Ezero(nv), f(nv, nt), s(nv, nt),
     &  runFion(nv, nt), 
     &  fion(nv, nt), sion(nv, nt),  Sspin(nv), Espin(nv),
     &  Fmag(nv, nt), Smag(nv, nt), Emag(nv, nt), Cmag(nv, nt), 
     &  cvion(nv, nt), hion(nv, nt), 
     &  felec(nv, nt), selec(nv, nt), cvelec(nv, nt), 
     &  gmu(nv, nt), exel(nv, nt), xgmu(nt),
     &  dVdT(nv, nt), dMudT(nv, nt),
     &  veclat(nv,9), latva(nt,9),
     &  vecCij(nv,36), ffdum(nv, nt),
     &  runv(NRUNX), runf(NRUNX), runp(NRUNX), elef(NRUNX),
     &  fdum(nt), sdum(nt), cdum(nt), f00(nv), STAT=istat)
      if (istat.ne.0) then 
        print *, "Cannot allocate main dynamic memory"
        stop
      endif 

      do i = 1, nt
        t(i) = t0
	t0 = t0 + dt
      enddo

      icalat = 1 
      nv = 0
 10   nv = nv+1
      if (spinf.eq."mag0") then
        read (5, *, end=100) temp, ee0, wkdir, Sspin(nv)
        Sspin(nv) = Sspin(nv)/dble(natom)
      else if (spinf.eq."mag1") then
        read (5, *, end=100) temp, ee0, wkdir, Sspin(nv), Espin(nv)
        Sspin(nv) = Sspin(nv)/dble(natom)
        Espin(nv) = Espin(nv)/dble(natom)
      else 
        read (5, *, end=100) temp, ee0, wkdir
      endif
      v(nv) = temp
      Ezero(nv) = ee0/dble(natom)
      if (mode(1:4).eq."vdos") then
        fn = trim(wkdir)//"/vdos.out"
	  if (doskey(1:4).eq."int0") then
            call sum_d_free_energy
     &        (nt, t, fn, fdum, sdum, cdum, natom, unitf)
	  else if (doskey(1:4).eq."int1") then
            call int_d_free_energy
     &        (nt, t, fn, fdum, sdum, cdum, natom, 0, unitf, izero)
	  else if (doskey(1:4).eq."int2") then
            call int_d_free_energy
     &        (nt, t, fn, fdum, sdum, cdum, natom, 1, unitf, izero)
	  endif
        do i = 1, nt
          if (t(i).eq.0.d0) f00(nv) = fdum(i)/dble(natom)
          fion(nv,i) = (ee0 + fdum(i))/dble(natom)
          sion(nv,i) = (sdum(i))/dble(natom)
          hion(nv,i) = fion(nv,i) + t(i)*sion(nv,i)
          cvion(nv,i) = (cdum(i))/dble(natom)
        enddo
      else if (mode(1:4).eq."fitf") then
        if (mode(1:4).eq."fitf") then
          fn = wkdir(1:lwk)//"/fitfc.log"
        endif
        call f_free_energy(nt, t, fn, fdum, sdum, cdum)
        do i = 1, nt
          if (t(i).eq.0.d0) f00(nv) = fdum(i)
          fion(nv,i) = ee0/dble(natom)+fdum(i)
          sion(nv,i) = sdum(i)
          hion(nv,i) = fion(nv,i) + t(i)*sion(nv,i)
          cvion(nv,i) = cdum(i)
        enddo
      endif

      do i=1, nt
        ffdum(nv, i) = fdum(i)/dble(natom)
      enddo

      v(nv) = v(nv)/dble(natom)
      if (mode(5:6).eq."_e") then
        fn = trim(wkdir)//"/fvib_ele"
        open (7, file=fn)
        do i = 1, nt
          read (7,*, err=41) tmp0, tmpf, tmps, tmpc, tmpgmu, tmpexel,
     &      tmpdVdT, tmpdMudT
          goto 42
 41       backspace 7
          read (7,*) tmp0, tmpf, tmps, tmpc
          tmpgmu = 0.0
          tmpexel = 0.0
          tmpdVdT = 0.0
          tmpdMudT = 0.0
 42       felec(nv,i) = tmpf/dble(natom)
          selec(nv,i) = tmps/dble(natom)
	  cvelec(nv,i) = tmpc/dble(natom)
	  exel(nv,i) = tmpexel/dble(natom)
	  gmu(nv,i) = tmpgmu
	  dVdT(nv,i) = tmpdVdT
	  dMudT(nv,i) = tmpdMudT
        enddo
        close (7)
      else
        do i = 1, nt
          felec(nv,i) = 0.e0
          selec(nv,i) = 0.e0
	  cvelec(nv,i) = 0.e0
	  exel(nv,i) = 0.e0
	  gmu(nv,i) = 0.e0
	  dVdT(nv,i) = 0.e0
	  dMudT(nv,i) = 0.e0
        enddo          
      endif

      if (icalat.eq.1) then
c        fn = wkdir(1:lwk)//"/"//contcar
        fn = trim(wkdir)//"/"//trim(contcar)
        open (7, file=fn, status="old", err=70)
        read (7,*)
        read (7,*) scale
        read (7,*) tmpv
        do i = 1, 9
          veclat(nv,i) = scale*tmpv(i)
        enddo
        close(7)
        goto 72
 70     icalat = 0
        write (6, 71) trim(fn)
 71     format (/"********CANNOT FIND file", A, 
     &    ", I will continue without calculating ",
     &    "thermal expansion tensor"/)
      endif

 72   if (trim(cijfile).ne."") then
        fn = trim(wkdir)//"/"//trim(cijfile)
        open (7, file=fn, status="old", err=80)
        read (7,*) tmpCij
        do i = 1, 36
          vecCij(nv,i) = tmpCij(i)
        enddo
c        write (0,*) tmpCij
        close(7)
        goto 10
 80     cijfile = ""
c        write (6, 81) trim(fn)
c 81     format (/"********CANNOT FINE file", A,
c     &    ", I will continue without calculating ",
c     &    "Cij(T)"/)
      endif

      goto 10
 100  nv = nv - 1

      if (modeFF.eq.1) then
        do i=1, nt
          call Least_Square(nv, 1, v, ffdum(1, i), tmpv)
c          write (0, *) tmpv(1), tmpv(2)
          do j=1, nv
            corr = tmpv(1) + v(j)*tmpv(2) - ffdum(j, i)
            runFion(j, i) = fion(j, i) + corr
c            fion(j, i) = fion(j, i) + corr
c            hion(j, i) = hion(j, i) + corr
          enddo
        enddo
      endif

      if (nv.ne.nv0) then
        write (0,*)
        write (0,*) "****i*** FETAL error, No. of volume=", nv, 
     &    " NOT=", nv0
        write (0,*) "****i*** FETAL error, No. of volume=", nv, 
     &    " NOT=", nv0
        write (0,*) "****i*** FETAL error, No. of volume=", nv, 
     &    " NOT=", nv0
        write (0,*)
        stop
      endif

      toeV = dCONST('Boltzman')/dconst('ElectronVolt')
      do i = 1, nt
        do j = 1, nv
          Fmag(j,i) = 0.d0
          Smag(j,i) = 0.d0
          Emag(j,i) = 0.d0
          Cmag(j,i) = 0.d0
          s(j,i) = sion(j,i)

          if (modeFF.eq.1) then
            f(j,i) = runFion(j,i) + Pressure*v(j)
          else
            f(j,i) = fion(j,i) + Pressure*v(j)
          endif
          if (spinf.eq."mag0") then
            tmp2j = Sspin(j)*(6.e0-Sspin(j))
c      write (0,*) "tmp2j=", tmp2j
            smag(j,i) = log(tmp2j+1.e0)*toeV
            s(j,i) = s(j,i) + smag(j,i)
            f(j,i) = f(j,i) - t(i)*smag(j,i)
          else if (spinf.eq."mag1") then
            tmp2j = Sspin(j)*(6.e0-Sspin(j))
            if (t(i).ne.0.0) then
              fx = tmp2j*exp(-Espin(j)/toeV/t(i))
	      Fmag(j,i) = -toeV*t(i)*log(1.d0+fx)
              Smag(j,i) = toeV*log(1.d0+fx)
     &          + Espin(j)/t(i)*fx/(1.d0+fx)
              Emag(j,i) = Fmag(j,i) + t(i)*Smag(j,i)
              Cmag(j,i) = toeV
     &          *(Espin(j)/toeV/t(i))**2
     &          *fx/(1.d0+fx)**2
            endif
c       write (0,*) Fmag(j,i), Emag(j,i), Smag(j,i), Cmag(j,i)
            s(j,i) = s(j,i) + Smag(j,i)
            f(j,i) = f(j,i) + Fmag(j,i)
          endif
          if (mode(5:6).eq."_e") then
            f(j,i) = f(j,i) + felec(j,i)
            s(j,i) = s(j,i) + selec(j,i)
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
c          write (0, *) "v_min=", t(i), vmin(i)
        else
          vmin(i) = v(1)
          fmin(i) = f(1,i)
        endif
      enddo
      do i = 1, nt
c          write (0, *) "v_min=", t(i), vmin(i)
        if (nv.gt.1) then
	    call DCSIEZ (nv, v, s(1,i), 1, vmin(i), smin(i))
	    call DCSIEZ (nv, v, fion(1,i), 1, vmin(i), fion_min(i))
	    call DCSIEZ (nv, v, sion(1,i), 1, vmin(i), sion_min(i))
	    call DCSIEZ (nv, v, gmu(1,i), 1, vmin(i), xgmu(i))
            if (modeFF.eq.1) then
	      call DCSIEZ (nv, v, runFion(1,i), 1, vmin(i), fcorr)
              fmin(i) = fmin(i) - fcorr +  fion_min(i)
            endif
        else
          smin(i) = s(1,i)
          fion_min(i) = fion(1,i)
          sion_min(i) = sion(1,i)
          xgmu(i) = gmu(1,i)
        endif
	  hmin(i) = fmin(i) + t(i)*smin(i)
	  hion_min(i) = fion_min(i) + t(i)*sion_min(i)
      enddo
      write (6, 500)
 500  format("# T(K)     volome         F(eV)   S(J/K)     ",
     &       "    H(J/K)      a(-6/K)     Cp(J/mol)         Cvion ",
     &       "        Cpion     Sion     Sele     Smag      ",
     &       "     Hion     T-D(K)           Cv  Bt(GPa)   ",
     &       "   dCp   Gamma")
      write (6, 510) Pressure*UtoGPa
 510  format("# Pressure = ", f12.4, "  GPa")
 	do i = 1, nt

	  vderv = DQDDER(1,t(i),nt,t,vmin,.true.)/vmin(i)/3.d0
c          kstrt = i-KITRP/2
c          kstrt = max(1, kstrt)
c          kstrt = min(nt-KITRP+1, kstrt)
c          call interp(t(kstrt),vmin(kstrt),
c     &      KITRP,t(i),dum,vderv,.true.)
c          vderv = vderv/vmin(i)/3.d0
          
	  Cp(i) = DQDDER(1,t(i),nt,t,hmin,.true.)
	    tmpgmu = 0.0
	    tmpexel = 0.0
	    tmpdVdT = 0.0
	    tmpdMudT = 0.0
          if (nv.gt.1) then
	    call DCSIEZ (nv, v, cvion(1,i), 1, vmin(i), Cv_ion)
	    call DCSIEZ (nv, v, sion(1,i), 1, vmin(i), S_ion)
	    call DCSIEZ (nv, v, selec(1,i), 1, vmin(i), S_elec)
	    call DCSIEZ (nv, v, Emag(1,i), 1, vmin(i), E_spin)
	    call DCSIEZ (nv, v, Smag(1,i), 1, vmin(i), S_spin)
	    call DCSIEZ (nv, v, Cmag(1,i), 1, vmin(i), Cv_spin)
	    call DCSIEZ (nv, v, cvelec(1,i), 1, vmin(i), Cv_ele)
	    call DCSIEZ (nv, v, gmu(1,i), 1, vmin(i), tmpgmu)
	    call DCSIEZ (nv, v, exel(1,i), 1, vmin(i), tmpexel)
	    call DCSIEZ (nv, v, dVdT(1,i), 1, vmin(i), tmpdVdT)
	    call DCSIEZ (nv, v, dMudT(1,i), 1, vmin(i), tmpdMudT)
	    tmpxgmu = DQDDER(1,t(i),nt,t,xgmu,.true.)
	    tmpexel = DQDDER(1,vmin(i),nv,v,gmu(1,i),.true.)
	  else
	    Cv_ion = cvion(1,i)
	    S_ion = sion(1,i)
            S_elec = selec(1,i)
            E_spin = Emag(1,i) 
            S_spin = Smag(1,i) 
            Cv_spin = Cmag(1,i) 
	    Cv_ele = cvelec(1,i)
	    tmpgmu = gmu(1,i)
	    tmpexel = 0.0
	    tmpdVdT = dVdT(1,i)
	    tmpdMudT = dMudT(1,i)
	    tmpxgmu = 0.0
	  endif
	  Cp_ion = ymDDER(t(i),nt,t,hion_min, 1)
	  Cv(i) = Cv_ion + Cv_ele + Cv_spin
c          write (0,*) "Cp(i)=",  Cv(i), Cv_ion, Cv_ele, Cv_spin
c      write(0,*) "Cv_spin=", Cv_spin

          if (t(i).eq.0.d0) then
            if (nv.eq.1) then
              dummy = f00(1)
            else
	      call DCSIEZ (nv, v, f00, 1, vmin(i), dummy)
            endif
            tdeb = dummy*8.d0/9.d0
     &        /dCONST('Boltzman')*dconst('ElectronVolt')
          else
            tdeb = debyet(Cv_ion, t(i))
          endif

	  if (nv.le.2) then
	    Bt = 0.d0
	  else
	    call DCSIEZ (nv, v, f(1,i),NRUNX, runv, runf)
            CALL ERSET (3, 0, -1)
            if (modeBt.eq.0) then
	      Bt = vmin(i)*DQDDER(2, vmin(i),NRUNX, runv, runf,.true.)
            else
	      Bt = vmin(i)*ymBDER(nv, v, f(1,i), vmin(i), 4)
            endif
            CALL ERSET (3, 2,  2)
	  endif
          bmin(i) = Bt

	  dCp = (vderv*3.d0)**2*Bt*t(i)*vmin(i)
	  if (Cv(i).gt.0.d0) then
            ggam = vderv*3.d0*Bt*vmin(i)/Cv(i)
            alpha(i) = vderv*3.d0
	  else
            alpha(i) = 0.d0
	    ggam = 0.d0
	  endif
	  if (abs(ggam).ge.100.d0) ggam = 0.0d0

c DO NOT USE the derivative for Cp
c          write (0,*) "Cp(i)=",  Cp(i)

          Cp(i) = Cv(i) + dCp
c          write (0,*) "Cp(i)=",  Cp(i), Cv(i), dCp
          if ( t(i).eq.0.0 ) then
             smin(i) = 0.0
             Cp(i) = 0.0
             Cv_ion = 0.0
             Cp_ion = 0.0
             S_ion = 0.0
             vderv = 0.0
             Cv(i) = 0.0
             dCp = 0.0
          endif

	  write (6, '(f6.1, (f11.6,f14.8,f9.4, f15.1, 
     &              f13.8, 1p3e14.6,
     &              0p3f9.4, 
     &              f15.1, f11.4, 
     &              1p3e14.6, 0pf9.4, f9.4,
     &              1p5e14.6))')
     &        t(i), vmin(i), fmin(i), smin(i)*UtoJ, hmin(i)*UtoJ,
     &        vderv*1.e6, Cp(i)*UtoJ, Cv_ion*UtoJ, Cp_ion*UtoJ, 
     &        S_ion*UtoJ, S_elec*UtoJ, S_spin*UtoJ,
     &        hion_min(i)*UtoJ, tdeb,  
     &        Cv(i)*UtoJ, Bt*UtoGPa, dCp*UtoJ, ggam, E_spin*UtoJ,
     &        tmpgmu, tmpexel*vderv*3.0*vmin(i)*1.e6, 
     &        tmpdVdT, tmpdMudT, tmpxgmu*1.e6
	enddo
      open (7, file=mode(1:len_trim(mode))//"_svib")
      open (8, file=mode(1:len_trim(mode))//"_fvib")
      open (9, file=mode(1:len_trim(mode))//"_fvib_plt")
      write (8, "(i6,i6,99f14.8)") nv, nt, (v(j), j=1,nv)
      do i = 1,nt
        write (7, "(f6.1,99f11.4)") t(i), (s(j,i)*UtoJ, j=1,nv)
        if (t(i).ne.0.d0) then
          deltaB = vmin(i)*t(i)*(alpha(i)*bmin(i))**2/Cv(i)
          CpCv = Cp(i)/Cv(i)
        else
          deltaB = 0.d0
          CpCv = 1.d0
        endif
        write (8, "((1p99e20.10))") t(i), 
     &    vmin(i), fmin(i), bmin(i), deltaB, CpCv,
     &    (f(j,i), s(j,i), Ezero(j), cvion(j,i)+cvelec(j,i), j=1,nv)
        if (nv.gt.1) then
	  call DCSIEZ (nv, v, f(1,i),NRUNX, runv, runf)
	  call DCSIEZ (nv, v, felec(1,i),NRUNX, runv, elef)
	  do j=1, NRUNX
	    runp(j) = -DQDDER(1, runv(j),NRUNX, runv, runf,.true.)
	  enddo
          write (9, "(f6.1,5f14.8)") 
     &      (t(i), runv(j), runf(j), elef(j), 
     &      runp(j)*UtoGPa, runf(j)+runv(j)*runp(j), 
     &      j=1,NRUNX,NRUNX/100)
        else
          write (9, "(f6.1,3f14.8)") 
     &      (t(i), v(j), f(j,i), felec(j,i), j=1,nv)
        endif
        write (9, "(/)")
      enddo
      close (7)
      close (8)
      close (9)
      if (nv.lt.4) icalat = 0
      if (icalat.eq.1) then
        do i = 1,nt
	  do j = 1, 9
	    call DCSIEZ (nv, v, veclat(1,j), 1, vmin(i), latva(i,j))
          enddo
        enddo
        open (7, file=trim(mode)//"_eij_dT")
        write (7, '("# T(K)",6("  e",i1,"(-6/K)"))')
     &    (i,i=1,6)
        open (9, file=trim(mode)//"_latij_dT")
        open (10, file=trim(mode)//"_abc_dT")

c        write (9, '("# T(K)",9("    a(",i1,",",i1,")"))')
c     &    ((i,j,j=1,3),i=1,3)
        write (9, '("# T(K)",9("     a(",i1,",",i1,")",
     &    " da/dT(-4/K)"))')
     &    ((i,j,j=1,3),i=1,3)
        write (10, '("# T(K)",3("          ",a1,
     &    " d",a1,"/dT(-6/K)"))') 'a', 'a', 'b', 'b', 'c', 'c'

        if (trim(cijfile).ne."") then
          open (8, file=trim(mode)//"_Cij")
          write (8, '("# T(K)",21("    C",i1,i1,"(T)", 
     &                            "    C",i1,i1,"(S)"), 
     &      "     Bv(T)     Gv(T)     Ev(T)")')
     &      ((i,j,i,j,j=i,6),i=1,6)
        endif
        do i = 1,nt
	  do j = 1, 9
	    tmpv(j) = DQDDER(1,t(i),nt,t,latva(1,j),.true.)
            revlat(j) = tmpv(j)
c            if (abs(latva(i,j)).ge.0.00001d0) then
c            else
c               revlat(j) = 0.d0
c            endif
          enddo
c          write (9, "(f6.1,9(f10.6))") t(i), 
c     &      (latva(i,j), j=1,9)
          write (9, "(f6.1,9(f11.6,1x,f11.6))") t(i),
     &      (latva(i,j), revlat(j)*1.e4, j=1,9)
          ta = latva(i,1)**2 + latva(i,2)**2 + latva(i,3)**2
          tb = latva(i,4)**2 + latva(i,5)**2 + latva(i,6)**2
          tc = latva(i,7)**2 + latva(i,8)**2 + latva(i,9)**2
          tat = (tmpv(1)*latva(i,1) + tmpv(2)*latva(i,2)
     &         + tmpv(3)*latva(i,3))/ta
          tbt = (tmpv(4)*latva(i,4) + tmpv(5)*latva(i,5)
     &         + tmpv(6)*latva(i,6))/tb
          tct = (tmpv(7)*latva(i,7) + tmpv(8)*latva(i,8)
     &         + tmpv(9)*latva(i,9))/tc
          write (10, "(f6.1,3(f11.6,1x,f11.6))") t(i),
     &      sqrt(ta),tat*1.e6, sqrt(tb),tbt*1.e6, sqrt(tc),tct*1.e6 

	  do j = 1, 9
            revlat(j) = latva(i,j)
          enddo
          call inv_m(revlat, revm)
          call mproduct(tmpv,revm,revlat)

c          call mproduct(revm,tmpv,revlat)
c          write (0, '("tmpv =", 9f8.3)') (tmpv(j)*1.e6,j=1,9)
c          write (0, '("revm =", 9f8.2)') (revm(j),j=1,9)
c          write (0, '("revl =", 9f8.3)') (revlat(j),j=1,9)
c          write (0, '("latv =", 9f8.3)') (latva(i,j),j=1,9)
          write (7, "(f6.1,6(f10.4))") t(i), 
     &        revlat(1)*1.e6, revlat(5)*1.e6, revlat(9)*1.e6,
     &       (revlat(2)+revlat(4))*0.5e6,
     &       (revlat(3)+revlat(7))*0.5e6,
     &       (revlat(6)+revlat(8))*0.5e6

          if (trim(cijfile).ne."") then
	    do j = 1, 36
	      call DCSIEZ (nv, v, vecCij(1,j), 1, vmin(i), tmpCij(j))
            enddo
            tmpv(1) = revlat(1)
            tmpv(2) = revlat(5)
            tmpv(3) = revlat(9)
            tmpv(6) = revlat(2)+revlat(4)
            tmpv(5) = revlat(3)+revlat(7)
            tmpv(4) = revlat(6)+revlat(8)
            do j=1,6
              tmp = 0.d0
              do k=1,6
                tmp = tmp + tmpv(k)*tmpCij((j-1)*6+k)
              enddo
              revm6(j) = tmp/UtoGPa
            enddo
            do j=1,6
              do k=1,6
                tmp = 0.d0
                if (Cv(i).ge.1.d-8)
     &            tmp = tmp + t(i)*vmin(i)*revm6(j)*revm6(k)/Cv(i)
                aseCij((j-1)*6+k) = tmpCij((j-1)*6+k)+tmp*UtoGPa
              enddo
            enddo
            C11 = tmpCij((1-1)*6+1)
            C22 = tmpCij((2-1)*6+2)
            C33 = tmpCij((3-1)*6+3)
            C12 = tmpCij((1-1)*6+2)
            C13 = tmpCij((1-1)*6+3)
            C23 = tmpCij((2-1)*6+3)
            C44 = tmpCij((4-1)*6+4)
            C55 = tmpCij((5-1)*6+5)
            C66 = tmpCij((6-1)*6+6)
            A = (C11 + C22 + C33)/3.d0
            B = (C12 + C13 + C23)/3.d0
            C = (C44 + C55 + C66)/3.d0
            Bv = (A + 2.d0*B)/3.d0
            Gv = (A - B + 3.d0*C)/5.d0
            Ev = 9.d0*Bv*Gv/(Gv+3.d0*Bv)
            write (8, "(f6.1,21(2(1x,f9.2)), 3(1x,f9.2))") t(i), 
     &        ((tmpCij((j-1)*6+k), aseCij((j-1)*6+k),k=j,6), j=1,6),
     &        Bv, Gv, Ev
          endif
        enddo
        close (7)
        close (9)
        close (10)
        if (trim(cijfile).ne."") close(8)
      endif

      open (7, file=mode(1:len_trim(mode))//"_Mix", form="unformatted")
      do i=1,nt
        hmin(i) = hmin(i) -Pressure*vmin(i)
        do j=1,nv
          f(j, i) = f(j,i) -Pressure*v(j)
        enddo
      enddo
      if (ihugoniot.eq.1)
     &  call Hugoniot(nt, nv, t, vmin, hmin, v, f, s, mode)
      write (7) nv, nt
      do i=1,nv
        do j=1,nt
          cvion(i,j) = cvion(i,j) + cvelec(i,j)
        enddo
      enddo
      write (7) v, t, vmin, f, s, cvion, sion, cvelec, selec
      close (7)
      call cputim ("Cost is:")
      end
c
      real*8 function Hug(x)
c
      implicit real*8(a-h,o-z)
      real*8, pointer:: vv(:), pp(:), ee(:)
      common/comchg/v0, p0, e0, vv, pp, ee
      common/comnnv/nnv
      call DCSIEZ (nnv, vv, pp, 1, x, p)
      call DCSIEZ (nnv, vv, ee, 1, x, e)
      dp = (p+p0)*(v0-x)
      de = e - e0
c      Hug = dp - 0.5e0*de
      Hug = 0.5e0*dp - de
c      Hug = dp - de
      end
c
      subroutine Hugoniot(nt, nv, t, vmin, emin, v, f, s, mode)
c
      implicit real*8(a-h,o-z)
      character*(*) mode
      real*8, allocatable, target:: 
     &  pH(:),  vH(:), fH(:), sH(:), eH(:), 
     &  vT(:), fT(:), pT(:), eT(:)
      real*8, pointer:: vv(:), pp(:), ee(:)
      common/comchg/v0, p0, e0, vv, pp, ee
      common/comnnv/nnv
      dimension t(nt), vmin(nt), emin(nt), v(nv), f(nv, nt), s(nv,nt)
      external Hug
      NRUNX = 101

      UtoGPa = dconst('ElectronVolt')*1.0d21
      allocate (pH(nt), vH(nt), fH(nt), sH(nt), eH(nt), 
     &  vT(NRUNX), fT(NRUNX), pT(NRUNX), eT(NRUNX), STAT=istat)

      do i=1, nt
        pH(i) = 0.e0;
        vH(i) = vmin(i);
        sH(i) = 0.e0;
        eH(i) = emin(i);
      enddo
      dv =(v(nv)-v(1))/dble(NRUNX-1)
      do j=1,NRUNX
        vT(j) =v(1) + dble(j-1)*dv
      enddo
      nnv = NRUNX

      p0 = 0.e0
      t0 = 300.e0
      call DCSIEZ (nt, t, vmin, 1, t0, v0)
      call DCSIEZ (nt, t, emin, 1, t0, e0)
c      write (0,*) "t0,e0", t0, e0
c
      vv => vT
      pp => pT
      ee => eT
      do i=1, nt
c	 write (0,*) t(i), -ymDDER(vmin(i), nv, v, f(1, i), 1)*UtoGPa
        if (t(i).gt.t0) then
c        write (0,*) "Hug called", t(i), nt
          do j=1,nv
            pT(j) = f(j,i) + t(i)*s(j,i)
          enddo
         call DCSIEZ (nv, v, pT, nnv, vv, eT)
         call DCSIEZ (nv, v, f(1,i), nnv, vv, fT)
         do j=1, nnv
	    pT(j) = -ymDDER(vT(j), nnv, vv, fT, 1)
c	    pT(j) = -DQDDER(1,v(j),nv,v,f(1,i),.true.)
          enddo
          A = vv(1)
          vA = Hug(A)
          do j=2, nnv
            B = vv(j)
            vB = Hug(B)
            if (vA*vB.le.zero) goto 100
            A = B
            vA = vB
          enddo
          goto 200

 100      ERRABS = 0.d0
          ERRREL = 1.d-12
          MAXFN = 999
c          write (0,*) "A,vA", A, vA, B, vB
          CALL DZBREN (Hug, ERRABS, ERRREL, A, B, MAXFN)
c          do k=1,32
c            C= 0.5d0*(A+B)
c            vC = Hug(C)
c            if (vC*vA.gt.0.e0) then
c              A = C
c              vA = vC
c            else
c              B = C
c              vB = vC
c            endif
c          enddo
          vH(i) = B
          call DCSIEZ (nnv, vv, pT, 1, B, pH(i))
          call DCSIEZ (nnv, vv, eT, 1, B, eH(i))
          if (B.ge.vv(nnv)*0.99e0) goto 300
        write (0, '(f8.1, f11.6,f14.8,f15.4)')
     &        t(i), vH(i), eH(i), pH(i)*UtoGPa
c        write (0, *) (eH(i) - e0) - 0.5e0*(pH(i) - p0) * (v0- vH(i))
c     &  , Hug(vH(i))
 200      continue
        endif
      enddo
 300  continue
      open (11, file=mode(1:len_trim(mode))//"_Hug")
      do i=1,nt
        if (vH(i).ge.vv(nnv)*0.99e0) goto 400
        write (11, '(f8.1, f11.6,f14.8,f15.4)')
     &        t(i), vH(i), eH(i), pH(i)*UtoGPa
      enddo
 400  continue
      end
c
      subroutine settim
c
      implicit real*8(a-h,o-z)
      character*(*) str
      real*4 tm(2),t0,t1,t00
      data zero/0.0/
      save t0
c      t00 = secnds(zero)
      call etime(tm, t00)
      t0 = tm(1) + tm(2)
      return
c
      entry cputim (str)
c      t1 = secnds(zero)
      call etime(tm, t00)
      t1 = tm(1) + tm(2)
      write (0, 100) t1-t0, t1-t00, str
      t0 = t1
 100  format (" Partial Time: ", f9.2," Secs. ", f9.2," (Total) in ", a)
      end
c
      subroutine sum_d_free_energy(nt, t, fndos, fdum, sdum, 
     &  cdum, natom, unitf)
c
      implicit real*8(a-h,o-z)
      common/lowtem/toobig
      character*1024 fndos
      dimension t(nt), fdum(nt), sdum(nt), cdum(nt)
      real*8, allocatable:: rdos(:), fdos(:)

      open (7, file=fndos)
      n = 0
      do while (.TRUE.)
        read (7, *, end=100) freq, dum
        n = n + 1
      enddo
100   allocate (rdos(n), fdos(n), STAT=istat)
      if (istat.ne.0) then 
        print *, "Cannot allocate dynamic memory for fdos"
        stop
      endif 

      rewind(7)
      do i = 1, n
        read (7, *) freq, dum
        rdos(i) = freq*unitf
        fdos(i) = dum/unitf
      enddo
200   close(7)

      w0 = 0.d0
      do i = 1, n
        w0 = w0 + fdos(i)
      enddo
c      write (0,*), "w0 = ", w0
      pp0 = 3.d0*dble(natom)/w0
      do i = 1, n
        fdos(i) = fdos(i)*pp0
      enddo
      w0 = 0.d0
      do i = 1, n
        w0 = w0 + fdos(i)
      enddo
c      write (0,*), "new w0 = ", w0

      h = dconst('Planck')
      do j = 1, nt
        temper = t(j)
	tkb = temper*dCONST('Boltzman')
	f0 = 0.d0
	s0 = 0.d0
	c0 = 0.d0
	  do i = 1, n
            ff = fdos(i)
	    hmu = h*rdos(i)
	    if (hmu.gt.0.d0.and.temper.gt.0.d0) then
	      x = hmu/tkb
	      xd2 = 0.5d0*x 
              if (x.ge.toobig) then
	        f0 = f0 + 0.5d0*hmu*ff
              else
		f0 = f0 + tkb*log(2.d0*sinh(xd2))*ff
		s0 = s0 + (x/(exp(x)-1.d0) - log(1.d0 - exp(-x)))*ff
		c0 = c0 + (xd2/sinh(xd2))**2*ff
              endif
	    else if (hmu.gt.0.d0.and.temper.eq.0.d0) then
	        f0 = f0 + 0.5d0*hmu*ff
	    endif
	  enddo
        fdum(j) = f0/dconst('ElectronVolt')
        sdum(j) = s0*dCONST('Boltzman')/dconst('ElectronVolt')
        cdum(j) = c0*dCONST('Boltzman')/dconst('ElectronVolt')
      enddo
      deallocate(fdos, rdos)
      end
c
      subroutine int_d_free_energy(nt, t, fndos, 
     &           fdum, sdum, cdum, natom, imode, unitf, izero)
c
      implicit real*8(a-h,o-z)
      common/rundis/NRUNX
      common/lowtem/toobig
      character*1024 fndos
      dimension t(nt), fdum(nt), sdum(nt), cdum(nt)
      real*8, allocatable:: fdos(:), rdos(:),
     &  rx(:), ados(:), a(:), fx(:), sx(:), cx(:) 
c      NRUNX = 10001

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
        read (7, *) freq, dum
        rdos(i) = freq*unitf
        fdos(i) = dum/unitf
c        read (7, *) rdos(i), fdos(i)  
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

c      azero = ados(1)
c      nM = NRUNX/100
c      nX = 1
c      do i=1, nM
c        if (ados(i).lt.azero) then
c          nX = i
c          azero = ados(i)
c        endif
c      enddo
c      if (izero.eq.0) then
c        do i=1, nX
c          ados(i) = azero
c        enddo
c      endif
c      open (7, file=trim(fndos)//"_Z")
c      do i = 1, NRUNX
c        if (izero.eq.0) then
c	  ados(i) = ados(i) - azero
c          if (ados(i).le.0.e0) ados(i) = 0.e0
c          write (7,*) rx(i), ados(i)
c        endif
c      enddo
c      close(7)

      do i = 1, NRUNX
        fx(i) = ados(i)
        if (ados(i).le.0.e0) ados(i) = 0.e0
      enddo

      call integr (ados, rx, NRUNX, a, 3)
      w0 = 3.d0*dble(natom)/a(NRUNX)
      write (0, '(a, f11.6, a, f11.6)') "******** Sum DOS = ",
     &    a(NRUNX), ",      Ideal DOS = ", 3.d0*dble(natom)
      if (imode.eq.1) then
        do i = 1, NRUNX
          ados(i) = ados(i)*w0
        enddo
      endif
      if (izero.ne.0) then
        call integr (fx, rx, NRUNX, a, 3)
        w0 = 3.d0*dble(natom)/a(NRUNX)
        write (0, '(a)') "******** Imaginary phonon corrected by -izero"
        do i = 1, NRUNX
          ados(i) = fx(i)*w0
        enddo
      endif

      call integr (ados, rx, NRUNX, a, 3)
      if (imode.eq.1) then
        write (0, '(a, f11.6)') 
     &    "******** renormalize to Ideal DOS ", a(NRUNX)
      endif

      h = dconst('Planck')
      do j = 1, nt
        temper = t(j)
	  tkb = temper*dCONST('Boltzman')
	  do i = 1, NRUNX
	    hmu = h*rx(i)
            if (izero.ne.0.and.hmu.le.0.d0) hmu = -hmu
	    fx(i) = 0.d0
	    sx(i) = 0.d0
	    cx(i) = 0.d0
	    if (hmu.gt.0.d0.and.temper.gt.0.d0) then
	      x = hmu/tkb
	      xd2 = 0.5d0*x 
              if (x.ge.toobig) then
	        fx(i) = fx(i) + 0.5d0*hmu*ados(i)
              else
		fx(i) = tkb*log(2.d0*sinh(xd2))*ados(i)
		sx(i) = (x/(exp(x)-1.d0) - log(1.d0 - exp(-x)))*ados(i)
		cx(i) = (xd2/sinh(xd2))**2*ados(i)
              endif
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
      subroutine f_free_energy(nt, t, fnfrq, fdum, sdum, cdum)
c
      implicit real*8(a-h,o-z)
      common/lowtem/toobig
      character*1024 fnfrq, line
      dimension t(nt), fdum(nt), sdum(nt), cdum(nt)
      open (7, file=fnfrq)
      line = ""
      do while (line.ne."Phonon frequencies:")
        read (7, '(a)', end=100) line
      enddo

      do i = 1, nt
        fdum(i) = 0.d0
        sdum(i) = 0.d0
        cdum(i) = 0.d0
      enddo
      n = 0
      read (7, '(a)', end=100) line
      Boltzman = dCONST('Boltzman')
      Planck = dconst('Planck')
      do while (line(1:3).ne."end")
        read (line, *) freq
        if (freq.gt.0.d0) then
          n = n + 1
          t0 = freq*Planck
          t05 = 0.5d0*t0
          do i = 1, nt
              temper = t(i)
              if (temper.gt.0.d0) then
                tkb = temper*Boltzman
                x = t0/tkb
                x05 = 0.5d0*x
                if (x.ge.toobig) then
                  fdum(i) = fdum(i) + t05
                else
                  sinh05 = sinh(x05)
                  c05 = x05/sinh05
                  fdum(i) = fdum(i) + tkb*log(2.d0*sinh05)
                  sdum(i) = sdum(i) + 
     &              (x/(exp(x)-1.d0) - log(1.d0 - exp(-x)))
                  cdum(i) = cdum(i) + c05*c05
                endif
              else
                fdum(i) = fdum(i) + t05
              endif
          enddo
        endif
        read (7, '(a)', end=100) line
      enddo
100   close(7)
      ftoeV = 1.d0/dconst('ElectronVolt')/dble(n)
      ctoeV = dCONST('Boltzman')/dconst('ElectronVolt')/dble(n)
      do i = 1, nt
        fdum(i) = 3.d0*fdum(i)*ftoeV
        sdum(i) = 3.d0*sdum(i)*ctoeV
        cdum(i) = 3.d0*cdum(i)*ctoeV
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

c      CALL DRCURV (nv, XDATA, f, NDEG, B, SSPOLY, STAT)
      CALL Least_Square(nv, NDEG, XDATA, f, B)
      x0 = vol0**(-2.d0/3.d0)
      x1 = (-2.d0/3.d0)*x0/vol0
      x2 = (-2.d0/3.d0)*(-5.d0/3.d0)*x0/vol0/vol0
      t1 = 0.d0
      t2 = 0.d0
      do i = 1, NDEG
        t1 = t1 + dble(i)*B(i+1)*x0**(i-1)
        t2 = t2 + dble(i-1)*dble(i)*B(i+1)*x0**(i-2)
      enddo
c      t1 = DQDDER(1, x0, nv, XDATA, f, .true.)
c      t2 = DQDDER(2, x0, nv, XDATA, f, .true.)
c      d1 = t1*x1
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
c	do B = A, 100000.1d0, 10.d0
        B = A
        do while ( B < 100000.1d0 )
	  vB = fdebye(B/temper)
	  if (vA*vB.le.zero) goto 100
	  A = B
	  vA = vB
          B = B + 10.d0
	enddo
        cdeb = (cv*dconst('ElectronVolt')/dconst('Boltzman'))
        debyet = (2.4e0*dconst("PI")**4/cdeb)**(1./3.)*temper
c	debyet = zero
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

	if (1.eq.1.or.vA*vB.gt.0.d0) then
          call interp(p,x,KITRP,zero,xzero,dum,.false.)
	  if (xzero.lt.x(1).or.xzero.gt.x(KITRP)) xzero=x(KITRP/2+1)
	else
	  ERRABS = 0.d0
	  ERRREL = 1.d-12
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

      subroutine mproduct(a,b,c)
      implicit real*8(a-h,o-z)
      dimension a(3,3), b(3,3), c(3,3)
      do i=1,3
        do j=1,3
          t = 0.d0
          do k=1,3
            t = t + a(i,k)*b(k,j)
          enddo
          c(i,j) = t
        enddo
      enddo
      end

      subroutine inv_m(m, mout)
      implicit real*8(a-h,o-z)
      real*8 m(0:2,0:2), mout(0:2,0:2)
      d = volume(m);
      mout(0,0)=(+m(1,1)*m(2,2)-m(1,2)*m(2,1))/d;
      mout(0,1)=(-m(0,1)*m(2,2)+m(0,2)*m(2,1))/d;
      mout(0,2)=(+m(0,1)*m(1,2)-m(0,2)*m(1,1))/d;
      mout(1,0)=(-m(1,0)*m(2,2)+m(1,2)*m(2,0))/d;
      mout(1,1)=(+m(0,0)*m(2,2)-m(0,2)*m(2,0))/d;
      mout(1,2)=(-m(0,0)*m(1,2)+m(0,2)*m(1,0))/d;
      mout(2,0)=(+m(1,0)*m(2,1)-m(1,1)*m(2,0))/d;
      mout(2,1)=(-m(0,0)*m(2,1)+m(0,1)*m(2,0))/d;
      mout(2,2)=(+m(0,0)*m(1,1)-m(0,1)*m(1,0))/d;
      end

      real*8 function volume(a)
      implicit real*8(a-h,o-z)
      dimension a(0:2,0:2)
      volume = a(0,0)*a(1,1)*a(2,2)
     &       + a(0,1)*a(1,2)*a(2,0)
     &       + a(0,2)*a(1,0)*a(2,1)
     &       - a(0,2)*a(1,1)*a(2,0)
     &       - a(0,0)*a(1,2)*a(2,1)
     &       - a(0,1)*a(1,0)*a(2,2)
      end


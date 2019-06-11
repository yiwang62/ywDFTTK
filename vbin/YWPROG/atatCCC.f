      program atat_of_phonon_elastic_constant
c      USE NUMERICAL_LIBRARIES
c
c -------- energy in Ryd, temperature in K, coordinate in a.u.
c
      implicit real*8(a-h,o-z)

      real*8, allocatable, target:: 
     &  t(:),  vmin(:), fmin(:), bmin(:),
     &  v(:), f(:), vcld(:), ecld(:),
     &  CCC(:, :), CCCcld(:, :), nstrn(:)
      character*1024 fn, fcld, strain_type*16, sgrp*8
      parameter (NSMOOTH = 10001)
      dimension vs(NSMOOTH), es(NSMOOTH)
c
      UtoGPa = dconst('ElectronVolt')*1.0d21
      fcld = " "
      narg = iargc()
      i = 1
      itst = 0
      do while (i.le.narg)
        call getarg(i, fn)
        if (fn(1:5).eq."-cold") then
          i = i + 1
          call getarg(i, fcld)
        else if (fn(1:5).eq."-test") then
          i = i + 1
          call getarg(i, fcld)
          i = i + 1
          call getarg(i, fn)
          read (fn, *) nt
          nv = nt
          itst = 1
        else if (fn(1:5).eq."-cub1") then
          sgrp = "cub1"
          nEC = 2
        else if (fn(1:5).eq."-cub2") then
          sgrp = "cub2"
          nEC = 2
        else if (fn(1:5).eq."-hex1") then
          sgrp = "hex1"
          nEC = 4
        else if (fn(1:5).eq."-hex2") then
          sgrp = "hex2"
          nEC = 4
        else if (fn(1:5).eq."-tet1") then
          sgrp = "tet1"
          nEC = 6
        else if (fn(1:5).eq."-ort1") then
          sgrp = "ort1"
          nEC = 10
        endif
        i = i + 1
      enddo

      IF (itst.eq.0) THEN
        read (5, *) fn
        open (7, file=fn)
        read (7, *) nv, nt
        close (7)
      ENDIF

      allocate (
     &  t(nt),  vmin(nt), fmin(nt), bmin(nt),
     &  v(nv), f(nv), vcld(nv), ecld(nv),
     &  CCC(nEC, nt), CCCcld(nv, nEC), nstrn(nEC),
     &  STAT=istat)
      if (istat.ne.0) then 
        print *, "Cannot allocate main dynamic memory"
        stop
      endif 

      do j = 1, nEC
        nstrn(j) = 0
        do i = 1, nt
          CCC(j, i) = 0.d0
        enddo
      enddo

      IF (itst.eq.0) THEN
        open (7, file=fn)
        read (7, *) i1, i2, v
        do i = 1, nt
          read (7,*) t(i), vmin(i), fmin(i), bmin(i)
        enddo
        close (7)
      ENDIF
c
c     semi phonon calculation
c
      if (fcld.eq." ") goto 1001
      do j = 1, nEC
        nstrn(j) = 1
      enddo
      open (7, file=fcld)
      do i = 1, nv
        read (7,*) vcld(i), ecld(i), (CCCcld(i,j), j=1, nEC)
        do j = 1, nEC
          CCCcld(i,j) = CCCcld(i,j)/UtoGPa
        enddo
      enddo
      close (7)
     
      IF (itst.eq.1) THEN
        do i = 1, nv
          t(i) = 0
          vmin(i) = vcld(i)
          v(i) = vcld(i)
        enddo
      ENDIF
 
      dv = (v(nv) - v(1))/dble(NSMOOTH-1)
      do i = 1, NSMOOTH
        vs(i) = v(1) + dble(i-1)*dv
      enddo
      call DCSIEZ (nv, v, ecld, NSMOOTH, vs, es)
      do i = 1,nt
        bmin(i) = vmin(i)*DQDDER(2, vmin(i),NSMOOTH,vs,es,.true.)
        do j = 1, nEC
          call DCSIEZ (nv, v, CCCcld(1,j), 1, vmin(i), CCC(j,i))
        enddo
      enddo
      goto 100
c
c     full phonon calculation
c
1001  continue
 10   read (5, *, end=100) fn, strain_type, ee
      if (sgrp(1:3).eq."cub") then
        if (index(strain_type,"Cp").gt.0) then
          key = 1
        else if (index(strain_type,"C44").gt.0) then
          key = 2
        endif
      else if (sgrp.eq."hex1") then
        if (index(strain_type,"Cp").gt.0) then
          key = 1
        else if (index(strain_type,"Cm").gt.0) then
          key = 2
        else if (index(strain_type,"C33").gt.0) then
          key = 3
        else if (index(strain_type,"C55").gt.0) then
          key = 4
        endif
      else if (sgrp.eq."hex2") then
        if (index(strain_type,"Cs").gt.0) then
          key = 1
        else if (index(strain_type,"C44").gt.0) then
          key = 2
        else if (index(strain_type,"C66").gt.0) then
          key = 3
        else if (index(strain_type,"R").gt.0) then
          key = 4
        endif
      else if (sgrp.eq."tet1") then
        if (index(strain_type,"C11").gt.0) then
          key = 1
        else if (index(strain_type,"C33").gt.0) then
          key = 2
        else if (index(strain_type,"C44").gt.0) then
          key = 3
        else if (index(strain_type,"Cxx").gt.0) then
          key = 4
        else if (index(strain_type,"Cyy").gt.0) then
          key = 5
        else if (index(strain_type,"Czz").gt.0) then
          key = 6
        endif
      else if (sgrp.eq."ort1") then
        if (index(strain_type,"C11").gt.0) then
          key = 1
        else if (index(strain_type,"C22").gt.0) then
          key = 2
        else if (index(strain_type,"C33").gt.0) then
          key = 3
        else if (index(strain_type,"C44").gt.0) then
          key = 4
        else if (index(strain_type,"C55").gt.0) then
          key = 5
        else if (index(strain_type,"C66").gt.0) then
          key = 6
        else if (index(strain_type,"Cxx").gt.0) then
          key = 7
        else if (index(strain_type,"Cyy").gt.0) then
          key = 8
        else if (index(strain_type,"Czz").gt.0) then
          key = 9
        endif
      endif

      nstrn(key) = nstrn(key) + 1

      open (7, file=fn)
      read (7, *) 
      do i = 1, nt
        read (7,*) d1,d2,d3,d4, f
        call DCSIEZ (nv, v, f, 1, vmin(i), fx)
        value = (fx - fmin(i))/vmin(i)/ee**2
        CCC(key, i) = CCC(key, i) + value
      enddo
      close (7)
      goto 10
 100  if (sgrp.eq."cub1") then
        write (6, '("# T V B C11 C12 C44")')
        do i = 1,nt
          Cpm = CCC(1, i)/nstrn(1)/6.d0
          C44 = CCC(2, i)/nstrn(2)/2.d0
          C11 = bmin(i) + 4.d0/3.d0*Cpm
          C12 = bmin(i) - 2.d0/3.d0*Cpm
          write (6, '(f7.1, f11.6, 4f9.3)') t(i), vmin(i),
     &      bmin(i)*UtoGPa, C11*UtoGPa, C12*UtoGPa,C44*UtoGPa
        enddo
      else if (sgrp.eq."cub2") then
        write (6, '("# T V B C11 C12 C44")')
        do i = 1,nt
          Cpm = CCC(1, i)/nstrn(1)/6.d0
          C44 = CCC(2, i)/nstrn(2)/6.d0
          C11 = bmin(i) + 4.d0/3.d0*Cpm
          C12 = bmin(i) - 2.d0/3.d0*Cpm
          write (6, '(f7.1, f11.6, 4f9.3)') t(i), vmin(i),
     &      bmin(i)*UtoGPa, C11*UtoGPa, C12*UtoGPa,C44*UtoGPa
        enddo
      else if (sgrp.eq."hex1") then
        write (6, '("# T V B C11 C12 C13 C33 C55")')
        do i = 1,nt
          Cp = CCC(1,i)/nstrn(1)
          Cm = CCC(2,i)/nstrn(2)
          C33 = CCC(3,i)/nstrn(3)
          C55 = CCC(4,i)/nstrn(4)
          C11 = 0.5d0*(Cp + Cm)
          C12 = 0.5d0*(Cp - Cm)
          C33 = 2.0d0*C33
          C55 = 0.5d0*C55
          C13 = 0.5d0*(4.5d0*bmin(i) - C11 -C12 - 0.5d0*C33)
          write (6, '(f7.1, f11.6, 7f9.3)') t(i), vmin(i),
     &      bmin(i)*UtoGPa,
     &      C11*UtoGPa, C12*UtoGPa, C13*UtoGPa, C33*UtoGPa, C55*UtoGPa
        enddo
      else if (sgrp.eq."hex2") then
        write (6, '("# T V B C11 C12 C13 C33 C44")')
        do i = 1,nt
          Cs = CCC(1,i)/nstrn(1)
          C44 = 0.5d0*CCC(2,i)/nstrn(2)
          C66 = 0.5d0*CCC(3,i)/nstrn(3)
          R = CCC(4,i)/nstrn(4)
          x = Cs
          y = R*Cs
          d = bmin(i)*Cs
          xy = x + y
          x2y = x -2.d0*y
          xyx2y = xy*x2y
          e = (6.d0*xy + 3.d0*x2y)/9.d0
          C13 = (d - xyx2y/9.d0)/e
          C33 = xy/3.d0 + C13
          A = x2y/3.d0 + 2.d0*C13
          C11 = 0.5d0*A + C66
          C12 = 0.5d0*A - C66
          write (6, '(f7.1, f11.6, 7f9.3)') t(i), vmin(i),
     &      bmin(i)*UtoGPa,
     &      C11*UtoGPa, C12*UtoGPa, C13*UtoGPa, C33*UtoGPa, C44*UtoGPa
        enddo
      else if (sgrp.eq."tet1") then
        write (6, '("# T V B C11 C12 C13 C33 C44 C66")')
        do i = 1,nt
          C11 = CCC(1,i)/nstrn(1)
          C33 = CCC(2,i)/nstrn(2)
          C44 = CCC(3,i)/nstrn(3)
          Cxx = CCC(4,i)/nstrn(4)
          Cyy = CCC(5,i)/nstrn(5)
          Czz = CCC(6,i)/nstrn(6)
          C11 = 2.d0*C11
          C33 = 2.d0*C33
          C44 = 0.5d0*C44
          C12 = (Cyy-4.d0*Cxx)/9.d0 + C11
          C13 = 0.5d0*(5.d0*C11-4.d0*C12+C33)- Cxx
          C66 = 0.5d0*(Czz - C11 - C12 +4.d0*C13 -2.d0*C33)
          write (6, '(f7.1, f11.6, 7f9.3)') t(i), vmin(i),
     &      bmin(i)*UtoGPa,
     &      C11*UtoGPa, C12*UtoGPa, C13*UtoGPa,
     &      C33*UtoGPa, C44*UtoGPa, C66*UtoGPa
        enddo
      else if (sgrp.eq."ort1") then
        write (6, '("# T V B C11 C12 C13 C33 C44 C66")')
        do i = 1,nt
          C11 = CCC(1,i)/nstrn(1)
          C22 = CCC(2,i)/nstrn(2)
          C33 = CCC(3,i)/nstrn(3)
          C44 = CCC(4,i)/nstrn(4)
          C55 = CCC(5,i)/nstrn(5)
          C66 = CCC(6,i)/nstrn(6)
          Cxx = CCC(7,i)/nstrn(7)
          Cyy = CCC(8,i)/nstrn(8)
          Czz = CCC(9,i)/nstrn(9)
        enddo
      endif
      end

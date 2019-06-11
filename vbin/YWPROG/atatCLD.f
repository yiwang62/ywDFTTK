      program atat_of_phonon_elastic_constant
c      USE NUMERICAL_LIBRARIES
c
c -------- energy in Ryd, temperature in K, coordinate in a.u.
c
      implicit real*8(a-h,o-z)

      real*8, allocatable, target:: 
     &  t(:),  vmin(:), fmin(:), bmin(:), 
     &  vdummy(:), f(:), vcld(:), ecld(:), bcld(:),
     &  CCC(:, :), CCCcld(:, :)
      character*1024 fn, fcld, sgrp*8, str*8
      parameter (NSMOOTH = 10001)
      dimension vs(NSMOOTH), es(NSMOOTH)
c
      UtoGPa = dconst('ElectronVolt')*1.0d21
      fcld = " "
      narg = iargc()
      itst = 0
      i = 1
      do while (i.le.narg)
        call getarg(i, fn)
        if (fn.eq."-cold") then
          i = i + 1
          call getarg(i, fcld)
        else if (fn.eq."-test") then
          itst = 1
        else if (fn.eq."-v0") then
          itst = 2
          i = i + 1
          call getarg(i, str)
          read (str,*) v0K
        else if (fn.eq."-cub") then
          sgrp = "cub"
          nEC = 3
        else if (fn.eq."-hex") then
          sgrp = "hex"
          nEC = 5
        else if (fn.eq."-tet") then
          sgrp = "tet"
          nEC = 6
        else if (fn.eq."-ort") then
          sgrp = "ort"
          nEC = 9
        else if (fn.eq."-mon") then
          sgrp = "mon"
          nEC = 13
        else if (fn.eq."-tri") then
          sgrp = "tri"
          nEC = 21
        else 
          fcld = fn
        endif
        i = i + 1
      enddo

      if (fcld.eq." ") then
          print *, "********** cold file is missed"
          stop
      endif

      allocate (CCCcld(nEC,1))
      open (7, file=fcld)
      nvcld = 0
      do while (1.eq.1)
        read (7,*, end=100) t1, t2, (CCCcld(i,1), i=1, nEC)
        nvcld = nvcld + 1
      enddo
 100  close (7)

      IF (itst.eq.0) THEN
        read (5, *) nv, nt
        rewind 5
      ELSE
        nt = nvcld
      ENDIF

      allocate (
     &  t(nt),  vmin(nt), fmin(nt), bmin(nt), bcld(nt),
     &  vdummy(nv), f(nvcld), vcld(nvcld), ecld(nvcld), 
     &  CCC(nEC, nt), CCCcld(nvcld, nEC),
     &  STAT=istat)
      if (istat.ne.0) then 
        print *, "Cannot allocate main dynamic memory"
        stop
      endif 

      do j = 1, nEC
        do i = 1, nt
          CCC(j, i) = 0.d0
        enddo
      enddo
c
c     semi phonon calculation
c
      open (7, file=fcld)
      do i = 1, nvcld
        read (7,*) vcld(i), ecld(i), (CCCcld(i,j), j=1, nEC)
        do j = 1, nEC
          CCCcld(i,j) = CCCcld(i,j)/UtoGPa
        enddo
      enddo
      close (7)
     
      IF (itst.eq.0) THEN
        read (5, *) i1, i2, v
        do i = 1, nt
          read (5,*) t(i), vmin(i), fmin(i), bmin(i)
        enddo
      ELSE IF (itst.eq.1) THEN
        do i = 1, nvcld
          t(i) = 0
          vmin(i) = vcld(i)
          bmin(i) = 0.d0
          bcld(i) = 0.d0
        enddo
      ELSE IF (itst.eq.2) THEN
        tmp = (vcld(nvcld)-v0K)/dble(nvcld-1)
        do i = 1, nvcld
          t(i) = 0
          vmin(i) = v0K + dble(i-1)*tmp
        enddo
      ENDIF
 
      if (nvcld.ge.3) then
        dv = (vcld(nvcld) - vcld(1))/dble(NSMOOTH-1)
        do i = 1, NSMOOTH
          vs(i) = vcld(1) + dble(i-1)*dv
        enddo
        call DCSIEZ (nvcld, vcld, ecld, NSMOOTH, vs, es)
        do i = 1,nt
          bcld(i) = vmin(i)*DQDDER(2, vmin(i),NSMOOTH,vs,es,.true.)
          do j = 1, nEC
            call DCSIEZ (nvcld, vcld, CCCcld(1,j), 1, vmin(i), CCC(j,i))
          enddo
        enddo
      else
        do i = 1,nt
          do j = 1, nEC
            CCC(j,i) = CCCcld(i,j)
          enddo
        enddo
      endif
c
c ---- T dependent EC
c
      if (sgrp.eq."cub") then
        call output (6, "# T(K) V(Ang^3) B(GPa) C11 C12 C44 
     &      Bc(GPa) E(GPa) G(GPa) Bm(GPa)",
     &      "(a1, a6, a11, 21a9)")
      else if (sgrp.eq."hex") then
        call output (6, "# T(K) V(Ang^3) B(GPa) C11 C12 C13 C33 C44
     &      Bc(GPa) E(GPa) G(GPa) Bm(GPa)",
     &      "(a1, a6, a11, 25a9)")
      else if (sgrp.eq."tet") then
        call output (6, "# T(K) V(Ang^3) B(GPa) C11 C12 C13 C33 C44 
     &      C66 Bc(GPa) E(GPa) G(GPa) Bm(GPa)",
     &      "(a1, a6, a11, 25a9)")
      else if (sgrp.eq."ort") then
        call output (6, "# T(K) V(Ang^3) B(GPa) C11 C12 C13 C22 C23 
     &      C33 C44 C55 C66 Bc(GPa) E(GPa) G(GPa) Bm(GPa)",
     &      "(a1, a6, a11, 25a9)")
      else if (sgrp.eq."mon") then
        call output (6, "# T(K) V(Ang^3) B(GPa) C11 C12 C13 C16 C22 
     &      C23 C26 C33 C36 C44 C45 C55 C66
     &      Bc(GPa) E(GPa) G(GPa) Bm(GPa)",
     &      "(a1, a6, a11, 25a9)")
      else if (sgrp.eq."tri") then
        call output (6, "# T(K) V(Ang^3) B(GPa) C11 C12 C13 C14 C15 
     &      C16 C22 C23 C24 C25 C26 C33 C34 C35 C36 C44 C45 C46 C55 
     &      C56 C66 Bc(GPa) E(GPa) G(GPa) Bm(GPa)",
     &      "(a1, a6, a11, 25a9)")
      endif
      do i = 1,nt
	if (sgrp.eq."cub") then
          call Modulus_cubic(CCC(1,i), nEC, sgrp, B, E, G)
        else
          call Modulus(CCC(1,i), nEC, sgrp, B, E, G)
        endif
        write (6, '(f7.1, f11.6, 21f9.3)')
     &    t(i), vmin(i), bmin(i)*UtoGPa, 
     &    (CCC(j,i)*UtoGPa, j=1,nEC), 
     &    bcld(i)*UtoGPa,B*UtoGPa, E*UtoGPa, G*UtoGPa
      enddo
      end
c
      subroutine Modulus(CCC, nEC, sgrp, Bv, Ev, Gv)
c
      implicit real*8(a-h,o-z)
      real *8 CCC(nEC)
      character*(*) sgrp
      if (sgrp.eq."cub") then
        C11 = CCC(1)
        C22 = CCC(1)
        C33 = CCC(1)
        C12 = CCC(2)
        C13 = CCC(2)
        C23 = CCC(2)
        C44 = CCC(3)
        C55 = CCC(3)
        C66 = CCC(3)
      else if (sgrp.eq."hex") then
        C11 = CCC(1)
        C22 = CCC(1)
        C33 = CCC(4)
        C12 = CCC(2)
        C13 = CCC(3)
        C23 = CCC(3)
        C44 = CCC(5)
        C55 = CCC(5)
        C66 = (C11 - C12)/2.d0
c        call output (6, "# T(K) V(Ang^3) B(GPa) C11 C12 C13 C33 C44",
      else if (sgrp.eq."tet") then
        C11 = CCC(1)
        C22 = CCC(1)
        C33 = CCC(4)
        C12 = CCC(2)
        C13 = CCC(3)
        C23 = CCC(3)
        C44 = CCC(5)
        C55 = CCC(5)
        C66 = CCC(6)
c        call output (6, "# T(K) V(Ang^3) B(GPa) C11 C12 C13 C33 C44
c     &      C66",
      else if (sgrp.eq."ort") then
        C11 = CCC(1)
        C22 = CCC(4)
        C33 = CCC(6)
        C12 = CCC(2)
        C13 = CCC(3)
        C23 = CCC(5)
        C44 = CCC(6)
        C55 = CCC(7)
        C66 = CCC(8)
c        call output (6, "# T(K) V(Ang^3) B(GPa) C11 C12 C13 C22 C23
c     &      C33 C44 C55 C66",
      else if (sgrp.eq."mon") then
        C11 = CCC(1)
        C22 = CCC(5)
        C33 = CCC(8)
        C12 = CCC(2)
        C13 = CCC(3)
        C23 = CCC(6)
        C44 = CCC(10)
        C55 = CCC(12)
        C66 = CCC(13)
c        call output (6, "# T(K) V(Ang^3) B(GPa) C11 C12 C13 C16 C22
c     &      C23 C26 C33 C36 C44 C45 C55 C66",
      else if (sgrp.eq."tri") then
        C11 = CCC(1)
        C22 = CCC(7)
        C33 = CCC(12)
        C12 = CCC(2)
        C13 = CCC(3)
        C23 = CCC(8)
        C44 = CCC(16)
        C55 = CCC(19)
        C66 = CCC(21)
c        call output (6, "# T(K) V(Ang^3) B(GPa) C11 C12 C13 C14 C15
c     &      C16 C22 C23 C24 C25 C26 C33 C34 C35 C36 C44 C45 C46 C55
c     &      C56 C66",
      endif
      A = (C11 + C22 + C33)/3.d0
      B = (C12 + C13 + C23)/3.d0
      C = (C44 + C55 + C66)/3.d0
      Bv = (A + 2.d0*B)/3.d0
      Gv = (A - B + 3.d0*C)/5.d0
      Ev = 9.d0*Bv*Gv/(Gv+3.d0*Bv)
      end
c
      subroutine Modulus_Cubic(CCC, nEC, sgrp, Bv, Ev, Gv)
c
      implicit real*8(a-h,o-z)
      real *8 CCC(nEC)
      real*8 pcoeff(4)
      complex*16 zero(3)
        C11 = CCC(1)
        C12 = CCC(2)
        C44 = CCC(3)
        pcoeff(4) = 1.d0
        pcoeff(3) = (5.d0*C11 + 4.d0*C12)/8.d0
        pcoeff(2) = -C44*(7.d0*C11 - 4.d0*C12)/8.d0
        pcoeff(1) = -C44*(C11-C12)*(C11+2.d0*C12)/8.d0
        call DZPLRC(3, pcoeff, zero)
        root = real(zero(1))
        do i = 2, 3
          root = max(root, real(zero(i)))
        enddo
        Gv = root
        Bv = (C11+2.d0*C12)/3.d0
        Ev = 9.d0*Bv*Gv/(3.d0*Bv+Gv)
      end
c
      subroutine output(iunit, line, fmt)
c
      character*(*) line, fmt
      character*80, allocatable:: substring(:)
      jlen = len(line)
      jsub = 0
      i = 1
      do while (i.le.jlen)
        if (ichar(line(i:i)).gt.32) then
          do while (i.le.jlen)
            if (ichar(line(i:i)).le.32) then
              jsub = jsub + 1
              goto 100
            endif
            i = i + 1
          enddo
          jsub = jsub + 1
        endif
 100    i = i + 1
      enddo
      if (jsub.eq.0) return
      allocate (substring(jsub))
      jsub = 0
      i = 1
      do while (i.le.jlen)
        if (ichar(line(i:i)).gt.32) then
          ib = i
          do while (i.le.jlen)
            if (ichar(line(i:i)).le.32) then
              jsub = jsub + 1
              substring(jsub) = line(ib:i-1)
              goto 200
            endif
            i = i + 1
          enddo
          jsub = jsub + 1
          substring(jsub) = line(ib:i-1)
        endif 
 200    i = i + 1
      enddo
      write (6, fmt) (trim(substring(i)), i=1,jsub)
      return
      end
      

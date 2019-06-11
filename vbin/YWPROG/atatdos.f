      program electron_free_energy_of_wang_yi
c      USE NUMERICAL_LIBRARIES
c
c -------- energy in eV, temperature in K
c
      implicit real*8(a-h,o-z)
	common/comcon/bkb, tbeta, temper
	real*8, pointer:: pe(:), pdos(:), fn(:), afn(:)
	common/comstk/pe, pdos, fn, afn
	common/comvar/tnele, ndosmx
      real*8, allocatable, target :: e(:), dos(:)
	character*80 arg
      
	t0 = 0.d0
	t1 = 1300.d0
	td = 10.d0
	i = iargc()
	if (i.ge.1) then
	  call getarg(1, arg)
	  read (arg, *) t0
	endif
	if (i.ge.2) then
	  call getarg(2, arg)
	  read (arg, *) t1
	endif
	if (i.ge.3) then
	  call getarg(3, arg)
	  read (arg, *) td
	endif

	ndosmx=10001
	xdn = -10.d0
	xup = 10.d0
	allocate (e(ndosmx), dos(ndosmx), fn(ndosmx), afn(ndosmx)) 
	pe => e
	pdos => dos
	tnele = getdos(e, dos, ndosmx, xdn, xup)

	open (7, file="fvib_ele")

	bkb = dCONST('Boltzman')/dconst('ElectronVolt')
      t = 0.d0
	tbeta = 1.e30
      call caclf(f0, s0)
c	write (0, '(" NELE = ", f16.10, " f0 = ", f16.10)') tnele,f0

      do t = t0, t1, td
	  temper = t
	  if (t.le.1.d-5) then
	    tbeta = 1.e30
	  else
	    tbeta = 1.d0/(t*bkb)
	  endif
	  call caclf(f,s)
	  write (7,*) t,f-f0,s 
	enddo
	close (7)
      end
c
      real*8 function getdos(e, dos, nx, xdn, xup)
c
	implicit real*8 (a-h,o-z)
	real*8 e(nx), dos(nx)
      real*8, allocatable :: ve(:), vdos(:), ados(:)
	character*256 line
	do i = 1, 5
	  read (5,*)
	enddo
	read (5,*) eup, edn, nvdos, ef
	eup = eup - ef
	edn = edn - ef
	vde = (eup-edn)/dble(nvdos-1)
	allocate (ve(nvdos), vdos(nvdos), ados(nx), STAT=istat)
	if (istat.ne.0) then
	  print *, "Cannot allocate ve(nvdos), vdos(nvdos), ados(nx)"
	  stop
	endif
	v0 = edn
	do i = 1, nvdos
	  read (5, '(a)') line
	  if (len_trim(line).ge.59) then
	    read (line, *) t, t, t, vdos(i), x
		vdos(i) = vdos(i) + x
	  else
	    read (line, *) t, t, vdos(i)
	  endif
	  ve(i) = v0 
	  v0 = v0 + vde
	enddo
	xdn = max(xdn, edn)
	xup = min(xup, eup)
	xde = (xup - xdn) / dble(nx-1)
	ntemp = nint(xdn/xde)
	xdn = xde * dble(ntemp)
	v0 = xdn
	do i = 1, nx
	  e(i) = v0
	  v0 = v0 + xde
	enddo
	if (1.eq.1) then
	  CALL DCSIEZ (nvdos, ve, vdos, nx, e, ados)
        do i = 1, nx
	    dos(i) = DQDDER(1, e(i), nx, e, ados, .TRUE.)
	  enddo
	else
        do i = 1, nx
	    tx = e(i) - edn
	    kx = tx/vde + 1
	    kx = max (kx,1)
	    kx = min (nvdos-1, kx)
	    ados(i) = vdos(kx) + (vdos(kx+1) - vdos(kx))/vde
     &            *(e(i) - ve(kx))
	  enddo 
        do i = 1, nx
          i0 = min(i, nx-1)
          dos(i) = (ados(i0+1) - ados(i0))/xde
	  enddo
	endif

	call integr(dos,e,nx,ados,3)
	tx = 0.d0
	CALL DCSIEZ (nx, e, ados, 1, tx, t)
	getdos = t
	deallocate (ve, vdos, ados)
	return
	end
c
      subroutine caclf(f,s)
c
	implicit real*8 (a-h,o-z)
	common/comcon/bkb, tbeta, temper
	call findg(g)
	u = cacle(g)
      s = cacls(g)
	f = u - temper*s
      return
	end
c
	real*8 function gfind(val)
c
      implicit real*8 (a-h,o-z)
	common/comcon/bkb, tbeta, temper
	real*8, pointer:: pe(:), pdos(:), fn(:), afn(:)
	common/comstk/pe, pdos, fn, afn
	common/comvar/tnele, ndosmx
	
c      call integr(pdos,pe,6000,afn,3)
	k = ndosmx
	do i=1, ndosmx
	  tc = tbeta*(pe(i)-val)
	  if (tc.le.-200.d0) then
	    fn(i) = pdos(i)
	  else if (tc.ge.200.d0) then
	    fn(i) = 0.d0
	    k = i
	    goto 10
	  else
	    fn(i) = pdos(i)/(exp(tc)+1.d0)
	  endif
	enddo
 10   call integr(fn,pe,k,afn,3)
	gfind = afn(k) - tnele
	return
	end   
c
      subroutine findg(val)
c
	implicit real*8 (a-h,o-z)
	external gfind
	A = -0.5d0
	B =  0.5d0
	vA = gfind(A)
	vB = gfind(B)
	if (vA*vB.lt.0.d0) goto 10
	do i=1,16
	  A = 1.2*A
	  B = 1.2*B
 	  vA = gfind(A)
	  vB = gfind(B)
	  if (vA*vB.lt.0.d0) goto 10
	enddo
	write (6,*) ' I CANNOT FIND gibbs energy (mu) '
	stop
c	    print *,' vA =', vA, 
c     &           ',  vB=', vB
 10     ERRABS =0.0d0
        ERRREL = 1.d-12
	  MAXFN = 999
        CALL DZBREN (gfind, ERRABS, ERRREL, A, B, MAXFN)
	  val = B
	return
	end
c
	real*8 function cacle(val)
c
      implicit real*8 (a-h,o-z)
	common/comcon/bkb, tbeta, temper
	real*8, pointer:: pe(:), pdos(:), fn(:), afn(:)
	common/comstk/pe, pdos, fn, afn
	common/comvar/tnele, ndosmx
	
	k = ndosmx
	do i=1, ndosmx
	  tc = tbeta*(pe(i)-val)
	  if (tc.le.-200) then
	    fn(i) = pdos(i)*pe(i)
	  else if (tc.ge.200) then
	    fn(i) = 0.d0
	    k = i
	    goto 10
	  else
	    fn(i) = pe(i)*pdos(i)/(exp(tc)+1.d0)
	  endif
	enddo
 10   call integr(fn,pe,k,afn,3)
	cacle = afn(k)
	return
	end   
c
	real*8 function cacls(val)
c
      implicit real*8 (a-h,o-z)
	common/comcon/bkb, tbeta, temper
	real*8, pointer:: pe(:), pdos(:), fn(:), afn(:)
	common/comstk/pe, pdos, fn, afn
	common/comvar/tnele, ndosmx
	data zero,one/0.d0,1.d0/
	k = ndosmx
	do i=1, ndosmx
	  tc = tbeta*(pe(i)-val)
	  if (tc.le.-200.d0) then
	    fn(i) = 0.d0
	  else if (tc.ge.200.d0) then
	    fn(i) = 0.d0
	    k = i
	    goto 10
	  else
	    tf = 1.d0/(exp(tc)+1.d0)
	    if (tf.eq.zero.or.tf.eq.one) then
	      fn(i) = 0.d0
	    else
	      fn(i) = pdos(i)*(tf*log(tf)+(1.d0-tf)*log(1.d0-tf))
	    endif
	  endif
	enddo
 10   call integr(fn,pe,k,afn,3)
	cacls = -afn(k)*bkb
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

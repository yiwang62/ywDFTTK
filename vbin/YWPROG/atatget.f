      program atat_of_phonon
c      USE NUMERICAL_LIBRARIES
c
c -------- energy in Ryd, temperature in K, coordinate in a.u.
c
      implicit real*8(a-h,o-z)
      real*8, allocatable:: t(:), v(:,:), a(:)
      real*8, allocatable:: t_tem(:), v_tem(:)
      integer, allocatable:: nc(:)
      character*1024, allocatable:: file_name(:)
      character*1024 fn, buf*8096
c
      nvn = 0
      narg = iargc()
	ipstep = 1
	iskip = 0
	i = 1
      do while (i.le.narg)
        call getarg(i, fn)
        if (fn.eq."-ipstep") then
          i = i + 1
	    call getarg(i, fn)
	    read (fn, *) ipstep
        else if (fn.eq."-iskip") then
          i = i + 1
	    call getarg(i, fn)
	    read (fn, *) iskip
        else if (.not.isblank(fn)) then
		if (is_p_number(fn).lt.0) nvn = nvn + 1
        endif
	  i = i + 1
      enddo
	
	ipstep = max(1,ipstep)

      allocate (file_name(nvn), nc(nvn*2))
      do i = 1, nvn*2
        nc(i) = 0
      enddo
      nvn = 0
	i = 1
      do while (i.le.narg)
	  call getarg(i, fn)
        if (fn.eq."-ipstep") then
          i = i + 1
        else if (fn.eq."-iskip") then
          i = i + 1
        else if (.not.isblank(fn)) then
          kx = is_p_number(fn)
          if (kx.gt.0) then
            nvn2 = (nvn + 1)*2
            if ( nc(nvn2).eq.0 ) then
              nc(nvn2) = kx
              nc(nvn2-1) = 1
            else
              nc(nvn2-1) = nc(nvn2)
              nc(nvn2) = kx
            endif
          else if (kx.lt.0) then
            nvn = nvn + 1
            file_name(nvn) = fn
          endif
        endif
        i = i + 1
	enddo
      if (nc(1).eq.0) then
        nc(1) = 1
        nc(2) = 2
      endif
      do i = 2, nvn
        if (nc(i*2).eq.0) then
          nc(i*2-1) = nc(i*2-3)  
          nc(i*2) = nc(i*2-2)  
        endif
      enddo
      kp = 0
      do i = 1, 2*nvn
        kp = max(kp, nc(i))
      enddo
      allocate (a(kp))

      do i = 1,nvn
        nx = no_f_line(file_name(i))
        allocate (t_tem(nx), v_tem(nx))
        open (7, file=file_name(i))
        do j = 1, iskip
	    read (7, *)
	  enddo
        nt = 0
        do while (.TRUE.)
          read (7, '(a)', end=200) buf
          if (.not.isblank(buf).and.buf(1:1).ne."#") then
            nc1 = nc(i*2-1)
            nc2 = nc(i*2)
            m = max(nc1, nc2)
	      read (buf, *, end=200) (a(j), j=1,m)
            nt = nt + 1
            t_tem(nt) = a(nc1)
            v_tem(nt) = a(nc2)
          endif
        enddo
200     close (7)
          if (i.eq.1) then
	    ntn = nt
            allocate (t(ntn), v(ntn, nvn))
	    do j = 1, ntn
	      t(j) = t_tem(j)
	    enddo
	  endif
	if (nt.eq.1) then
          v(1,i) = v_tem(1)
	else
          call DCSIEZ (nt, t_tem, v_tem, ntn, t, v(1,i))
	endif
        deallocate (t_tem, v_tem)
      enddo
c
      jb = 1
      if (ntn.eq.1) jb = 2
      do i = 1,ntn
        if (mod(i,ipstep).eq.1.or.ipstep.eq.1) 
     &      write (6, '(99f18.8)') t(i), (v(i,j),j=jb,nvn)
      enddo
      end
c
      integer function no_f_line(fn)
c
      character*(*) fn, line*8
      open (7, file=fn, err=200, status="old")
      nt = 0
      do while (.TRUE.)
        read (7, '(a)', end=100) line
        if (.not.isblank(line).and.line(1:1).ne."#") nt = nt + 1
      enddo
 100  no_f_line = nt
      close(7)
	return
 200  print *, "********ERROR! CANNOT open file", trim(fn)
      stop
      end
c
      integer function is_p_number(str)
c
      character*(*) str
      read (str, *, err=100) ix
      goto 200
 100  ix = -1
 200  is_p_number = ix
      return
      end
c
      logical function isblank(str)
c
      character*(*) str
      do i=len(str), 1, -1
        if (ichar(str(i:i)).gt.32) goto 10
      enddo
 10   if (i.eq.0) then
        isblank = .TRUE.
        return
      endif
      if (str(1:1).eq."#") then
c        write (6, '(a)') str(1:i)
        isblank = .TRUE.
        return
      endif
      isblank = .FALSE.
      return
      end
      


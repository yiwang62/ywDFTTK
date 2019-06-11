      implicit real*8 (a-h,o-z)
      character*80 str
      real*8, allocatable, target:: x(:), y(:), v(:,:)
      integer, allocatable, target:: nat(:)
      character*8, allocatable, target:: sym(:)

      cpara = 1.0d0
      nofile = 0
      ntype = 0
      n = iargc()

      i = 0
      do while ( i<n )
        i = i + 1
        call getarg(i,str)
        if (str.eq."-cpara") then
          i = i + 1
          call getarg(i,str)
          read (str,*) cpara
        else if (str.eq."-n") then
          i = i + 1
          call getarg(i,str)
          read (str,*) nofile
        else if (str.eq."-ntype") then
          i = i + 1
          call getarg(i,str)
          read (str,*) ntype
          allocate (nat(ntype), sym(ntype), STAT=istat)
          if (istat.ne.0) then 
            print *, "2:pcoor Cannot allocate main dynamic memory"
            stop
          endif 
          do j = 1, ntype
            nat(j) = 0
            write (sym(j), '(a,i1)') "Element",j
          enddo
        else if (str.eq."-nat") then
          do j = 1, ntype
            i = i + 1
            call getarg(i,str)
            read (str,*) nat(j)
          enddo
        else if (str.eq."-atom") then
          do j = 1, ntype
            i = i + 1
            call getarg(i,sym(j))
          enddo
        endif
      enddo

      open(33, file='str_t.m.p.out')

c      write (0,*) "nofile =", nofile

      allocate (x(nofile), y(nofile), v(nofile, 3), STAT=istat)
      if (istat.ne.0) then 
        print *, "1:pcoor Cannot allocate main dynamic memory"
        stop
      endif 

      do i = 1, nofile
        read (5, *) x(i), str
        open (10+i, file=str)
      enddo

      call copy (nofile)
      call readx(y, nofile)
      call dcsiez(nofile,x,y,1,cpara,alat) 
      write (6,10) alat
 10   format (3f17.12, " ",a)
      do i = 1, 3
        call readv(v, nofile, *100)
        call dcsiez(nofile,x,v(1,1),1,cpara,b1) 
        call dcsiez(nofile,x,v(1,2),1,cpara,b2) 
        call dcsiez(nofile,x,v(1,3),1,cpara,b3) 
        write (6,10) b1,b2,b3
        write (33,10) b1*alat,b2*alat,b3*alat
      enddo
      call copy (nofile)
      call copy2 (nofile)
      write (33,'("1 0 0"/"0 1 0"/"0 0 1")')
      do i = 1, ntype
        do j = 1, nat(i)
          call readvv(v, nofile, *100)
          call dcsiez(nofile,x,v(1,1),1,cpara,b1) 
          call dcsiez(nofile,x,v(1,2),1,cpara,b2) 
          call dcsiez(nofile,x,v(1,3),1,cpara,b3) 
          write (6,10) b1,b2,b3,sym(i)(1:len_trim(sym(i)))
          write (33,10) b1,b2,b3,sym(i)(1:len_trim(sym(i)))
        enddo
      enddo
 100  continue
      end
c
      subroutine copy (n)
c
      implicit real*8 (a-h,o-z)
      character*80 str
      read (11, 10) str
 10   format (a)
      i = len_trim(str)
      write (6, 10) str(1:i)
      do i=12,10+n
        read (i,*)
      enddo
      end
c
      subroutine copy2(n)
c
      implicit real*8 (a-h,o-z)
      character*80 str
      read (11, 10) str
 10   format (a)
      i = len_trim(str)
      write (6, 10) str(1:i)
      do i=12,10+n
        read (i,*)
      enddo
      if (str(1:1).eq.'S'.or.str(1:1).eq.'s') call copy
      end
c
      subroutine readx(x, nofile)
c
      implicit real*8 (a-h,o-z)
      dimension x(nofile)
      do i=11,10+nofile
        read (i,*) x(i-10)
      enddo
      end
c
      subroutine readv(v,nofile,*)
c
      implicit real*8 (a-h,o-z)
      dimension v(nofile,3)
      do i=11,10+nofile
        read (i,*,end=100) (v(i-10,k),k=1,3)
      enddo
      return
 100  continue
      return 1
      end
c
      subroutine readvv(v,nofile,*)
c
      implicit real*8 (a-h,o-z)
      dimension v(nofile,3)
      do i=11,10+nofile
        read (i,*,end=100) (v(i-10,k),k=1,3)
      enddo
      do i = 1,3
        v0 = v(1,i)
        do j = 2, 5
          vd = v(j,i) - v0
          if (vd.ge.0.3d0) v(j,i) = v(j,i) - 1.0d0
          if (vd.le.-0.3d0) v(j,i) = v(j,i) + 1.0d0
          vd = v(j,i) - v0
          if (vd.ge.0.3d0.or.vd.le.-0.3d0) then
            print *, ' Something wrong in CONTCAR'
            stop
          endif
        enddo
      enddo
      return
 100  continue
      return 1
      end

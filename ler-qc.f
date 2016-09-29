      program qconst
      implicit none
      real*8 dr,qc,qcd,qcu,q,s,z,upd
      integer*4 k,i,j
      character(len=3) numb
      parameter(dr=0.08)

      do k = 25,50
      write(numb,'(I3.3)') k
      open(k,file=numb)
      enddo

      do j = 25,50
      qc=j*dr
      qcu=qc+0.0001d0
      qcd=qc-0.0001d0
      open(4,file='qsz7.txt',action='read',status='old')
      do i = 1,300107
      read(4,*) q,s,z,upd
      if (q .gt. qcd .and. q .lt. qcu) then
      write(j,*) q,s,z,upd
      endif
      enddo
      close(4)
      enddo

      end program
      

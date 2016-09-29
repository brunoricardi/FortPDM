      program read_diagonal
      implicit none
      real*8 dr,ref,r1,r2,cs,rho,refd,refu,up
      real*8  pi,lamb,tau,fc
      integer*4 i,j
      parameter(dr=0.04d0,pi=acos(-1.d0))
      parameter(lamb=12.1192393d0,tau=1.d0/40.d0)
      fc = (4.d0*pi*lamb*tau)**(1.5d0)

      open(3,file='ep7.txt')
      do j = 1,200
      ref=j*dr
      refd=ref-0.0001d0
      refu=ref+0.0001d0
      open(4,file='d7.txt',status='old',action='read')      
      do i = 1,19946
      read(4,*) r1,r2,cs,rho
      if ( (r1.gt.refd .and. r1.lt.refu) .and. 
     + (r2.gt.refd .and. r2.lt.refu) .and. 
     + (cs .gt. 0.99d0 .and. cs .lt. 1.01d0)) then
      up=-log(fc*rho)
      write(3,*) r1,up
      endif
      enddo
      close(4)
      enddo 
      close(3)


      end program




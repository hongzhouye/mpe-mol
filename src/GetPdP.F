      Subroutine GetPdP(N,No,Nstr,M,U,v,h0,wi,P,dP,E)
!     Build Density and dP/dv for a given v
!     M is the number of auxilliary basis functions (only tested for M=N so far)
      Use FCIUtils
      Use omp_lib
      Implicit None
      Integer*4 N,No,Nstr,M,Nv,N2,Max1,Max2,N0,i,j,k,ii,jj,kk,IfAll
      Real*8 U(N,N,N,N), wi(N,N,M)
      Real*8 v(M),h0(N,N),E,P(M)
!      Real*8 Tmp1(M),Tmp2(M),h(N,N)
      Real*8 dP(M,M),Fac, Eps,t1,t2,t1cpu,t2cpu
!      Real*8 X(Nstr,NStr),X0(NStr,NStr)
      Real*8,Dimension(:,:),ALLOCATABLE :: h,X,X0
      Real*8,Dimension(:),ALLOCATABLE :: Tmp1,Tmp2

      ALLOCATE(h(N,N))
      ALLOCATE(X(Nstr,Nstr))
      ALLOCATE(X0(Nstr,Nstr))
      ALLOCATE(Tmp1(M))
      ALLOCATE(Tmp2(M))

!      write ( *, '(a,i8)' ) '  The number of processors available = ',
!     \                       omp_get_num_procs()
!      write ( *, '(a,i8)' ) '  The number of threads available    = ',
!     \                       omp_get_max_threads()

      call omp_set_num_threads(4)
!      write(6,*)'threads serial',OMP_get_num_threads()

      P=0d0; dP=0d0
      h=h0
      Do i=1,M
         h=h+v(i)*wi(:,:,i)
      End Do

!     Solve for Ground State + Density
      N0=12
      Nv=N-No; N2=NcR(N,2); N0=Min(N0,NStr)
      Max1=NcR(N-1,No-1); Max2=NcR(N-2,No-2)
      X0=0d0

      Call FCIman(N,No,Nv,M,N2,Nstr,Max1,Max2,N0,h,h0,U,wi,X0,P,E)

!     Compute Density Derivatives by finite differences
      Eps=5d-4

      t1=omp_get_wtime()
      call cpu_time(t1cpu)

!$OMP PARALLEL DO FIRSTPRIVATE(h) PRIVATE(X,Tmp1,Tmp2,E)
      Do i=1,M
!          write(6,*)'threads parallel',OMP_get_num_threads()
          h=h+eps*wi(:,:,i); X=X0
          Call FCIman(N,No,Nv,M,N2,Nstr,Max1,Max2,N0,h,h0,U,wi,X,Tmp1,E)
          h=h-2d0*eps*wi(:,:,i); X=2d0*X0-X
          Call FCIman(N,No,Nv,M,N2,Nstr,Max1,Max2,N0,h,h0,U,wi,X,Tmp2,E)
          h=h+eps*wi(:,:,i)
          dP(:,i)=(Tmp1-Tmp2)/2d0/eps
      End Do
!$OMP END PARALLEL DO

      t2=omp_get_wtime()
      call cpu_time(t2cpu)

!      write(*,*) 'Fortran CPU time elapsed', t2cpu-t1cpu
!      write(*,*) 'OpenMP Walltime elapsed', t2-t1

!     Get Correct Density on exit
      Call FCIman(N,No,Nv,M,N2,Nstr,Max1,Max2,N0,h,h0,U,wi,X0,P,E)
      X=X0

      DEALLOCATE(h)
      DEALLOCATE(X)
      DEALLOCATE(X0)
      DEALLOCATE(Tmp1)
      DEALLOCATE(Tmp2)

      End subroutine GetPdP

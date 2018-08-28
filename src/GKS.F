      Subroutine GKS(N,No,M,h,h0,U,wi,C,P,E,diis)
      Implicit None
      Integer*4 N,No,M,i,j,k,l,kk,jter,diis
      Real*8 v(M),h0(N,N),wi(N,N,M),U(N,N,N,N),E
      Real*8 P(M),Tmp1(M),Tmp2(M),h(N,N)
      Real*8 Eps, PP(N,N),PPold(N,N),Resid
      Real*8 C(N,No),G(N,N),F(N,N),eigs(N),eigv(N,N)
      Real*8 errsave(N,N,diis),Fsave(N,N,diis),err(N,N)
      Real*8 cvec(diis+1),vecmut(diis+1,diis+1),bvec(diis+1)

!     Obtain initial density matrix
      PP=0d0
      Do i=1,No
          PP=PP+2d0*spread(C(:,i),dim=1,ncopies=N)
     \       *spread(C(:,i),dim=2,ncopies=N)
      end do

!     Set initial DIIS parameters
      Fsave=0d0; errsave=0d0; cvec=0d0; vecmut=0d0; bvec=0d0
      do i=1,diis
        vecmut(i,diis+1)=-1d0; vecmut(diis+1,i)=-1d0
      end do
      bvec(diis+1)=-1d0

!     Construct initial Fock matrix
      G=0d0
      do j=1,N; do k=1,N
        G=G+PP(j,k)*(U(:,:,k,j)-0.5d0*U(:,j,k,:))
      end do; end do
      F=h+G

      jter=0; Resid=1d0
      Do while (Resid.gt.1d-8.and.jter.lt.500)
          jter=jter+1

!     Diagonalize Fock matrix
          Call Diagonalize(N,F,eigs,eigv,kk)

!     Obtain new density matrix
          PPold=PP; PP=0d0
          Do i=1,No
            PP=PP+2d0*spread(eigv(:,i),dim=1,ncopies=N)
     \         *spread(eigv(:,i),dim=2,ncopies=N)
          end do

!     Construct new Fock matrix
          G=0d0
          do j=1,N; do k=1,N
            G=G+PP(j,k)*(U(:,:,k,j)-0.5d0*U(:,j,k,:))
          end do; end do
          F=h+G

!     Obtain HF energy
          E=0d0
          Do j=1,N; do k=1,N
            E=E+0.5d0*PP(k,j)*(h(j,k)+F(j,k))
          end do; end do

!     DIIS: save Fock matrix and err vector
!     Use the density matrix that constructs F for computing error vector
          err=MatMul(F,PP)-MatMul(PP,F)
          Do i=1,diis-1
            Fsave(:,:,i)=Fsave(:,:,i+1)
          end do
          Fsave(:,:,diis)=F
          Do i=1,diis-1
            errsave(:,:,i)=errsave(:,:,i+1)
          end do
          errsave(:,:,diis)=err

!     Determine DIIS coefficient
          vecmut(1:diis,1:diis)=0d0; cvec=0d0
          do i=1,diis; do j=1,diis
            do k=1,N; do l=1,N
              vecmut(i,j)=vecmut(i,j)+errsave(k,l,i)*errsave(k,l,j)
            end do; end do
          end do; end do
          if(jter.lt.diis)then
            Call GetLSQ(jter+1,jter+1,-bvec((diis+1-jter):),
     \                  vecmut((diis+1-jter):,(diis+1-jter):),
     \                  0d0,cvec((diis+1-jter):))
          else
            Call GetLSQ(diis+1,diis+1,-bvec,vecmut,0d0,cvec)
          end if

!     Construct new Fock matrix
          if(jter.gt.1)then
            F=0d0
            Do i=1,diis
              F=F+Fsave(:,:,i)*cvec(i)
            end do
          end if

!     Check density matrix convergence
          Resid=0d0
          do i=1,N; do j=1,N
            Resid=Resid+(PP(i,j)-PPold(i,j))**2
          end do; end do
          Resid=sqrt(Resid)/N

          if(jter.eq.500)then
            write(6,*)'#Warning Maximum Cycles Reached!!',jter,Resid,E
          end if
      End do
!      write(6,*)'#GKSjter,Resid,EHF',jter,Resid,E

!     Get correct orbitals, RI density and energy on exit
      C=eigv(:,1:No)
      P=0d0
      Do i=1,No
        do j=1,N; do k=1,N
          P=P+2d0*C(j,i)*C(k,i)*wi(j,k,:)
        end do; end do
      end do
      E=0d0
      Do j=1,N; do k=1,N
        E=E+0.5d0*PP(k,j)*(2d0*h0(j,k)+G(j,k))
      end do; end do

      End subroutine

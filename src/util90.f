      Subroutine GenDiagonalize(N,H,S,E,C)
!     Solve Symmetric Generalized Eigenvalue Problem H.C = E S.C
      Implicit None
      Integer*4 i,j,k,N,LWork, Info
      Real*8 H(N,N),S(N,N),E(N),C(N,N)
      Real*8 T(N,N),Hortho(N,N),Fac
      Real*8 VL(N,N),VR(N,N),Sig(N),Work(30*N+30)
     
!     Build T=S^-1/2
      LWork=30*N+30; T=S
      Call DGeSVD('A','A',N,N,T,N,Sig,VL,N,VR,N,Work,LWork,Info)
      T=0d0
      Do k=1,N
         If(Sig(k)>1d-8) Then
            Fac=1d0/Sqrt(Sig(k))
            Do i=1,N; Do j=1,N
               T(i,j)=T(i,j)
     $              +VL(j,k)*Fac*VR(k,i)
            End Do; End Do
         End If
      End Do
      
!     HOrtho=S^{-1/2}.H.S^{-1/2}
      HOrtho=MatMuL(T,MatMul(H,T))

!     Diagonalize HOrtho
      Call DSyEV('V','U',N,HOrtho,N,E,Work,LWork,i)

!     Convert Coefficients back to nonorthogonal basis
!      C=MatMul(T,HOrtho)
!     Leave Coefficients in Orthogonal basis
      C=HOrtho

      End Subroutine

      Subroutine SVDInverse(N,A,Ainv)
!     Compute the Pseudo-Inverse of A by SVD 
      Implicit None
      Integer*4 i,j,k,N,LWork,Info
      Real*8 A(N,N),Ainv(N,N)
      Real*8 VL(N,N),VR(N,N),S(N),Work(30*N+30)
      Real*8 Zero, One

      Zero=0d0; One=1d0; LWork=30*N+30

      Ainv=A
      Call DGeSVD('A','A',N,N,Ainv,N,S,VL,N,VR,N,Work,LWork,Info)
      Ainv=Zero
      Do i=1,N; Do j=1,N; Do k=1,N
         If(S(k)>1d-8) Ainv(i,j)=Ainv(i,j)
     $        +VL(j,k)/S(k)*VR(k,i)
      End Do; End Do; End Do

      End Subroutine

      Subroutine GetDet(N,A,Det)
C     Computes the Determinant of a general NxN matrix
      Implicit None
      Integer*4 N,LWork, Info, i
      Real*8 A(N,N),B(N,N),Det
      Real*8 Vl(N,N),Vr(N,N),Wi(N),Wr(N),Work(10*N+10)
      Complex*16 W(N),III, X
     
      LWork=10*N+10; Info=0; III=(0d0,1d0)
      B=A

      Call DGeEV('V','V',N,B,N,Wr,Wi,Vl,N,Vr,N,Work,LWork,Info)
 
      W=Wr+III*Wi; X=(1d0,0d0)
C     Compute Determinant
      Do i=1,N; X=X*W(i); End Do
      
      Det=X
      End Subroutine

      Integer*4 Function MyMod(N,M)
C     Computes N mod M with my preferred rounding.
      Integer*4 M,N,I

      I=N
      Do While (.not.(I.le.M.and.I.gt.0))
         If(N.gt.M)I=I-M
         If(N.le.0)I=I+M
      End Do

      MyMod=I

      End Function

      Subroutine ExpMat(N,D,U)
!     U=Exp(D)
      Implicit None
      Integer*4 i,j,k,N
      Real*8 U(N,N), D(N,N), D1(N,N), D2(N,N), D3(N,N)

!     Build exp(D/4096) to near-machine precision
      U=0d0; Do i=1,N; U(i,i)=1d0; End Do
      D1=D/4096d0; D2=MatMul(D1,D1); D3=MatMul(D2,D1)
      U=U+D1+D2/2d0+D3/6d0
      
!     Multiply repeatedly to get up to full U
      Do i=1,6
         D1=MatMul(U,U); 
         U=MatMul(D1,D1)
      End Do

      End Subroutine

      Subroutine LogQ(N,U,D)
!     Take the Log of a rotation matrix
!        U=exp(D)  with D real and skew-symmetric
!        Note: this doesn't deal with rotation+inversion
      Implicit None
      Integer*4 i,j,k,N, LWork, iUnit
      Real*8 U(N,N), D(N,N), Tmp(N,N), W(N,N), LogW(N,N)
      Real*8 WR(N),WI(N), VL(N,N), VR(N,N)
      Real*8 Work(10*N+10),eps

!     Diagonalize U
      Tmp=U; LWork=10*N+10; eps=1d-8
      Call DGeEV('V','V',N,Tmp,N,WR,WI,VL,N,VR,N,Work,LWork,i)

!     Check for eigenvalues of -1
      Do i=1,N; 
         If(Abs(WI(i)).lt.1d-8.and. WR(i).lt.0d0) 
     $        Write(6,*)'##Error! Inversion in U Matrix',i,WR(i)
      End Do

!     Fix Up Normalization of Eigenvectors
      Do i=1,N
         If(Abs(WI(i)).gt.eps) VL(:,i)=VL(:,i)*Sqrt(2d0)
         If(Abs(WI(i)).gt.eps) VR(:,i)=VR(:,i)*Sqrt(2d0)
      End Do

!     Take Log of W in canonical form (series of 2x2 real blocks)
!     Unless WR is 1, in which case the order gets messed up.
      LogW=0d0
      i=1
      Do While(i<N)
         If(Abs(WI(i)).gt. eps) Then
            LogW(i,i+1)=sign(acos(WR(i)),WI(i)) 
            LogW(i+1,i)=-LogW(i,i+1)
            i=i+2
         Else
            LogW(i,i)=0d0
            i=i+1
         End If
      End Do

!     D=V.LogW.V^{dagger}
      D=MatMul(VL,MatMul(LogW,Transpose(VL)))

!!     Check that This worked
!      Call ExpMat(N,D,Tmp)
!      Do i=1,N
!         Do j=1,N
!            If(Abs(U(i,j)-Tmp(i,j)).gt.1d-6)
!     $           Write(6,*)i,j,U(i,j),Tmp(i,j)
!         End Do
!      End Do

      End Subroutine
      Subroutine Diagonalize(N,A,Eigs,U,ierr)
      Implicit None
      Integer N,ierr,LWork
      Real*8 A(N,N),Eigs(N),U(N,N),Work(15*N+10)
! Get Eigenvalues and eigenvectors of a Real Symmetric Matrix

      LWork=15*N+10; U=A
      Call DSyEV('V','U',N,U,N,Eigs,Work,LWork,ierr)

      end subroutine

      Subroutine DgeMatMul(L,M,N,Y,Z,X,alpha,beta)
      Implicit None
!  X(L,N)= alpha*Y(L,M)*Z(M,N)+beta*(X(L,N)

      Integer*4 L,M,N
      Real*8 Y(L,M),Z(M,N),X(L,N),alpha,beta

!  ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )

      Call Dgemm('N','N',L,N,M,alpha,Y(1,1),L,Z(1,1),M,beta,X(1,1),L)

      Return
      End Subroutine


      Subroutine Axb(N,LDA,A,b,M,ierr)
      Implicit None
! Solves Ax=b in the most stable way possible
      Integer*4 N,LDA,M,Ipiv(N),IWork(N),ierr
      Real*8 A(LDA,N),b(N,M)
      Real*8 AF(N,N),R(N),C(N),X(N,M),Rcond,Ferr,Berr,Work(4*N)
      Character*1 EQ
      Real*8 zero,one,two,three,four
      Data zero,one,two,three,four /0.0d0,1.0d0,2.0d0,3.0d0,4.0d0/

! Clear out all the temporary arrays
      Ipiv=0; IWork=0
      AF=zero; R=zero; C=zero; X=zero; Rcond=zero; Ferr=zero 
      Berr=zero; Work=zero


!  Solve  Ax=b
!    In the first spot:
!     'N' - no equilibration
!     'E' - equilibrate A

      Call DGeSVX('N','N',N,M,A,LDA,AF,N,Ipiv,EQ,R,
     $               C,B,N,X,N,Rcond,Ferr,Berr,Work,IWork,ierr)

! Put the result back in b
      b=x
      
      Return
      End

      Real*8 Function factorial(n)
      Implicit None
! factorial - computes n!

      Integer*4 i,n

      factorial=1.0d0
      Do i=1,n
        factorial=factorial*i
      End Do

      Return
      End Function







      Subroutine GetLSQ(M,N,y,A,alpha,dx)
      Implicit None
!     Solves A.dx = y in a least squares sense (Using an SVD of A)
      Integer*4 i,j,k,N,M,LWork,info
      Real*8 y(M),A(M,N),dx(N),Ainv(M,N)
      Real*8 Tmp(M,N),U(M,M),VT(N,N),Sig(N),Work(10*M+50)
      Real*8 alpha,Beta
      Real*8 Zero,One,Two,Half

      Zero=0d0; One=1d0; Two=2d0; Half=0.5d0
      LWork=10*M+50
!      alpha=1d-6

      Tmp=A;Sig=Zero
      Call DGeSVD('A','A',M,N,A,M,Sig,U,M,VT,N,Work,LWork,info)
!      Write(6,107)Sig
! 107  Format(7f10.6)
!     Use Rational Function Aproximation to Inverse
      A=Zero; Ainv=Zero
      Do k=1,Min(N,M)
         Beta=Sig(k)/(Alpha+Sig(k)*Sig(k))
!         Beta=0d0
!         If(Sig(k).gt.Alpha)Beta=1d0/Sig(k)
         Do i=1,M; Do j=1,N
            A(i,j)=A(i,j)+U(i,k)*Sig(k)*VT(k,j)
            Ainv(i,j)=Ainv(i,j)+U(i,k)*VT(k,j)*Beta
         End Do; End Do
      End Do
!      Do i=1,N; Do j=1,M
!         Write(6,*)'#Check:',j,i,Tmp(j,i),A(j,i)
!      End Do; End Do

      dx=-MatMul(y,Ainv)
      
      End Subroutine

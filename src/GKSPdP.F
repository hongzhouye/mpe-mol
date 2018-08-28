      Subroutine GKSPdP(N,No,M,U,v,h0,wi,P,dP,E)
      Implicit None
      Integer*4 N,No,M,i,j,k
      Real*8 v(M),h0(N,N),wi(N,N,M),U(N,N,N,N),E
      Real*8 P(M),Tmp1(M),Tmp2(M),h(N,N)
      Real*8 dP(M,M), Eps
      Real*8 X(N,No),X0(N,No)

      P=0d0; dP=0d0
      h=h0
      Do i=1,M
         h=h+v(i)*wi(:,:,i)
      End Do

!     Solve for Ground State + Density
!     Default DIIS subspace is set to be 15
      X0=0d0
      Call GKS(N,No,M,h,h0,U,wi,X0,P,E,15)

!     Compute Density Derivatives by finite differences
      Eps=5d-4
      Do i=1,M
         h=h+eps*wi(:,:,i); X=X0
         Call GKS(N,No,M,h,h0,U,wi,X,Tmp1,E,15)
         h=h-2d0*eps*wi(:,:,i); X=2d0*X0-X
         Call GKS(N,No,M,h,h0,U,wi,X,Tmp2,E,15)
         h=h+eps*wi(:,:,i)

         dP(:,i)=(Tmp1-Tmp2)/2d0/eps
      End Do

!     Get Correct Density on exit
      Call GKS(N,No,M,h,h0,U,wi,X0,P,E,15)
      X=X0

      End subroutine GKSPdP

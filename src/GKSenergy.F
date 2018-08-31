      Subroutine GKSenergy(N,No,M,v,hh,VV,wi,Rho,Energy)
      Implicit None
      Integer*4 i,j,k,N,No,iter,M
      Real*8 v(M),Rho(M),wi(N,N,M)
      Real*8 Energy,C(N,N),Energynew
      Real*8 hh(N,N), VV(N,N,N,N)
      Real*8 Diff(M),dRho(M,M),dv(M),NewRho(M)
      Real*8 Norm,Resid,eps

!     Iterate to convergence
      Norm=1d0; Resid=1d0; iter=0
      Do While(Resid.gt.1d-8.and.iter.lt.100)
         iter=iter+1
         Call GKSPdP(N,No,M,VV,v,hh,wi,NewRho,dRho,Energy)

!     Get Step Direction
         Do i=1,M; Diff(i)=NewRho(i)-Rho(i); End Do
         Resid=0d0
         Do i=1,M; Resid=Resid+Diff(i)**2; End Do; Resid=Sqrt(Resid)
         if(Resid.gt.5d-3)then
           Call GetLSQ(M,M,Diff,dRho,1d-5,dv)
         else if(Resid.gt.1d-5.and.Resid.lt.5d-3) then
           Call GetLSQ(M,M,Diff,dRho,1d-12,dv)
         else
           Call GetLSQ(M,M,Diff,dRho,1d-18,dv)
         end if

!     Step along direction
         Norm=0d0; Do i=1,M; Norm=Norm+dv(i)**2; End Do; Norm=Sqrt(Norm)
         eps=1d0
         if(Resid.lt.1d-5)eps=2d0
         If(Norm.gt.eps)dv=dv/Norm*eps

!     First-order correction to energy when density is not perfectly matched
         Energynew=Energy
         Do i=1,M
           Energynew=Energynew-v(i)*(Rho(i)-NewRho(i))
         End Do

         v=v+dv

         Write(6,*)'##Convergence',iter,Resid,Norm,Energy,Energynew
         Call flush(6)
      End Do

!     First-order correction to energy when density is not perfectly matched
      v=v-dv
      Energynew=Energy
      Do i=1,M; Energynew=Energynew-v(i)*(Rho(i)-NewRho(i)); End Do

      Write(6,*)'##GKS Convergence',iter,Resid,Norm,Energy,Energynew
      Call flush(6)

      End subroutine GKSenergy

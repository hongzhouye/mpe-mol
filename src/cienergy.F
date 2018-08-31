      Subroutine cienergy(N,No,M,v,hh,VV,wi,Rho,Energy)
      Use FCIUtils
      Use omp_lib
      Implicit None
      Integer*4 i,j,k,l,jj,kk,N,No,ii,iter,NStr,IfAll,M
      Real*8 v(M),Rho(M),wi(N,N,M)
      Real*8 Energy,h(N,N),C(N,N),Energynew
      Real*8 vHF(M),hh(N,N), VV(N,N,N,N)
      Real*8 Diff(M),dRho(M,M),dv(M),NewRho(M)
      Real*8 Norm, Max, Resid,Trace, Value,Fac,eps,Pred,Tmp
      Real*8 tstart, tend, ostart, oend

!     Iterate to convergence
      Norm=1d0; Resid=1d0; iter=0
      Do While(Resid.gt.1d-7.and.iter.lt.50)
         iter=iter+1
         Nstr=NcR(N,No)

         Call CPU_time(tstart)
         ostart = omp_get_wtime()

         Call GetPdP(N,No,Nstr,M,VV,v,hh,wi,NewRho,dRho,Energy)

         Call CPU_time(tend)
         oend = omp_get_wtime()

         Write(6,'(A,F10.3,A)') "  CPU time for GetPdP:", 
     \     tend-tstart, " sec"
         Write(6,'(A,F10.3,A)') "  Wall time for GetPdP:", 
     \     oend-ostart, " sec"
         Write(6,*)
         Call flush(6)
      
!     Get Step Direction
         Call CPU_time(tstart)

         Do i=1,M; Diff(i)=NewRho(i)-Rho(i); End Do
         Resid=0d0
         Do i=1,M; Resid=Resid+Diff(i)**2; End Do; Resid=Sqrt(Resid)
         if(Resid.gt.1d-3)then
           Call GetLSQ(M,M,Diff,dRho,1d-6,dv)
         else if(Resid.gt.5d-4.and.Resid.lt.1d-3) then
           Call GetLSQ(M,M,Diff,dRho,1d-12,dv)
         else
           Call GetLSQ(M,M,Diff,dRho,1d-15,dv)
         end if

         Call CPU_time(tend)

         Write(6,'(A,F10.3,A)') "  CPU time for LnSrch:", 
     \     tend-tstart, " sec"
         Write(6,*)
         Call flush(6)

!     Step along direction
         Norm=0d0; Do i=1,M; Norm=Norm+dv(i)**2; End Do; Norm=Sqrt(Norm)
         eps=1d0
         if(Resid.lt.1d-5)eps=2d0
         If(Norm.gt.eps)dv=dv/Norm*eps

!     First-order correction to energy when density is not perfectly matched
         Energynew=Energy
         Do i=1,M
           Energynew=Energynew-2d0*v(i)*(Rho(i)-NewRho(i))
         End Do

         v=v+dv

         Write(6,*)'##Convergence',iter,Resid,Norm,Energy,Energynew
         Call flush(6)
      End Do

      v=v-dv
      Write(6,*)'##CI Convergence',iter, Resid,Norm,Energy,Energynew
      Call flush(6)

      End subroutine cienergy

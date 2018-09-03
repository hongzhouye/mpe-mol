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
      Real*8 threshs(3), lbds(3), epss(2)
      logical file_present

! @@HY: read regularization and step length from file "ls_params.txt"
!       if file not exist, use default values

      Inquire(File="ls_params.txt", Exist=file_present )
      If(file_present) then
        Write(6, '(A)') "  LinSrch parameter file detected!"
        Write(6, '(A)') "  Reading parameters from file:"
        Open(Unit=21, File="ls_params.txt")
        Read(21, *) threshs
        Read(21, *) lbds
        Read(21, *) epss
        Close(21)
      Else
        Write(6, '(A)') "  LinSrch parameter file not found!"
        Write(6, '(A)') "  Using default parameters:"
        threshs = (/ 1d-3, 5d-4, 1d-5 /)
        lbds = (/ 1d-6, 1d-12, 1d-15 /)
        epss = (/ 2d0, 1d0 /)
      End If

      Write(6, '(A, E8.1, E8.1, E8.1)') "  threshs : ", threshs
      Write(6, '(A, E8.1, E8.1, E8.1)') "  lbds    : ", lbds
      Write(6, '(A, E8.1, E8.1)')       "  epss    : ", epss
      Call flush(6)

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
         Call flush(6)
      
!     Get Step Direction
         Do i=1,M; Diff(i)=NewRho(i)-Rho(i); End Do
         Resid=0d0
         Do i=1,M; Resid=Resid+Diff(i)**2; End Do; Resid=Sqrt(Resid)
         if(Resid.gt.threshs(1))then
           Call GetLSQ(M,M,Diff,dRho,lbds(1),dv)
         else if(Resid.gt.threshs(2).and.Resid.lt.threshs(1)) then
           Call GetLSQ(M,M,Diff,dRho,lbds(2),dv)
         else
           Call GetLSQ(M,M,Diff,dRho,lbds(3),dv)
         end if

!     Step along direction
         Norm=0d0; Do i=1,M; Norm=Norm+dv(i)**2; End Do; Norm=Sqrt(Norm)
         eps=epss(2)
         if(Resid.lt.threshs(3))eps=epss(1)
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

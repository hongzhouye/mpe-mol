      Subroutine GKSenergy(N,No,M,v,hh,VV,wi,Rho,Energy,vgks_name)
      Implicit None
      Integer*4 i,j,k,N,No,iter,M
      Real*8 v(M),Rho(M),wi(N,N,M)
      Real*8 Energy,C(N,N),Energynew
      Real*8 hh(N,N), VV(N,N,N,N)
      Real*8 Diff(M),dRho(M,M),dv(M),NewRho(M)
      Real*8 Norm,Resid,eps
      Real*8 threshs(3), lbds(3), epss(2)
      Character*200 vgks_name
      logical file_present

! @@HY: read regularization and step length from file "ls_params_gks.txt"
!       if file not exist, use default values

      Inquire(File="ls_params_gks.txt", Exist=file_present )
      If(file_present) then
        Write(6, '(A)') "  LinSrch parameter file detected!"
        Write(6, '(A)') "  Reading parameters from file:"
        Open(Unit=21, File="ls_params_gks.txt")
        Read(21, *) threshs
        Read(21, *) lbds
        Read(21, *) epss
        Close(21)
      Else
        Write(6, '(A)') "  LinSrch parameter file not found!"
        Write(6, '(A)') "  Using default parameters:"
        Write(6, '(A)') 
     \       "  [Warning] The default is generally good as"
     \    // " the GKS inversion is less singular compared to the FCI"
     \    // " inversion, which justifies the use of smaller lbd such"
     \    // " as 1d-12. However, for hard cases such as Cl- (large"
     \    // " anion), even the GKS inversion becomes singular and we"
     \    // " need larger lbd's. Recommended values:"
        Write(6, '(A)') "    threshs = (/ 1d-4, 1d-6, 1d-4 /)"
        Write(6, '(A)') "    lbds    = (/ 1d-5, 1d-7, 1d-9 /)"
        Write(6, '(A)') "    epss    = (/ 2d0, 1e0 /) [End Warning]"

        threshs = (/ 5d-3, 1d-5, 1d-5 /)
        lbds = (/ 1d-5, 1d-12, 1d-18 /)
        epss = (/ 2d0, 1d0 /)
      End If

      Write(6, '(A, E8.1, E8.1, E8.1)') "  threshs : ", threshs
      Write(6, '(A, E8.1, E8.1, E8.1)') "  lbds    : ", lbds
      Write(6, '(A, E8.1, E8.1)')       "  epss    : ", epss
      Write(6, *)
      Call flush(6)

!     Iterate to convergence
      Norm=1d0; Resid=1d0; iter=0
!      Do While(Resid.gt.1d-8.and.iter.lt.100)
      Do While(Resid.gt.1d-8.and.iter.lt.50)
         iter=iter+1
         Call GKSPdP(N,No,M,VV,v,hh,wi,NewRho,dRho,Energy)

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
           Energynew=Energynew-v(i)*(Rho(i)-NewRho(i))
         End Do

         v=v+dv

! @@HY: write vgks to file at each iteration!
!       make it easy for restarting the calculation!
         Open(Unit=37, File=trim(vgks_name))
         Write(37, *) v
         Close(37)

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

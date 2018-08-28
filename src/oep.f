!     @@HY: Added new variable 'FuncType' and 'Kxx'      
!     FuncType = 0 --> EXX (HF)
!            1 --> LDA
!            2 --> PBE
!            3 --> PBE0 (Kxx will be set to 0.25)
!            4 --> BLYP 
!            13 --> PBE0 w/ user specified amount of EXX
!
!     Kxx    portion of EXX, must be specified for FuncType=13
!            Kxx will be set to 0.25 for FuncType=3,
!                               0.00 for 0<FuncType<3,
!                               1.00 for FuncType=0. 

      SUBROUTINE oep(data_path,coeff,vx,nocc,modims,auxdims,auxmomo,
     \ oneeint,oneekin,oneenuc,VV,DFTmo,FuncType,Kxx,Energy)

! ------------------------------------------------
! @@HY 08/26/18
!                 OEP CHEAT-SHEET
! input:
!     coeff      [nbas, nocc] density expanded in MO basis  
!     auxmomo    [nbas, nbas, naux] (ij|P)
!     VV         [nbas**4] (ij|kl)
!     DFTmo      [nbas, nbas] MO coeff
!
! output:
!     vx         [naux] oep potential
!     energy     oep energy
! ------------------------------------------------

      Character*200 data_path
      Integer*4 nocc,modims,auxdims,Natom
      Real*8 coeff(modims,nocc), vx(auxdims),Ek,Energy
      Real*8 auxmomo(modims,modims,auxdims), oneeint(modims,modims)
      Real*8 oneekin(modims,modims), oneenuc(modims,modims) 
      Integer*4 i,j
      Real*8 nuc_cons,alpha,beta,tmp1d(modims),W
      Real*8 VV(modims,modims,modims,modims),DFTmo(modims,modims)
      Integer*4 FuncType
      Real*8 Kxx

      alpha=1d0; beta=0d0
!     Compute the external potential energy (constant)
      nuc_cons=0d0
      do i=1,nocc
          tmp1d=0d0
          Call DGemm("N","N",1,modims,modims,alpha,
     \               coeff(:,i),1,oneenuc,modims,beta,tmp1d,1)
          Call DGemm("N","N",1,1,modims,-2d0,
     \               coeff(:,i),1,tmp1d,modims,1d0,nuc_cons,1)
      end do

!     Newton's Optimization to search for OEP
      Call newton(vx,coeff,nocc,modims,auxdims,auxmomo,
     \ oneeint,oneekin,oneenuc,1d-12,1000)

!     Final Lagrangian and kinetic energy
      Call func(vx,coeff,nocc,modims,auxdims,auxmomo,
     \ oneeint,oneekin,oneenuc,nuc_cons,W)

      Call kinetic(vx,coeff,nocc,modims,auxdims,auxmomo,
     \ oneeint,oneekin,oneenuc,nuc_cons,Ek)

!     W-T is 0 for exact density matching
      Write(6,*)'#final W T',W,Ek

!     Compute DFT energy using grids
      open(unit=21,File=trim(data_path) // '/grid_dimension.txt')
      read(21,*)Natom
      close(21)

      Call DFTEnergy(data_path,vx,coeff,nocc,modims,auxdims,auxmomo,
     \ oneeint,VV,DFTmo,Natom,FuncType,Kxx,Energy)

      End subroutine oep

!     Compute the DFT energy

!     @@HY: added new input variable 'FuncType' and 'Kxx'
!     FuncType = 0 --> EXX (HF)
!            1 --> LDA
!            2 --> PBE
!            3 --> PBE0 (Kxx will be set to 0.25)
!            13 --> PBE0 w/ user specified amount of EXX
!
!     Kxx    portion of EXX
!            Kxx will be set to 0.25 for FuncType=3, and 0.0 for FuncType<3. 
      Subroutine DFTEnergy(data_path,vx,coeff,nocc,modims,auxdims,
     \  auxmomo,oneeint,VV,DFTmo,Natom,FuncType,Kxx,Energy)

      Character*200 data_path
      Integer*4 nocc,modims,auxdims,i,j,k,l,a,b,P,ierr
      Integer*4 Npoints(Natom),Nptot,Nbasis,Natom
      Integer*4 FuncType,Func_ids(2)
      Real*8 coeff(modims,nocc),vx(auxdims),alpha,beta
      Real*8 auxmomo(modims,modims,auxdims), oneeint(modims,modims)
      Real*8 F(modims,modims),tmp2d(modims,modims),tmp1d(modims)
      Real*8 eigs(modims),eigv(modims,modims),rdm1(modims,modims)
      Real*8 VV(modims,modims,modims,modims),DFTmo(modims,modims)
      Real*8,allocatable :: Grid(:,:), aogrid(:,:)
      Real*8 Kxx,Energy,E1e,EJ,Exx,Exc(3),tmp0d

      Real*8, parameter :: PI = 4.d0*atan(1.d0)

!     Read numbers of grid points, atomic basis functions and atoms
      open(unit=21,File=trim(data_path) // '/grid_dimension.txt')
      read(21,*)Natom,Npoints,Nbasis
      close(21)

      Nptot=0
      do i=1,Natom
        Nptot=Nptot+Npoints(i)
      end do

!     Read grid coordinates and weights
      open(unit=21,File=trim(data_path) // '/grid.txt')
      allocate(Grid(Nptot,4))
      do i=1,4
        read(21,*) Grid(:,i)
      end do
      close(21)

!     Read AO basis values on grids
      open(unit=21,File=trim(data_path) // '/ao_grid.txt')
      allocate(aogrid(Nptot*Nbasis,4))
      do i=1,4
        read(21,*) aogrid(:,i)
      end do
      close(21)

!     Construct Fock matrix using OEP potential
!     @@HY: for the total system OEP, this is not the same
!     as computed from the exact Fock matrix due to the finite
!     accuracy of OEP
      tmp2d=0d0
      do P=1,auxdims
          tmp2d=tmp2d+vx(P)*auxmomo(:,:,P)
      end do
      F=oneeint+tmp2d
      Call Diagonalize(modims,F,eigs,eigv,ierr)

!     Compute 1e energy using OEP orbitals
      E1e = 0d0
      Do i=1,modims; Do j=1,modims
        Do k=1,nocc
          E1e = E1e + eigv(i,k)*oneeint(i,j)*eigv(j,k)
        End Do
      End Do; End Do
      E1e = E1e * 2d0

!     Compute Coulomb repulsion using OEP orbitals
      EJ = 0d0
      Do a=1,nocc; Do b=1,nocc
        Do i=1,modims; Do j=1,modims
        Do k=1,modims; Do l=1,modims
          EJ = EJ + eigv(i,a)*eigv(j,a)*eigv(k,b)*eigv(l,b)*VV(i,j,k,l)
        End Do; End Do; 
        End Do; End Do;
      End Do; End Do
      EJ = EJ * 2d0

!     Compute EXX using OEP orbitals
      EXX = 0d0
      Do a=1,nocc; Do b=1,nocc
        Do i=1,modims; Do j=1,modims
        Do k=1,modims; Do l=1,modims
          EXX = EXX + eigv(i,a)*eigv(j,a)*eigv(k,b)*eigv(l,b)*
     \      VV(i,l,k,j)
        End Do; End Do; 
        End Do; End Do;
      End Do; End Do
      EXX = EXX * (-1d0)

!     Construct rdm1 in AO basis 
      Do i=1,modims; Do j=1,i
        rdm1(i,j) = 0d0
        Do k=1,modims; Do l=1,modims
          tmp0d = 0d0
          Do a=1,nocc
            tmp0d = tmp0d + eigv(k,a)*eigv(l,a)
          End Do
          rdm1(i,j) = rdm1(i,j) + DFTmo(i,k)*DFTmo(j,l)*tmp0d
        End Do; End Do
        rdm1(i,j) = rdm1(i,j) * 2d0
        rdm1(j,i) = rdm1(i,j)
      End Do; End Do

!     Call C function 'eval_xc'
      Exc = 0d0
      If(FuncType .eq. 0) then
        Write(*,*) "Requesting EXX; skipping C function 'eval_xc'" 
        Kxx = 1d0
        Go to 311
      Else If(FuncType .eq. 1) then
        Write(*,*) "Requesting LDA; calling C function 'eval_xc'"
        Func_ids(1) = 1
        Func_ids(2) = 7
        Kxx = 0d0
      Else If(FuncType .eq. 2) then
        Write(*,*) "Requesting PBE; calling C function 'eval_xc'"
        Func_ids(1) = 101
        Func_ids(2) = 130
        Kxx = 0d0
      Else If(FuncType .eq. 3) then
        Write(*,*) "Requesting PBE0; calling C function 'eval_xc'"
        Func_ids(1) = 101
        Func_ids(2) = 130
        Kxx = 2.5d-1
      Else If(FuncType .eq. 4) then
        Write(*,*) "Requesting BLYP; calling C function 'eval_xc'"
        Func_ids(1) = 106
        Func_ids(2) = 131
        Kxx = 0d0
      Else If(FuncType .eq. 13) then
        Write(*,*) "Requesting PBE0*; calling C function 'eval_xc'"
        Func_ids(1) = 101
        Func_ids(2) = 130
      Else
        Write(*,'(A,I4)') "[ERROR] Requesting unknown func ", FuncType
        Call EXIT(1)
      End If
      Ierr = eval_xc(Nptot, modims, nocc, 2, rdm1, Grid, 
     \  aogrid, Func_ids, Exc)

311   Continue

!     Collect various contributions to energy
      Energy = (1d0-Kxx)*Exc(1)+Exc(2)+E1e+EJ+Kxx*EXX

      Write(*,*) "  Energy components"
      Write(*,*) "---------------------"
      Write(*,'(A,F12.6)') "  Kxx        = ", Kxx
      Write(*,'(A,F12.6)') "  E1e        = ", E1e
      Write(*,'(A,F12.6)') "  EJ         = ", EJ
      Write(*,'(A,F12.6)') "  raw EXX    = ", EXX
      Write(*,'(A,F12.6)') "  scaled EXX = ", EXX*Kxx
      write(*,'(A,F12.6)') "  raw Ex     = ", Exc(1)
      write(*,'(A,F12.6)') "  scaled Ex  = ", Exc(1)*(1d0-Kxx)
      write(*,'(A,F12.6)') "  Ec         = ", Exc(2)
      Write(*,'(A,F12.6)') "  EDFT total = ", Energy
      Write(*,*) "---------------------"

      deallocate(Grid)
      deallocate(aogrid)

      End Subroutine DFTEnergy

!     Compute the Lagragian
      Subroutine func(vx,coeff,nocc,modims,auxdims,auxmomo,
     \ oneeint,oneekin,oneenuc,nuc_cons,W)
      
      Integer*4 nocc,modims,auxdims,i,j,k,P
      Real*8 coeff(modims,nocc), vx(auxdims),nuc_cons,alpha,beta
      Real*8 auxmomo(modims,modims,auxdims), oneeint(modims,modims)
      Real*8 oneekin(modims,modims), oneenuc(modims,modims),W
      Real*8 F(modims,modims),tmp2d(modims,modims),tmp1d(modims)
      Real*8 tmp1d2(auxdims),eigs(modims),eigv(modims,modims)

      alpha=1d0; beta=0d0
!     Construct Fock matrix and diagonalize
      tmp2d=0d0
      do P=1,auxdims
          tmp2d=tmp2d+vx(P)*auxmomo(:,:,P)
      end do
      F=oneeint+tmp2d
      Call Diagonalize(modims,F,eigs,eigv,k)

      W=nuc_cons
      do i=1,nocc
          tmp1d=0d0
          Call DGemm("N","N",1,modims,modims,alpha,
     \               eigv(:,i),1,oneeint,modims,beta,tmp1d,1)
          Call DGemm("N","N",1,1,modims,2d0,
     \               eigv(:,i),1,tmp1d,modims,1d0,W,1)
      end do

      do i=1,nocc
        tmp1d2=0d0
        do P=1,auxdims
          tmp1d=0d0
          Call DGemm("N","N",1,modims,modims,alpha,
     \               eigv(:,i),1,auxmomo(:,:,P),modims,beta,tmp1d,1)
          Call DGemm("N","N",1,1,modims,alpha,
     \               eigv(:,i),1,tmp1d,modims,beta,tmp1d2(P),1)
        end do
        Call DGemm("N","N",1,1,auxdims,2d0,vx,1,tmp1d2,auxdims,1d0,W,1)
      end do

      do i=1,nocc
        tmp1d2=0d0
        do P=1,auxdims
          tmp1d=0d0
          Call DGemm("N","N",1,modims,modims,alpha,
     \               coeff(:,i),1,auxmomo(:,:,P),modims,beta,tmp1d,1)
          Call DGemm("N","N",1,1,modims,alpha,
     \               coeff(:,i),1,tmp1d,modims,beta,tmp1d2(P),1)
        end do
        Call DGemm("N","N",1,1,auxdims,-2d0,vx,1,tmp1d2,auxdims,1d0,W,1)
      end do

      W=-W

      End subroutine func


!     Compute the gradient of Lagrangian (density matching)
      Subroutine gradient(vx,coeff,nocc,modims,auxdims,auxmomo,
     \ eigs,eigv,deriv)

      Integer*4 nocc,modims,auxdims,i,j,k,P
      Real*8 coeff(modims,nocc), vx(auxdims),alpha,beta
      Real*8 auxmomo(modims,modims,auxdims), oneeint(modims,modims)
      Real*8 oneekin(modims,modims), oneenuc(modims,modims),W
      Real*8 F(modims,modims),tmp2d(modims,modims),tmp1d(modims)
      Real*8 tmp1d2(auxdims),eigs(modims),eigv(modims,modims)
      Real*8 deriv(auxdims)

      alpha=1d0; beta=0d0
      deriv=0d0
      do i=1,nocc
        do P=1,auxdims
          tmp1d=0d0
          Call DGemm("N","N",1,modims,modims,alpha,
     \               eigv(:,i),1,auxmomo(:,:,P),modims,beta,tmp1d,1)
          Call DGemm("N","N",1,1,modims,2d0,
     \               eigv(:,i),1,tmp1d,modims,1d0,deriv(P),1)
        end do
      end do

      do i=1,nocc
        do P=1,auxdims
          tmp1d=0d0
          Call DGemm("N","N",1,modims,modims,alpha,
     \               coeff(:,i),1,auxmomo(:,:,P),modims,beta,tmp1d,1)
          Call DGemm("N","N",1,1,modims,-2d0,
     \               coeff(:,i),1,tmp1d,modims,1d0,deriv(P),1)
        end do
      end do

      deriv=-deriv

      End subroutine gradient


!     Compute the hessian of Lagrangian (needed for Newton's method)
      Subroutine hessian(vx,coeff,nocc,modims,auxdims,auxmomo,
     \ eigs,eigv,hess)

      Integer*4 nocc,modims,auxdims,i,j,k,P,Q,a
      Real*8 coeff(modims,nocc), vx(auxdims),alpha,beta
      Real*8 auxmomo(modims,modims,auxdims), oneeint(modims,modims)
      Real*8 oneekin(modims,modims), oneenuc(modims,modims),W
      Real*8 F(modims,modims),tmp2d(modims,modims),tmp1d(modims)
      Real*8 tmp1d2(auxdims),eigs(modims),eigv(modims,modims)
      Real*8 hess(auxdims,auxdims),der2(modims,auxdims)
      Real*8 der1(modims,modims,auxdims)

      alpha=1d0; beta=0d0
      hess=0d0
!      do i=1,nocc
!        do P=1,auxdims
!          tmp2d=0d0
!          Call DGemm("N","N",modims,modims,modims,alpha,
!     \             auxmomo(:,:,P),modims,eigv,modims,beta,tmp2d,modims)
!          Call DGemm("N","N",1,modims,modims,alpha,
!     \             eigv(:,i),1,tmp2d,modims,beta,der2(:,P),1)
!        end do
!        do a=nocc+1,modims
!          do Q=1,auxdims; do P=1,auxdims
!            hess(P,Q)=hess(P,Q)+4d0*der2(a,P)*der2(a,Q)
!     $                    /(eigs(i)-eigs(a))
!          end do; end do
!        end do
!      end do

      do P=1,auxdims
        tmp2d=0d0
        Call DGemm("N","N",modims,modims,modims,alpha,
     \             auxmomo(:,:,P),modims,eigv,modims,beta,tmp2d,modims)
        Call DGemm("N","N",modims,modims,modims,alpha,transpose(eigv),
     \             modims,tmp2d,modims,beta,der1(:,:,P),modims)
      end do
      do i=1,nocc
        do a=nocc+1,modims
          do Q=1,auxdims; do P=1,auxdims
            hess(P,Q)=hess(P,Q)+4d0*der1(i,a,P)*der1(i,a,Q)
     $                    /(eigs(i)-eigs(a))
          end do; end do
        end do
      end do

      hess=-hess

      End subroutine hessian


!     Newton's method
      Subroutine newton(vx,coeff,nocc,modims,auxdims,auxmomo,
     \ oneeint,oneekin,oneenuc,gtol,maxiter)

      Integer*4 nocc,modims,auxdims,i,j,k,P,maxiter
      Real*8 coeff(modims,nocc), vx(auxdims),alpha,beta
      Real*8 auxmomo(modims,modims,auxdims), oneeint(modims,modims)
      Real*8 oneekin(modims,modims), oneenuc(modims,modims),W
      Real*8 F(modims,modims),tmp2d(modims,modims),tmp1d(modims)
      Real*8 tmp1d2(auxdims),eigs(modims),eigv(modims,modims)
      Real*8 hess(auxdims,auxdims),deriv(auxdims),gtol
      Real*8 norm_grad,dvx(auxdims),step,eps

!     Construct Fock matrix and diagonalize
      tmp2d=0d0
      do P=1,auxdims
          tmp2d=tmp2d+vx(P)*auxmomo(:,:,P)
      end do
      F=oneeint+tmp2d
      Call Diagonalize(modims,F,eigs,eigv,k)

      Call gradient(vx,coeff,nocc,modims,auxdims,auxmomo,
     \              eigs,eigv,deriv)
      norm_grad=0d0
      do i=1,auxdims; norm_grad=norm_grad+deriv(i)**2; end do
      norm_grad=sqrt(norm_grad)

      jter=0
      do while (norm_grad.gt.gtol .and. jter.lt.maxiter)
          jter=jter+1

          Call hessian(vx,coeff,nocc,modims,auxdims,auxmomo,
     \                 eigs,eigv,hess)

!         Get the search direction
          if (norm_grad.gt.1d-3) then
              Call GetLSQ(auxdims,auxdims,deriv,hess,1d-6,dvx)
          else if (norm_grad.gt.1d-5.and.norm_grad.lt.1d-3) then
              Call GetLSQ(auxdims,auxdims,deriv,hess,1d-12,dvx)
          else
              Call GetLSQ(auxdims,auxdims,deriv,hess,1d-18,dvx)
          end if

!         Trust radius
          step=0d0
          do i=1,auxdims; step=step+dvx(i)**2; end do; step=sqrt(step)
          eps=0.5d0
          if(norm_grad.lt.1d-5)eps=2d0
          if (step.gt.eps) then
              dvx=dvx/step*eps
          end if

!         Update potential
          vx=vx+dvx

!         Update Fock matrix
          tmp2d=0d0
          do P=1,auxdims
              tmp2d=tmp2d+vx(P)*auxmomo(:,:,P)
          end do
          F=oneeint+tmp2d
          Call Diagonalize(modims,F,eigs,eigv,k)

          Call gradient(vx,coeff,nocc,modims,auxdims,auxmomo,
     \                  eigs,eigv,deriv)
          norm_grad=0d0
          do i=1,auxdims; norm_grad=norm_grad+deriv(i)**2; end do
          norm_grad=sqrt(norm_grad)

          if ((jter/100)*100.eq.jter) then
              write(6,*)'#jter,norm_grad,step',jter,norm_grad,step
          end if
      end do

      write(6,*)'#jter,norm_grad,step',jter,norm_grad,step

      End subroutine newton


!     Compute kinetic energy
      Subroutine kinetic(vx,coeff,nocc,modims,auxdims,auxmomo,
     \ oneeint,oneekin,oneenuc,nuc_cons,Ek)

      Integer*4 nocc,modims,auxdims,i,j,k,P
      Real*8 coeff(modims,nocc), vx(auxdims),nuc_cons,alpha,beta
      Real*8 auxmomo(modims,modims,auxdims), oneeint(modims,modims)
      Real*8 oneekin(modims,modims), oneenuc(modims,modims),Ek
      Real*8 F(modims,modims),tmp2d(modims,modims),tmp1d(modims)
      Real*8 tmp1d2(auxdims),eigs(modims),eigv(modims,modims)

      alpha=1d0; beta=0d0
!     Construct Fock matrix and diagonalize
      tmp2d=0d0
      do P=1,auxdims
          tmp2d=tmp2d+vx(P)*auxmomo(:,:,P)
      end do
      F=oneeint+tmp2d
      Call Diagonalize(modims,F,eigs,eigv,k)

      Ek=0d0
      do i=1,nocc
          tmp1d=0d0
          Call DGemm("N","N",1,modims,modims,alpha,
     \               eigv(:,i),1,oneekin,modims,beta,tmp1d,1)
          Call DGemm("N","N",1,1,modims,2d0,
     \               eigv(:,i),1,tmp1d,modims,1d0,Ek,1)
      end do

      Ek=-Ek

      End Subroutine kinetic


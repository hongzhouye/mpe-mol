!     @@HY: added new variable 'data_path'

      Subroutine RunMPE(N,No,M,Alpha,data_path)
! ------------------------------------------------
! @@HY 08/24/18
!                 MPE CHEAT-SHEET
! General:
!     N      number of basis funcitons
!     No     number of electron pairs
!     M      number of auxiliary basis functions
!     h      1e integrals in MO basis (h = hkin + hnuc)
!     VVc    2e integrals in MO basis in chemists' notation
!     VVp    same thing, in physicists' notation
! SAH:
!     pair   N*No (could be either read from file "pairs.txt"
!            or computed from scratch)
!     alpha  degree of localization
! ------------------------------------------------
      Implicit None
      Character*200 data_path
      Integer*4 N,No,M,i,j,k,l,P,LineNum,ios,a,b,c,d
      Integer*4 FuncType
      Integer*4 maxrecs
      Real*8 h(N,N),hkin(N,N),hnuc(N,N),PQ(M,M),PQ_inv(M,M)
      Real*8 auxmomo(N,N,M),VVc(N,N,N,N),VVp(N,N,N,N)
      Real*8 pair(N,No),Alpha,Rho(M),Energy,tmp,denstwo(N,2)
      Real*8 total(N,No),EMPE,Ecorr1,Ecorr2
      Real*8 vx(M),vci(M),EDFA,ECI,vxsave(M,No)
      Real*8 EXX1(No),EXX2(No,No),ECI1(No),ECI2(No,No),DFTmo(N,N)
      Real*8 vGKS(M),EGKS,EGKS1(No),EGKS2(No,No)
      Real*8 Kxx
      Character(len=1) junk
      logical file_present

!     Read kin+nuc matrix in MO basis
      Open(Unit=12,File=trim(data_path) // '/oneeint.txt')
      Read(12,*)h
      Close(12)

!     Read kinetic matrix in MO basis
      Open(Unit=13,File=trim(data_path) // '/oneekin.txt')
      Read(13,*)hkin
      Close(13)

!     Read nuclear matrix in MO basis
      Open(Unit=14,File=trim(data_path) // '/oneenuc.txt')
      Read(14,*)hnuc
      Close(14)

!     Read MO coefficients
      Open(Unit=22,File=trim(data_path) // '/coeff.txt')
      Read(22,*)DFTMO
      Close(22)

!     Read (P | Q)
      Open(Unit=15,File=trim(data_path) // '/PQ.txt')
      Read(15,*)PQ
      Close(15)

!     Read (P | Q)^(-1)
      Open(Unit=16,File=trim(data_path) // '/PQ_inv.txt')
      Read(16,*)PQ_inv
      Close(16)

!     Read (ij | P), ij in MO basis
      Open(Unit=17,File=trim(data_path) // '/362.0.txt')
      Do i=1,N
        Do P=1,M
          Do j=1,N
            Read(17,*)auxmomo(i,j,P)
          End do
        End do
      End do
      Close(17)

!     Read two-electron integrals
!     First determine how many lines in 809, then read from 22
      maxrecs=100000000; LineNum=0
      Open(Unit=18, File=trim(data_path) // '/809.0.txt')
      do i=1,maxrecs
        read(18,*,IOSTAT=ios) junk
        if (ios/=0) Exit
        if (j==maxrecs) then
          write(*,*)'Error: Max number of lines exceeded'
          stop
        endif
        LineNum=LineNum+1
      end do
      Close(18)

!     Chemist's Notation
      Open(Unit=18,File=trim(data_path) // '/809.0.txt')
      Open(Unit=19,File=trim(data_path) // '/22.0.txt')
      VVc=0d0
      Do i=1,LineNum
        read(18,*)a,b,c,d
        if(i.gt.1)then
          read(19,*)VVc(a,b,c,d)
          VVc(b,a,c,d) = VVc(a,b,c,d)
          VVc(a,b,d,c) = VVc(a,b,c,d)
          VVc(b,a,d,c) = VVc(a,b,c,d)
          VVc(c,d,a,b) = VVc(a,b,c,d)
          VVc(d,c,a,b) = VVc(a,b,c,d)
          VVc(c,d,b,a) = VVc(a,b,c,d)
          VVc(d,c,b,a) = VVc(a,b,c,d)
        endif
      enddo
      Close(18)
      Close(19)

!     Physicist's Notation (used in FCI routine)
      Open(Unit=18,File=trim(data_path) // '/809.0.txt')
      Open(Unit=19,File=trim(data_path) // '/22.0.txt')
      VVp=0d0
      Do i=1,LineNum
        read(18,*)a,b,c,d
        if(i.gt.1)then
          read(19,*)tmp
          VVp(a,c,b,d) = tmp
          VVp(c,a,d,b) = tmp
          VVp(b,c,a,d) = tmp
          VVp(c,b,d,a) = tmp
          VVp(a,d,b,c) = tmp
          VVp(d,a,c,b) = tmp
          VVp(b,d,a,c) = tmp
          VVp(d,b,c,a) = tmp
        endif
      enddo
      Close(18)
      Close(19)

!     Do Self-attractive Hartree decomposition first
      write(6,*)'**********************************'
      write(6,*)'-------SAH decomposition-------'

      inquire(FILE=trim(data_path) // '/pairs.txt',EXIST=file_present) !Check if file exists
      if (file_present) then
        open(unit=21,File=trim(data_path) // '/pairs.txt')
        read(21,*) alpha
        write(6,*) 'Reading pair densities from a file. Alpha= ', alpha
        read(21,*) pair
        close(21)
      else
        Write(6,'(A,F7.3)') "Running SAH with alpha = ", Alpha
        pair=0d0
! @@HY: using the HF occupied orbitals as initial guess?
        Do i=1,No
            pair(i,i)=1d0
        End do
        Call partition(No,N,M,PQ,PQ_inv,auxmomo,h,hkin,alpha,pair)
        open(unit=21,File=trim(data_path) // '/pairs.txt')
          write(21,*) alpha
          write(21,*) pair
        close(21)
      end if

      Open(Unit=21, File=trim(data_path) // "/xc_func.txt")
      Read(21,*) FuncType, Kxx
      Close(21)

!     DFA-OEP for total density
      write(6,*)'-----DFA-OEP for total density-----'
      vx=0d0; total=0d0 
      Do i=1,No
        total(i,i)=1d0
      End do
      Call oep(data_path,total,vx,No,N,M,auxmomo,h,hkin,hnuc,VVc,
     \         DFTMO,FuncType,Kxx,EDFA)
      write(6,*)'#EDFA',EDFA
      EMPE=EDFA

!     DFA-OEP for one pair
      write(6,*)'-------MPE1-------'
      Do i=1,No
        vx=0d0
        Call oep(data_path,pair(:,i),vx,1,N,M,auxmomo,h,hkin,hnuc,VVc,
     \           DFTMO,FuncType,Kxx,EXX1(i))
        write(6,*)'#pair,EDFA1',i,EXX1(i)

        vGKS=vx; Rho=0d0
!        vGKS=0.d0
!       Note Rho is one-pair density
        do j=1,N; do k=1,N
            Rho=Rho+2d0*pair(j,i)*pair(k,i)*auxmomo(j,k,:)
        end do; end do
        Call GKSenergy(N,1,M,vGKS,h,VVc,auxmomo,Rho,EGKS1(i))
        write(6,*)'#pair,EGKS1',i,EGKS1(i)

        vci=vGKS; Rho=0d0
!        vci=vx
!       Note Rho is one-electron density
        do j=1,N; do k=1,N
            Rho=Rho+pair(j,i)*pair(k,i)*auxmomo(j,k,:)
        end do; end do
        Call cienergy(N,1,M,vci,h,VVp,auxmomo,Rho,ECI1(i))
        write(6,*)'#pair,ECI1,Ecorr1',i,ECI1(i),ECI1(i)-EXX1(i)
        EMPE=EMPE+ECI1(i)-EXX1(i)
        write(6,*)'------------------'

      End do
      write(6,*)'###EDFA,EMPE1',EDFA,EMPE

!     EXX-OEP for two pairs
! @@HY: hacky way to skip MPE2
      Go to 137

      write(6,*)'-------MPE2-------'
      Do i=1,No; do j=i+1,No
        vx=0d0
        denstwo(:,1)=pair(:,i); denstwo(:,2)=pair(:,j)
        Call oep(denstwo,vx,2,N,M,auxmomo,h,hkin,hnuc,VVc,
     \           DFTMO,EXX2(i,j)) 
        write(6,*)'#pair i j,EDFA2',i,j,EXX2(i,j)

        vGKS=vx; Rho=0d0
!        vGKS=0d0
        do l=1,N; do k=1,N
            Rho=Rho+2d0*pair(l,i)*pair(k,i)*auxmomo(l,k,:)
     \             +2d0*pair(l,j)*pair(k,j)*auxmomo(l,k,:)
        end do; end do
        Call GKSenergy(N,2,M,vGKS,h,VVc,auxmomo,Rho,EGKS2(i,j))
        write(6,*)'#pair i j,EGKS2',i,j,EGKS2(i,j)

        vci=vGKS; Rho=0d0
!        vci=vx
        do l=1,N; do k=1,N
            Rho=Rho+pair(l,i)*pair(k,i)*auxmomo(l,k,:)
     \             +pair(l,j)*pair(k,j)*auxmomo(l,k,:)
        end do; end do
        Call cienergy(N,2,M,vci,h,VVp,auxmomo,Rho,ECI2(i,j))
        Ecorr2=ECI2(i,j)-EXX2(i,j)-(ECI1(i)-EXX1(i))
     \            -(ECI1(j)-EXX1(j))
        write(6,*)'#pair i j,ECI2,Ecorr2',i,j,ECI2(i,j),Ecorr2
        EMPE=EMPE+Ecorr2
        write(6,*)'------------------'
      End do; end do
      write(6,*)'###EMPE2',EMPE

!     FCI-OEP for total density
      write(6,*)'-----FCI-OEP for total density-----'
      Rho=0d0
      Do i=1,No
        Rho=Rho+auxmomo(i,i,:)
      end do
      vci=0d0; ECI=0d0
!      Call cienergy(N,No,M,vci,h,VVp,auxmomo,Rho,ECI)
      write(6,*)'###FCI-OEP',ECI
      write(6,*)'**********************************'

137   Continue

      End subroutine


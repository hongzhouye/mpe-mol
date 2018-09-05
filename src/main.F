      Implicit None
      Integer*4 N,No,M,pair_num
      Real*8 Alpha
      Character*200 arg_buf,data_path,grid_path

      If(Iargc() .lt. 3) then
        Write(6,'(A)') "Usage: alpha data_path grid_path [pair_num]"
        Write(6, *)
        Write(6,'(A)') "  Arg4 (pair_num) is optional:"
        Write(6,'(A)') "  Case1: not give --> perform normal MPE1 for"
     \    // " all pairs sequentially."
        Write(6,'(A)') "  Case2: positive --> perform normal MPE1 for"
     \    // " that pair only."
        Write(6,'(A)') "  Case3: negative --> perform KS-OEP for that"
     \    // " pair only (skip GKS/FCI)."
        Write(6,'(A)') "  Case4: 0        --> perform KS-OEP for all"
     \    // " pairs (skip GKS/FCI)."
        Call EXIT(1)
      End If

      Call GetArg(1, arg_buf)
      Read(arg_buf,*) Alpha
      Call GetArg(2, data_path)
      Call GetArg(3, grid_path)
      If(Iargc() .gt. 3) Then
        Call GetArg(4, arg_buf)
        Read(arg_buf,*) pair_num
      Else
        pair_num = 0
      End If

      Open(Unit=11,File=trim(data_path) // '/dimensions.txt')
      Read(11,*)N,No,M
      Close(11)

      Call RunMPE(N,No,M,Alpha,data_path,grid_path,pair_num)

      End

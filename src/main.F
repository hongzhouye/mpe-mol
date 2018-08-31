      Implicit None
      Integer*4 N,No,M
      Real*8 Alpha
      Character*200 arg_buf,data_path,grid_path

      If(Iargc() .lt. 3) then
        Write(6,'(A)') "Usage: alpha data_path grid_path"
        Call EXIT(1)
      End If

      Call GetArg(1, arg_buf)
      Read(arg_buf,*) Alpha
      Call GetArg(2, data_path)
      Call GetArg(2, grid_path)

      Open(Unit=11,File=trim(data_path) // '/dimensions.txt')
      Read(11,*)N,No,M
      Close(11)

      Call RunMPE(N,No,M,Alpha,data_path,grid_path)

      End

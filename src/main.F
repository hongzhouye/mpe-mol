      Implicit None
      Integer*4 N,No,M
      Real*8 Alpha
      Character*200 data_path,arg_buf

      If(Iargc() .lt. 2) then
        Write(*,'(A)') "Usage: alpha data_path"
        Call EXIT(1)
      End If

      Call GetArg(1, arg_buf)
      Read(arg_buf,*) Alpha
      Call GetArg(2, data_path)

      Open(Unit=11,File=trim(data_path) // '/dimensions.txt')
      Read(11,*)N,No,M
      Close(11)

      Call RunMPE(N,No,M,Alpha,data_path)

      End

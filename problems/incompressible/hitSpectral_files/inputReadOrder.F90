read(fileInput,*) 
read(fileInput,*) 
read(fileInput,*) 
read(fileInput,*) 
read(fileInput,*) 
read(fileInput,*) 
read(fileInput,*    ,ERR=2000) Nx,Ny,Nz
read(fileInput,*)
read(fileInput,*    ,ERR=2000) RESTART
read(fileInput,*)
read(fileInput,*    ,ERR=2000) useHit3dinit
read(fileInput,*)
read(fileInput,'(A)',ERR=2000) HIT3dinitDir
read(fileInput,*)
read(fileInput,'(A)',ERR=2000) RESTARTDir
read(fileInput,*)
read(fileInput,*    ,ERR=2000) NT
read(fileInput,*)
read(fileInput,*    ,ERR=2000) useCFL
read(fileInput,*)
read(fileInput,*    ,ERR=2000) CFL
read(fileInput,*)
read(fileInput,*    ,ERR=2000) DT
read(fileInput,*)
read(fileInput,*    ,ERR=2000) REY
read(fileInput,*)
read(fileInput,*    ,ERR=2000) Dealias
read(fileInput,*)
read(fileInput,*    ,ERR=2000) TSTEP_DUMP
read(fileInput,*)
read(fileInput,'(A)',ERR=2000) OutputDir
read(fileInput,*)
read(fileInput,*    ,ERR=2000) TSTEP_RESTART
read(fileInput,*)
read(fileInput,*    ,ERR=2000) FFT_PLAN_EXHAUSTIVE
read(fileInput,*)
read(fileInput,*    ,ERR=2000) RunIDX
read(fileInput,*)
read(fileInput,'(A)',ERR=2000) eof

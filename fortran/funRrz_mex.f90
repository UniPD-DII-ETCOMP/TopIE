!---- never change ------

#include "fintrf.h"      

subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    implicit none

! mwPointer mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
!---- end never change ------
	  
	  
! mwSize stuff for mexing
      mwSize mo,no,siz
	  
! mwPointer stuff for mexing
      mwPointer mxCreateDoubleMatrix
      mwPointer mxGetPr
	  mwPointer N_node_pr,N_face_pr,rho_pr,N_GAUSS_pr,Matrix_P0_pr,F1_pr
	  mwPointer N_edge_pr,C1_pr,Fac_Ed_loc_pr,Aed_pr,Led_pr,N_thread_pr,R_pr
      mwPointer m, n
      mwPointer mxGetM, mxGetN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!integer (normal not 4/8)
      integer flag,i,j,ii,jj,icount
      integer mxIsNumeric 
      integer*4 ComplexFlag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
! fortran subroutine arguments
	  real*8,allocatable,dimension(:,:) :: Matrix_P0, F1_r, R, C1_r, Fac_Ed_loc_r
      integer*8,allocatable,dimension(:,:) :: F1,C1, Fac_Ed_loc
	  real*8,allocatable,dimension(:) ::  Aed, Led, rho
	  !integer*8,allocatable,dimension(:) :: 
!	  integer*8,allocatable,dimension(:) :: 
	  real*8 :: N_node_r, N_face_r, N_edge_r, N_thread_r, N_GAUSS_r
      integer*8 :: N_node, N_face, N_edge, N_thread, N_GAUSS 
	  character*80 msg
      logical debu
       
      debu = .true. ! .true. o .false. per attivare o disattivare il debug
	  if(debu) open(unit=66,file='logRrz_axi.txt',status='unknown')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
if (nrhs .ne. 12) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_auto:nInput', &
                                '12 input arguments required')
elseif (nlhs .ne. 1) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_auto:nOutput', &
                                '1 output argument required')
endif	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Check to see INPUTS are numeric.
	  do ii = 1,12
        if (mxIsNumeric(prhs(ii)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:L_axi_vv:NonNumeric', &
                                'Inputs must be numeric.')
        endif
	  enddo
	  if(debu) write(66,*) 'check numeric done'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Check that input #1 is INTEGER and fetch it
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 1 must be scalar, N_node.')
      endif	  
      siz = m*n
      N_node_pr = mxGetPr(prhs(1))
      call mxCopyPtrToReal8(N_node_pr, N_node_r, siz) 
      N_node=int(N_node_r,8) 
	  if(debu) write(66,*) 'input 1' 	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Check that input #2 is INTEGER and fetch it
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 2 must be scalar, N_face.')
      endif	  
      siz = m*n
      N_face_pr = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(N_face_pr, N_face_r, siz) 
      N_face=int(N_face_r,8) 
	  if(debu) write(66,*) 'input 2' 	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!     Check that input #3 is INTEGER vector and fetch it
      m = mxGetM(prhs(3))
      n = mxGetN(prhs(3))
      if(m .ne. 1 .or. n .ne. N_face) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 3: rho')
      endif	  
      siz = m*n
      rho_pr = mxGetPr(prhs(3))
      allocate(rho(N_face))
      call mxCopyPtrToReal8(rho_pr, rho, siz)  
	  if(debu) write(66,*) 'input 3'	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Check that input #4 is INTEGER and fetch it
      m = mxGetM(prhs(4))
      n = mxGetN(prhs(4))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 4 must be scalar, N_GAUSS.')
      endif	  
      siz = m*n
      N_GAUSS_pr = mxGetPr(prhs(4))
      call mxCopyPtrToReal8(N_GAUSS_pr, N_GAUSS_r, siz) 
      N_GAUSS=int(N_GAUSS_r,8) 
	  if(debu) write(66,*) 'input 4'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	       
!	  Check that input #5 is REAL matrix and fetch it
      m = mxGetM(prhs(5))
      n = mxGetN(prhs(5))
	  if(debu) write(66,*) m
	  if(debu) write(66,*) n
      if(m .ne. N_node .or. n .ne. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 5: Matrix_P0')
      endif	  
      siz = m*n
      Matrix_P0_pr = mxGetPr(prhs(5))
	  allocate(Matrix_P0(N_node,3))
      call mxCopyPtrToReal8(Matrix_P0_pr, Matrix_P0, siz) 
	  if(debu) write(66,*) 'input 5'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
!     Check that input #6 is INTEGER matrix and fetch it
      m = mxGetM(prhs(6))
      n = mxGetN(prhs(6))
      if(m .ne. 4 .or. n .ne. N_face) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 6: F1')
      endif	  
      siz = m*n
      F1_pr = mxGetPr(prhs(6))
      allocate(F1_r(4,N_face))
      allocate(F1(4,N_face))
      call mxCopyPtrToReal8(F1_pr, F1_r, siz) 
      F1 = int(F1_r,8)
      deallocate(F1_r)
      if(debu) write(66,*) 'input 6'	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Check that input #7 is INTEGER and fetch it
      m = mxGetM(prhs(7))
      n = mxGetN(prhs(7))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 7 must be scalar, N_edge.')
      endif	  
      siz = m*n
      N_edge_pr = mxGetPr(prhs(7))
      call mxCopyPtrToReal8(N_edge_pr, N_edge_r, siz) 
      N_edge=int(N_edge_r,8) 
	  if(debu) write(66,*) 'input 7'   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
!     Check that input #8 is INTEGER matrix and fetch it
      m = mxGetM(prhs(8))
      n = mxGetN(prhs(8))
      if(m .ne. 4 .or. n .ne. N_face) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 8: C1')
      endif	  
      siz = m*n
      C1_pr = mxGetPr(prhs(8))
      allocate(C1_r(4,N_face))
      allocate(C1(4,N_face))
      call mxCopyPtrToReal8(C1_pr, C1_r, siz) ! da double precision a reale  
      C1 = int(C1_r,8)
      deallocate(C1_r)
      if(debu) write(66,*) 'input 8'	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
!     Check that input #9 is INTEGER matrix and fetch it
      m = mxGetM(prhs(9))
      n = mxGetN(prhs(9))
      if(m .ne. 4 .or. n .ne. N_face) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 9: Fac_Ed_loc')
      endif	  
      siz = m*n
      Fac_Ed_loc_pr = mxGetPr(prhs(9))
      allocate(Fac_Ed_loc_r(4,N_face))
      allocate(Fac_Ed_loc(4,N_face))
      call mxCopyPtrToReal8(Fac_Ed_loc_pr, Fac_Ed_loc_r, siz)  
      Fac_Ed_loc = int(Fac_Ed_loc_r,8)
      deallocate(Fac_Ed_loc_r)
      if(debu) write(66,*) 'input 9'		  	 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!     Check that input #10 is INTEGER vector and fetch it
      m = mxGetM(prhs(10))
      n = mxGetN(prhs(10))
      if(m .ne. 1 .or. n .ne. N_edge) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 10: Aed')
      endif	  
      siz = m*n
      Aed_pr = mxGetPr(prhs(10))
      allocate(Aed(N_edge))
      call mxCopyPtrToReal8(Aed_pr, Aed, siz) 
	  if(debu) write(66,*) 'input 10'	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!     Check that input #11 is INTEGER vector and fetch it
      m = mxGetM(prhs(11))
      n = mxGetN(prhs(11))
      if(m .ne. 1 .or. n .ne. N_edge) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 11: Led')
      endif	  
      siz = m*n
      Led_pr = mxGetPr(prhs(11))
      allocate(Led(N_edge))
      call mxCopyPtrToReal8(Led_pr, Led, siz)  
	  if(debu) write(66,*) 'input 11'	 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Check that input #12 is INTEGER and fetch it
      m = mxGetM(prhs(12))
      n = mxGetN(prhs(12))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 12 must be scalar, N_thread.')
      endif	  
      siz = m*n
      N_thread_pr = mxGetPr(prhs(12))
      call mxCopyPtrToReal8(N_thread_pr, N_thread_r, siz)
      N_thread=int(N_thread_r,8) 
	  if(debu) write(66,*) 'input 12'   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      	  	  
! call the computational subroutine
      allocate(R(N_edge,N_edge))
      call funRrz(N_node,N_face,rho,N_GAUSS,Matrix_P0,F1,N_edge,C1,Fac_Ed_loc,Aed,Led,N_thread,R)
       if(debu) write(66,*) 'created matrices'
      deallocate(rho,Matrix_P0,F1,C1,Fac_Ed_loc,Aed,Led)
	  if(debu) write(66,*) 'deallocated unuseful stuffs'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
! Create a matrix for the return argument 1
      mo=N_edge
      no=N_edge
	  ComplexFlag = 0
      plhs(1) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
! Load the output 1 into a MATLAB array.
      R_pr = mxGetPr(plhs(1))
      siz=mo*no
      call mxCopyReal8ToPtr(R, R_pr, siz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      if(debu) write(66,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      deallocate(R)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      if(debu) write(66,*) 'close logRrz_axi.txt'
      if(debu) close(66)
      return
      end


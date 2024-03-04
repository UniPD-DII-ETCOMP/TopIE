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
	  mwPointer N_node_pr,N_edge_pr,N_GAUSS1_pr,N_GAUSS2_pr,N_thread_pr,Matrix_P0_pr,G1_pr,P_pr,Aed_pr
      mwPointer m, n
      mwPointer mxGetM, mxGetN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!integer (normal not 4/8)
      integer flag,i,j,ii,jj,icount
      integer mxIsNumeric 
      integer*4 ComplexFlag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
! fortran subroutine arguments
	  real*8,allocatable,dimension(:,:) :: Matrix_P0, G1_r, P
      integer*8,allocatable,dimension(:,:) :: G1
	  real*8,allocatable,dimension(:) :: Aed
!	  integer*8,allocatable,dimension(:) :: 
	  real*8 :: N_node_r, N_edge_r, N_GAUSS1_r, N_GAUSS2_r, N_thread_r
      integer*8 :: N_node, N_edge, N_GAUSS1, N_GAUSS2, N_thread
	  character*80 msg
      logical debu
       
      debu = .true. ! .true. o .false. per attivare o disattivare il debug
	  if(debu) open(unit=66,file='logL_axi_ss.txt',status='unknown')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
if (nrhs .ne. 8) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_auto:nInput', &
                                '8 input arguments required')
elseif (nlhs .ne. 1) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_auto:nOutput', &
                                '1 output argument required')
endif	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Check to see INPUTS are numeric.
	  do ii = 1,8
        if (mxIsNumeric(prhs(ii)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:L_axi_vv:NonNumeric', &
                                'Inputs must be numeric.')
        endif
	  enddo
	  if(debu) write(66,*) 'sono arrivato al check numerico'
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
      call mxCopyPtrToReal8(N_node_pr, N_node_r, siz) ! da double precision a reale
      N_node=int(N_node_r,8) ! da reale a intero
	  if(debu) write(66,*) 'input 1' 	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Check that input #2 is INTEGER and fetch it
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 2 must be scalar, N_edge.')
      endif	  
      siz = m*n
      N_edge_pr = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(N_edge_pr, N_edge_r, siz) ! da double precision a reale
      N_edge=int(N_edge_r,8) ! da reale a intero
	  if(debu) write(66,*) 'input 2' 	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Check that input #3 is INTEGER and fetch it
      m = mxGetM(prhs(3))
      n = mxGetN(prhs(3))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 3 must be scalar, N_GAUSS1.')
      endif	  
      siz = m*n
      N_GAUSS1_pr = mxGetPr(prhs(3))
      call mxCopyPtrToReal8(N_GAUSS1_pr, N_GAUSS1_r, siz) ! da double precision a reale
      N_GAUSS1=int(N_GAUSS1_r,8) ! da reale a intero
	  if(debu) write(66,*) 'input 3' 	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Check that input #4 is INTEGER and fetch it
      m = mxGetM(prhs(4))
      n = mxGetN(prhs(4))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 4 must be scalar, N_GAUSS2.')
      endif	  
      siz = m*n
      N_GAUSS2_pr = mxGetPr(prhs(4))
      call mxCopyPtrToReal8(N_GAUSS2_pr, N_GAUSS2_r, siz) ! da double precision a reale
      N_GAUSS2=int(N_GAUSS2_r,8) ! da reale a intero
	  if(debu) write(66,*) 'input 4' 	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Check that input #5 is INTEGER and fetch it
      m = mxGetM(prhs(5))
      n = mxGetN(prhs(5))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 5 must be scalar, N_thread.')
      endif	  
      siz = m*n
      N_thread_pr = mxGetPr(prhs(5))
      call mxCopyPtrToReal8(N_thread_pr, N_thread_r, siz) ! da double precision a reale
      N_thread=int(N_thread_r,8) ! da reale a intero
	  if(debu) write(66,*) 'input 5' 	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	       
!	  Check that input #6 is REAL matrix and fetch it
      m = mxGetM(prhs(6))
      n = mxGetN(prhs(6))
	  if(debu) write(66,*) m
	  if(debu) write(66,*) n
      if(m .ne. N_node .or. n .ne. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 6: Matrix_P0')
      endif	  
      siz = m*n
      Matrix_P0_pr = mxGetPr(prhs(6))
	  allocate(Matrix_P0(N_node,3))
      call mxCopyPtrToReal8(Matrix_P0_pr, Matrix_P0, siz) ! da double precision a reale
	  if(debu) write(66,*) 'input 6'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
!     Check that input #7 is INTEGER matrix and fetch it
      m = mxGetM(prhs(7))
      n = mxGetN(prhs(7))
      if(m .ne. 2 .or. n .ne. N_edge) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 7: F1')
      endif	  
      siz = m*n
      G1_pr = mxGetPr(prhs(7))
	  allocate(G1_r(2,N_edge))
      allocate(G1(2,N_edge))
      call mxCopyPtrToReal8(G1_pr, G1_r, siz) ! da double precision a reale  
	  G1 = int(G1_r,8)
	  deallocate(G1_r)
	  if(debu) write(66,*) 'input 7'	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
!     Check that input #8 is INTEGER matrix and fetch it
      m = mxGetM(prhs(8))
      n = mxGetN(prhs(8))
      if(m .ne. 1 .or. n .ne. N_edge) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 8: Aed')
      endif	  
      siz = m*n
      Aed_pr = mxGetPr(prhs(8))
      allocate(Aed(N_edge))
      call mxCopyPtrToReal8(Aed_pr, Aed, siz) ! da double precision a reale  
	  if(debu) write(66,*) 'input 8'	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      	  	  
! call the computational subroutine.
      allocate(P(N_edge,N_edge))
      call funPphiphi3_ss(N_node,N_edge,N_thread,N_GAUSS1,N_GAUSS2,Matrix_P0,G1,Aed,P)
       if(debu) write(66,*) 'ho chiamato la subroutine e creato le matrici'
      deallocate(Matrix_P0,G1,Aed)
	  if(debu) write(66,*) 'ho deallocato le cose inutili'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
! Create a matrix for the return argument 1
      mo=N_edge
      no=N_edge
	  ComplexFlag = 0
      plhs(1) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
! Load the output 1 into a MATLAB array.
      P_pr = mxGetPr(plhs(1))
      siz=mo*no
      call mxCopyReal8ToPtr(P, P_pr, siz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      if(debu) write(66,*) 'ho convertito la matrice per matlab'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      deallocate(P)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      if(debu) write(66,*) 'ho deallocato la matrice e ora chiudo il logPaxi_ss.txt'
      if(debu) close(66)
      return
      end


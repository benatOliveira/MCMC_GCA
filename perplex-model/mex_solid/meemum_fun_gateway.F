#include <fintrf.h> 
       SUBROUTINE mexFunction(nlhs, plhs, nrhs, prhs)
!!!pointers to input data: prhs and output: plhs
       mwPointer plhs(*), prhs(*)
!!!number of inputs on right and left hand side
       INTEGER nlhs, nrhs
!!!pointers to intrinsic functions of the matlab library
       mwPointer mxgetpr,mxgetscalar!!, mxCreateString
       mwPointer mxgetm, mxgetn, mxCreateDoubleMatrix
       mwPointer mxCreateCharMatrixFromStrings
C     Array information:
       mwSize v_size, cblk_size !!, pname_size, props_size, psys_size

C     Arguments for computational routine:
       include 'perplex_parameters.h'
       double precision cblk_i(k5)!! k5 = 12
       double precision v_i(l2)!!l2 = 5
       double precision props_o(k5,12+k5) 
       double precision psys_o(1,10+k5) 
c       character*14 names(k5)
       character*14 pname_o(k5)
       double precision iPs
       integer props_r, props_c, psys_r, psys_c!!, pname_r, pname_c 
       integer props_s, psys_s, iPs_s!!, pname_s 

C     Get the sizes of the input arrays.
       v_row = mxgetm(prhs(1))
       v_col = mxgetn(prhs(1))
       v_size = v_row*v_col
       CALL mxCopyPtrToReal8(mxGetPr(prhs(1)), v_i, v_size)

       cblk_row = mxgetm(prhs(2))
       cblk_col = mxgetn(prhs(2))
       cblk_size = cblk_row*cblk_col
       CALL mxCopyPtrToReal8(mxGetPr(prhs(2)), cblk_i, cblk_size)

c Call the FORTRAN routine
       CALL meemum_fun(pname_o,props_o,psys_o,iPs,v_i,cblk_i)       
c Creating null matrices as default output
c       pname_r = 1
c       pname_c = k5
c       pname_s = pname_r*pname_c
       plhs(1) = mxCreateCharMatrixFromStrings(k5, pname_o)

       props_r = size(props_o,1)
       props_c = size(props_o,2)
       props_s = props_r*props_c
       plhs(2)= mxCreateDoubleMatrix(props_r, props_c, 0) 
      
       psys_r = size(psys_o,1)
       psys_c = size(psys_o,2)     
       psys_s = psys_r*psys_c
       plhs(3)= mxCreateDoubleMatrix(psys_r, psys_c, 0)
 
       plhs(4) = mxCreateDoubleMatrix(1, 1, 0)
c Copy output data from local variables to matlab output variables
c       CALL mxCopyCharacterToPtr(pname_o,mxGetData(plhs(1)),pname_s)       
       CALL mxCopyReal8ToPtr(props_o, mxGetPr(plhs(2)), props_s)
       CALL mxCopyReal8ToPtr(psys_o, mxGetPr(plhs(3)), psys_s)
       CALL mxCopyReal8ToPtr(iPs, mxGetPr(plhs(4)), 1)
       return
       end

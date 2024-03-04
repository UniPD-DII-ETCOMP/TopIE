function [Matrix_Rm] = fun_construct_Rm(param,N,Matrix_P0_m,...
                                        F1_m,C1_m,Fac_Ed_loc_m,...
                                        Aed_m,Led_m,...
                                        NGAUSS,nThread)
try
    Matrix_Rm=funRrz_for(N.node_mag,N.face_mag,param',...
                        NGAUSS,Matrix_P0_m,F1_m,...
                        N.edge_mag,C1_m,Fac_Ed_loc_m,Aed_m.',Led_m.',...
                        nThread);
    Matrix_Rm=sparse(Matrix_Rm);
catch ME
    warning(' FORTAN MEX-FILE is not supported, try to re-mex it in fortran/make_Rrz')
    warning([' MESSAGE ERROR:' ME.message]);
    return
end 
end
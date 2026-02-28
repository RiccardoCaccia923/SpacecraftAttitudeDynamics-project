function N = nonOrtMat(sigma_no)
    
    n_xy = sigma_no * randn();
    n_xz = sigma_no * randn();
    n_yx = sigma_no * randn();
    n_yz = sigma_no * randn();
    n_zx = sigma_no * randn();
    n_zy = sigma_no * randn();
    
    N = [ 1    n_xy n_xz;
          n_yx 1    n_yz;
          n_zx n_zy 1   ];
end
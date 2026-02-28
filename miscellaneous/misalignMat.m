function Rmisalign = misalignMat(sigma_theta)
e_x = sigma_theta*randn();
e_y = sigma_theta*randn();
e_z = sigma_theta*randn();

Rx = [1 0 0;
      0 cos(e_x) -sin(e_x);
      0 sin(e_x)  cos(e_x)];

Ry = [cos(e_y) 0 sin(e_y);
      0        1 0;
     -sin(e_y) 0 cos(e_y)];

Rz = [cos(e_z) -sin(e_z) 0;
      sin(e_z)  cos(e_z) 0;
      0         0        1];
  
Rmisalign = Rz*Ry*Rx;
end
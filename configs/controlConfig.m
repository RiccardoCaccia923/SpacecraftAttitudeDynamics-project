function control = controlConfig(data)
    
    Jx = data.Imatrix(1,1);
    Jy = data.Imatrix(2,2);
    Jz = data.Imatrix(3,3);

    % Wheel Parameters
    wc_zWheel = 0.02;     
    PM_zWheel = 80;       
    
    % Z axis (Yaw - Wheel)
    [control.Kp_z, control.Kd_z] = PDcontrollerGains(Jz, wc_zWheel, PM_zWheel);

    % Magneorquers Parameters
    wc_xyMag = 0.03;      
    PM_xyMag = 80;        

    % X axis (Roll - Magnetorquers)
    [control.Kp_x, control.Kd_x] = PDcontrollerGains(Jx, wc_xyMag, PM_xyMag);
    
    
    % Y axis (Pitch - Magnetorques)
    [control.Kp_y, control.Kd_y] = PDcontrollerGains(Jy, wc_xyMag, PM_xyMag);
    control.Kp_y=control.Kp_y/4;
    control.Kp_x=control.Kp_x/4;
    control.Kp_z=control.Kp_z/2;

    control.Kdet = -(4*pi/data.T)*(1+sin(deg2rad(data.i)))*min([Jx,Jy,Jz])*100000;
end
function si = gen_si(info)
    sys = info.sys;
    Phi = blkdiag(info.T*info.epsilon^2*eye(sys.nx), -eye(info.T));
    data_matrix = [eye(sys.nx)             info.Xp;...
                   zeros(sys.nx,sys.nx)   -info.Xm;...
                   zeros(sys.nu,sys.nx)   -info.Um];
    si = data_matrix*Phi*data_matrix';
end
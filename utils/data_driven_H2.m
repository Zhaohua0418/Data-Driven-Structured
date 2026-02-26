function out = data_driven_H2(info,pattern,SI)
    sys = info.sys;
%     SI = gen_si(info);
    lambda = 1;
    prev_X = eye(sys.nx);
    while 1
        cvx_begin sdp
            variables X(sys.nx,sys.nx) Y(sys.nx,sys.nx) K(sys.nu,sys.nx) W(sys.ny,sys.ny) Z_plus(sys.nx,sys.nx) a b
            minimize trace(W) + lambda * trace(Z_plus)
            subject to
                    [X-sys.G*sys.G'-b*eye(sys.nx)   zeros(sys.nx,sys.nx)    zeros(sys.nx,sys.nu)    zeros(sys.nx,sys.nx);...
                     zeros(sys.nx,sys.nx)'          zeros(sys.nx,sys.nx)    zeros(sys.nx,sys.nu)    eye(sys.nx,sys.nx);...
                     zeros(sys.nx,sys.nu)'          zeros(sys.nx,sys.nu)'   zeros(sys.nu,sys.nu)    K;...
                     zeros(sys.nx,sys.nx)'          eye(sys.nx,sys.nx)      K'                      Y] - a*blkdiag(SI,zeros(sys.nx,sys.nx))>=0;
    
                    [W          sys.C+sys.D*K;...
                     (sys.C+sys.D*K)'   Y]>=0;
    
                    Y <= inv(prev_X) - inv(prev_X)*(X-prev_X)*inv(prev_X) + Z_plus;
    
                    a>=0; b>0;
                    X>=0; W>=0; Z_plus >= 0;
                    
                    K.*(1-pattern) == 0;
        cvx_end
        h2 = sqrt(trace(W));
        if norm(X-prev_X,"fro") < 0.01 && norm(X-inv(Y),"fro") < 0.01
             out = struct;
             out.X = X;
             out.Y = Y;
             out.K = K;
             out.W = W;
             out.Z_plus = Z_plus;
             out.a = a;
             out.b = b;
             out.h2 = h2;
             break
        end
        lambda = lambda * 2;
        prev_X = X;
    end
end
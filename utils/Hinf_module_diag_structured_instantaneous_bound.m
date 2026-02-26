classdef Hinf_module_diag_structured_instantaneous_bound
    
    properties
        A,B,G,C,D,H,nx,nu,ny,nd,pattern
    end
    
    methods
        function obj = Hinf_module_diag_structured_instantaneous_bound(sys,pattern)
            obj.A = sys.A;
            obj.B = sys.B;
            obj.G = sys.G;
            obj.C = sys.C;
            obj.D = sys.D;
            obj.H = sys.H;
            obj.nx = sys.nx;
            obj.nu = sys.nu;
            obj.ny = sys.ny;
            obj.nd = sys.nd;
            obj.pattern = pattern;
        end
        
        function out = model_based_Hinf(obj)
            cvx_begin sdp
                variables P(obj.nx,obj.nx) L(obj.nu,obj.nx) hinf_square
                minimize hinf_square
                subject to
                        [P                      obj.A*P+obj.B*L             obj.G                   zeros(obj.nx,obj.ny);...
                         (obj.A*P+obj.B*L)'     P                           zeros(obj.nx,obj.nd)    (obj.C*P+obj.D*L)';...
                         obj.G'                 zeros(obj.nd,obj.nx)        eye(obj.nd)             obj.H';...
                         zeros(obj.ny,obj.nx)   obj.C*P+obj.D*L             obj.H                   hinf_square*eye(obj.ny)]>0;
                        P.*(1-eye(obj.nx,obj.nx)) == 0;
                        L.*(1-obj.pattern) == 0;
            cvx_end
            hinf = sqrt(hinf_square);
            K = L*inv(P);
            out = struct;
            out.K = K;
            out.P = P;
            out.L = L;
            out.hinf = hinf;
        end

        function out = data_driven_Hinf(obj,info)
            SI = gen_si(info);
            cvx_begin sdp
                variables P(obj.nx,obj.nx) L(obj.nu,obj.nx) R(obj.nu+obj.nx,obj.nu+obj.nx) tao(1,info.T) b hinf_square
                minimize hinf_square
                subject to
                        sum_over = 0;
                        for t = 1:info.T
                            left = [eye(obj.nx)            info.Xp(:,t);...
                                  zeros(obj.nx)           -info.Xm(:,t);...
                                  zeros(obj.nu,obj.nx)    -info.Um(:,t);...
                                  zeros(obj.nx)            zeros(obj.nx,1)];
                            mid = blkdiag(info.epsilon^2*eye(obj.nx),-1);
                            sum_over = sum_over + tao(1,t) * left * mid * left';
                        end

                        upper_left = [P-obj.G*obj.G'-b*eye(obj.nx)   zeros(obj.nx,obj.nu+obj.nx);...
                                      zeros(obj.nu+obj.nx,obj.nx)       -R];
                        upper_right = [obj.G*obj.H';...
                                       R*[obj.C';...
                                          obj.D']];
                        bottom_right = [hinf_square*eye(obj.ny) - [obj.C obj.D]*R*[obj.C';...
                                                                                   obj.D'] - obj.H*obj.H'];
                        [upper_left     upper_right;...
                         upper_right'   bottom_right] - sum_over >=0;

                        upper_left_2 = R-[P L';...
                                          L zeros(obj.nu)];
                        upper_right_2 = [zeros(obj.nx);...
                                         L];
                        [upper_left_2       upper_right_2;...
                         upper_right_2'     P]>=0;
                        P.*(1-eye(obj.nx,obj.nx)) == 0;
                        L.*(1-obj.pattern) == 0;
                        tao>=0;
                        b>0;
                        P>=0; R>=0;
            cvx_end
            hinf = sqrt(hinf_square);
            K = L*inv(P);
            out = struct;
            out.P = P;
            out.L = L;
            out.K = K;
            out.R = R;
            out.tao = tao;
            out.b = b;
            out.hinf = hinf;
        end
    end
end


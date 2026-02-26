classdef Hinf_module_my_unstructured
    properties
        A,B,G,C,D,H,nx,nu,ny,nd
    end
    
    methods
        function obj = Hinf_module_my_unstructured(sys)
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
            cvx_begin sdp
                variables P(obj.nx,obj.nx) L(obj.nu,obj.nx) Y(obj.nu+obj.nx,obj.nu+obj.nx) hinf_square tao(1,info.T) 
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

                        upper_left = [P-obj.G*obj.G'   zeros(obj.nx,obj.nu+obj.nx);...
                                      zeros(obj.nu+obj.nx,obj.nx)       -Y];
                        upper_right = [obj.G*obj.H';...
                                       Y*[obj.C';...
                                          obj.D']];
                        bottom_right = [hinf_square*eye(obj.ny) - [obj.C obj.D]*Y*[obj.C';...
                                                                                   obj.D'] - obj.H*obj.H'];
                        [upper_left     upper_right;...
                         upper_right'   bottom_right] - sum_over>=0;

                        upper_left_2 = Y-[P L';...
                                          L zeros(obj.nu)];
                        upper_right_2 = [zeros(obj.nx);...
                                         L];
                        [upper_left_2       upper_right_2;...
                         upper_right_2'     P]>=0;
                        tao >= 0;
                        P>=0; Y>=0;
            cvx_end
            hinf = sqrt(hinf_square);
            K = L*inv(P);
            out = struct;
            out.P = P;
            out.L = L;
            out.K = K;
            out.Y = Y;
            out.hinf = hinf;
        end
    end
end


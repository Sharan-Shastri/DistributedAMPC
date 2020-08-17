function cont = prestab_controller(sys,cont)
% function to find a robust control gain matrix
% Computes K such that it is prestabilizing, and that the terminal set is
% RPI and of the form alpha_bar*X0. Any set smaller than this will also
% work as a terminal set. 

warning('Using hardcoded vertices for H_theta. Check their validity')

Ac = NaN*ones(sys.n,sys.n,sys.p);
for k  = 1:size(sys.H_theta_v,2)
    Ac(:,:,k) = sys.A0+ sum(bsxfun(@times,sys.Ap,reshape(sys.H_theta_v(:,k),[1,1,sys.p])),3);     
end 


%% Finding a control gain which is robustly stabilizing
Pinv = sdpvar(sys.n,sys.n);
Y = sdpvar(sys.m,sys.n,'full');
Cons_K = [];
lambda_c = 0.9;
for k = 1:size(sys.H_theta_v,2)
   % contractive
    Cons_K = [Cons_K; 
                [lambda_c*Pinv           Pinv*Ac(:,:,k)'+Y'*sys.B0';
                  Ac(:,:,k)*Pinv+sys.B0*Y     Pinv]>=0];
end
for i = 1:sys.nc
% Constraint satisfaction    
Cons_K = [Cons_K; 
            [1 sys.F(i,:)*Pinv+sys.G(i,:)*Y
             (sys.F(i,:)*Pinv+sys.G(i,:)*Y)' Pinv]>=0 ];
end
        
options = sdpsettings('verbose',0);
diagnostics = optimize(Cons_K,-log(det(Pinv)),options); % 
if diagnostics.problem
   error(diagnostics.info) 
end
%%
% Ensure that the gain is stabilizing for all plants
cont.K = value(Y)*inv(value(Pinv));


% compute polyhedron using ellipsoid data
P = value(Pinv)^-1;
[L,D] = ldl(P);
sDinv = sqrt(inv(D));
d = diag(sDinv);

% polytope 1: outer bound of ellipsoid
% |y| < d
% x < Linv d
H_x = [eye(sys.n); -eye(sys.n)];
h_x = [L\d;L\d];

X_1 = Polyhedron(H_x,h_x);
X_1.minVRep();
xc1 = (Ac(:,:,1)+sys.B0*cont.K)*X_1.V(1,:)';

% max(X_1.A*xc1-X_1.b)
% Not so good performance


% % polytope 2: inscribed polytope in ellipsoid: results in small X0, but
% less vertices
% verts = [L'\sDinv -L'\sDinv];
% X_2 = Polyhedron(verts');
% X_2.minHRep();
% for i = 1:size(verts,2);
%     xc1 = (Ac(:,:,1)+sys.B0*cont.K)*verts(:,i);
%     max(X_2.A*xc1-X_2.b);
% end

cont.H_x = X_1.A;
cont.h_x = X_1.b;
end

%%  LMI with P and K included together

% solve nominal MPC for DGUs

% clc; clear;
% sys = system_desc();

Np = 5;


x_nk = sdpvar(sys.n,Np+1,'full');
u_nk = sdpvar(sys.m,Np,'full');
cost_bounds = sdpvar(Np,1);

cons = [x_nk(:,1) == sys.x0];

for i = 1:Np
    cons = [cons, ...
        x_nk(:,i+1) == sys.A0*x_nk(:,i) + sys.B0*u_nk(:,i),...
        sys.F*x_nk(:,i) + sys.G*u_nk(:,i) <= sys.vec_1_cons, ...
        norm(x_nk(:,i),'inf')+ norm(u_nk(:,i),'inf') <= cost_bounds(i)]; 
end

cons = [cons, x_nk(:,end) == zeros(sys.n,1)];

options = sdpsettings('solver','gurobi','verbose',1,'debug',0);
optimize(cons, sum(cost_bounds));


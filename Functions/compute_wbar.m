function cont =  compute_wbar(sys,cont)
%% Computation of w_bar
nHw = size(sys.H_w,1);
nHx = size(cont.H_x,1);
% max Hx_i w, s.t. H_w w <= h_w
% min h_w' li, s.t. H_w' li == Hx_i, li>=0
lam = sdpvar(nHw,nHx,'full');
Constraints = [lam>=0];
J = 0;
for i = 1:nHx
    Constraints = [Constraints,...
        sys.H_w'*lam(:,i)==cont.H_x(i,:)'];
    J = J+sys.h_w'*lam(:,i);
end
options = sdpsettings('solver','gurobi','verbose',0);
diagnostics = optimize(Constraints,J,options);
if diagnostics.problem
   error(['Calculating w_bar:',diagnostics.info]);
end

for i = 1:nHx
    cont.w_bar(i,1) = value(sys.h_w'*lam(:,i));
end
%% compute other constants
cont.f_bar = max((sys.F+sys.G*cont.K)*cont.x_v,[],2);

cont.alpha_bar = 1/max(cont.f_bar);
% Compute minimum value of alpha
den = zeros(size(sys.H_theta_v,2),1);
for k = 1:size(sys.H_theta_v,2)
        den(k) = max(cont.w_bar./(cont.h_x-cont.H_x*(sys.A0+sum(bsxfun(@times,sys.Ap,reshape(sys.H_theta_v(:,k),[1,1,sys.p])),3)+sys.B0*cont.K)*cont.x_v),[],'all');  
%             min((cont.h_x-cont.w_bar)./abs(cont.H_x*(sys.A0+sum(bsxfun(@times,sys.Ap,reshape(sys.H_theta_v(:,k),[1,1,sys.p])),3)+sys.B0*cont.K)*cont.x_v(:,j))));    
end
cont.alpha_min = max(den);

if cont.alpha_min >cont.alpha_bar
    error('The uncertainty and constriants form an infeasible problem. Relax the problem, or reduce uncertainty');
end
end
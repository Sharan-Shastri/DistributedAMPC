function cont = computeHx(sys,cont)
%% Compute Invariant Set 

%Closed loop given feedback controller K

A_cl = NaN*ones(sys.n,sys.n,sys.p);
for k  = 1:size(sys.H_theta_v,2)
    A_cl(:,:,k) = sys.A0+ sum(bsxfun(@times,sys.Ap,reshape(sys.H_theta_v(:,k),[1,1,sys.p])),3) + sys.B0*cont.K;     
end 

% Starting Polytope :
% F -- State Constraints 
% G -- Input constraints
% A_s=[sys.F(1:8,:);sys.G(9:12,:)*cont.K];
% b_s=sys.vec_1_cons;

A_s= cont.H_x;
b_s=cont.h_x;


x0=zeros(size(b_s));
Aeq = [];
beq = [];
lb = [];
ub = [];
% options=optimoptions(@linprog,'Display','off');
% Desired lambda-contractivity factor (Invariant -> lambda=1)
lambda_desired=0.98;
i=1;
% For every halfspace
while (i<=size(A_s,1))
   %i
   %size(A_s,1)
   a=A_s(i,:);
   b=b_s(i,:);
   % For every vertex of Theta
   for j=1%:sys.p
       % Compute if at the next time step the constraint still holds
       [x,fval]=linprog((-a*A_cl(:,:,j))',A_s,b_s,Aeq,beq,lb,ub);
       c_i(j)=-fval-lambda_desired*b;
   end
   for j=1%:sys.p
       % If it does not -> add new halfspace
       if c_i(j)>0
           A_s=[A_s;1/lambda_desired*a*A_cl(:,:,j)];
           b_s=[b_s;b]; 
       end
   end
   i=i+1;
end

% Use MPT to get minimal amount of halfspaces
X_0=Polyhedron(A_s,ones(size(A_s,1),1));
X_0.minVRep();
cont.H_x=X_0.A;
% Since z=ones(.,1), h_x will also be a vector of ones
cont.h_x=X_0.b;
cont.x_v = X_0.V';
cont.nx_v = length(cont.x_v); % number of vertices
cont.nHx = size(cont.H_x,1);

% Check if X_0 satisfies state and input constraints

end
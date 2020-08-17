function sys = system_desc()

%% Define system parameters

sys.Ts = 1e-4; % Time constant = 1.1e-4s

par.V1 = 100;
par.C1 = 2.2e-3;
par.L1 = 1.8e-3;
par.R1 = 0.2;

par.V2 = 100;
par.C2 = 2.2e-3;
par.L2 = 1.8e-3;
par.R2 = 0.2;

par.Rc_min = 0.045;
par.Rc_true = 0.05;
par.Rc_max = 0.055;

par.Lc = 1.8e-6;

par.V_ref1 = 50;
par.V_ref2 = 50;
par.Il1 = 10;
par.Il2 = 10;

warning('Using true line resistance to calculate reference current')
I_ref1 = (par.V_ref1-par.V_ref2)/par.Rc_true + par.Il1;
I_ref2 = (par.V_ref2-par.V_ref1)/par.Rc_true + par.Il2;
d_ref1 = (par.V_ref1 + par.Il1*par.R1)/par.V1;
d_ref2 = (par.V_ref2 + par.Il2*par.R1)/par.V2;

% Define continuous time system
sys_C.Acont_min = [ -1/(par.Rc_min*par.C1)   1/par.C1 1/(par.Rc_min*par.C1) 0;
    -1/par.L1  -par.R1/par.L1 0 0
    1/(par.Rc_min*par.C2) 0 -1/(par.Rc_min*par.C1) 1/par.C2
    0 0 -1/par.L2 -par.R2/par.L2];

sys_C.Acont_max = [ -1/(par.Rc_max*par.C1)   1/par.C1 1/(par.Rc_max*par.C1) 0;
    -1/par.L1  -par.R1/par.L1 0 0
    1/(par.Rc_max*par.C2) 0 -1/(par.Rc_max*par.C1) 1/par.C2
    0 0 -1/par.L2 -par.R2/par.L2];

sys_C.Acont_true = [ -1/(par.Rc_true*par.C1)   1/par.C1 1/(par.Rc_true*par.C1) 0;
    -1/par.L1  -par.R1/par.L1 0 0
    1/(par.Rc_true*par.C2) 0 -1/(par.Rc_true*par.C1) 1/par.C2
    0 0 -1/par.L2 -par.R2/par.L2]; 

sys_C.Bcont = [0 0 ; par.V1/par.L1 0; 0 0; 0 par.V2/par.L2];

% Scale the matrices: x_bar = Sinv*x
S = diag([50,10,50,10]);
sys_C.Acont_min = inv(S)*sys_C.Acont_min*S;
sys_C.Acont_max = inv(S)*sys_C.Acont_max*S;
sys_C.Acont_true = inv(S)*sys_C.Acont_true*S;
sys_C.Bcont = inv(S)*sys_C.Bcont;

% zero order hold discretization
[A0_min,B0_min] = c2d(sys_C.Acont_min,sys_C.Bcont,sys.Ts);
[A0_max,B0_max] = c2d(sys_C.Acont_max,sys_C.Bcont,sys.Ts);
[A0_true,B0_true] = c2d(sys_C.Acont_true,sys_C.Bcont,sys.Ts);

% parameters to be identified: theta1: (R_12*C_1)^-1, (R_12*C_2)^-1
% use structure to identify parameters: parameterized as A_min + theta*Ap. 
p_var = max((A0_max - A0_min),[],'all');
sys.theta_true = max(A0_true-A0_min,[],'all')/p_var;
sys.theta_true = [sys.theta_true, sys.theta_true];

% define matrices for parameter variation
sys.A0 = A0_min;
sys.Ap(:,:,1) = [p_var 0 -p_var 0
                 zeros(3,4)];
sys.Ap(:,:,2) = [zeros(2,4)
                 -p_var 0 p_var 0
                 zeros(1,4)];

sys.B0 = B0_min;

% define dimensions
sys.n = size(sys.Ap,1);
sys.m = size(sys.B0,2);
sys.p = size(sys.Ap,3);

% define bounds on theta: H_theta*theta <= h_theta
sys.H_theta = [eye(sys.p);
               -eye(sys.p)];
sys.h_theta = [1;
               1;
               0;
               0];

sys.H_theta_v = [1 1; 1 0; 0 1; 0 0]';


% Generate H_theta with a predetermined complexity
% sys = boundedComplexity(sys);
sys.nHtheta = length(sys.h_theta);

% define disturbance bounds: H_w*w <= h_w
sys.w_bound = 0.0001;
warning('The additive disturbance is not discretized along with the model');
sys.H_w = [eye(sys.n);-eye(sys.n)]/sys.w_bound;
sys.h_w = [ones(sys.n,1);ones(sys.n,1)];
sys.nHw = length(sys.h_w);


% define state and input constraints: F*x + G*u <= vec_1_cons
% F*S*x_bar + Gu <= vec_1_cons
sys.F = [ 1/(par.V1-par.V_ref1) 0 0 0
          -1/par.V_ref1 0 0 0
          0 1/I_ref1 0 0 
          0 -1/I_ref1 0 0
          0 0 1/(par.V2-par.V_ref2) 0
          0 0 -1/par.V_ref2 0 
          0 0 0 1/I_ref2
          0 0 0 -1/I_ref2
         zeros(4,4)];     
sys.F = sys.F*S;          

sys.G = [zeros(8,2); 
         1/(1-d_ref1) 0
         -1/d_ref1 0
         0 1/(1-d_ref2)
         0 -1/d_ref2];
sys.nc = size(sys.F,1);
sys.vec_1_cons = ones(sys.nc,1);

sys.x0 = [0,-0.7,0,0]';
end

% 
% par.V1 = 10;
% par.C1 = 2.2e-3;
% par.L1 = 1.8e-3;
% par.R1 = 0.2;
% par.I1 = 0.05;
% V_ref = 5;
% I_load = 0.5 ;
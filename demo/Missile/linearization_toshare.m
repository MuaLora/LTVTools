%% Linearization of Missile Model
% This script constructs matrices for the linearization of the missile 
% model presented in the 2000 ACC paper. The states are x=[alpha;q],the  
% input is delta, and the outputs are y=[q;Fz/m]. These outputs are the 
% variables that the three-loop autopilot feeds back. 
% Transfer functions are computed at each point in the gridded state space.
% As in the Ohlemeyer paper, the tf from delta to q has one zero while the
% tf from delta to acceleration has two.

% physical quantities
rho = 0.0023769; % slugs/ft^3
d = 0.5;
Ad = 0.25*pi*d^2;
m = 250;
Iyy = 49.4;
x0 = [0;0;3050];

V = 3050;

% construct state derivatives to find equilibrium (delta,q) given (alpha,V)
p1 = @(a) -21*a - 24.5*a^3;
c1 = -2*6/pi;
CD = 0.3;
p2 = @(a) -15.75*a - 36.75*a^3;
c2 = -2*42/pi;
Cmtd = -720;

R = @(a,V) cos(a)*rho*V*Ad*p1(a)/(2*m);
S = @(a,V) cos(a)*rho*V*Ad*c1/(2*m);
T = @(a,V) (1/4)*sin(a)*cos(a)*rho*V^2*Ad;
G = @(a,V) sin(a)*cos(a)*CD;
H = @(a,V) rho*V^2*Ad*d*p2(a)/(2*Iyy);
J = @(a,V) rho*V^2*Ad*d*c2/(2*Iyy);
L = @(a,V) rho*V*Ad*d^2*Cmtd/(4*Iyy);

% equilibrium  q
q = @(a,V,deltaP) -R(a,V)-S(a,V)*deltaP-T(a,V)*(a+deltaP)^2-G(a,V);

% coeffients of quadratic function to be solved for equlibrium delta
quad2 = @(a,V) [-L(a,V)*T(a,V) (J(a,V)-L(a,V)*(T(a,V)*2*a+S(a,V)))...
    H(a,V)-L(a,V)*(R(a,V)+G(a,V)+T(a,V)*a^2)];

r = @(a,V) roots(quad2(a,V));

% derivatives with respect to alpha
p1p = @(a) -21 - 73.5*a^2;
p2p = @(a) -15.75 - 110.25*a^2;

Gp = @(a,V) (cos(a)^2-sin(a)^2)*CD;
Rp = @(a,V) (-sin(a)*p1(a)+cos(a)*p1p(a))*rho*V*a/(2*m);
Sp = @(a,V) -sin(a)*rho*V*Ad/(2*m)*c1;
Tp = @(a,V) (cos(a)^2-sin(a)^2)*rho*V^2*Ad/4;
Hp = @(a,V) rho*V^2*Ad*d/(2*Iyy)*p2p(a);

% Jacobian linearization
gradf1 = @(a,V,deltaP) [Gp(a,V) + Rp(a,V) + Sp(a,V)*deltaP ...
    + Tp(a,V)*(a+deltaP)^2 + 2*(a+deltaP)*T(a,V), 1];
gradf2 = @(a,V) [Hp(a,V) L(a,V)];
Alin = @(a,V,deltaP) [gradf1(a,V,deltaP); gradf2(a,V)];

f1delta = @(a,V,deltaP) S(a,V)+2*T(a,V)*(a+deltaP);
f2delta = @(a,V) J(a,V);
Blin = @(a,V,deltaP) [f1delta(a,V,deltaP); f2delta(a,V)];

Clin = @(a,V)[0 1; rho*V^2*Ad*p1p(a)/(2*m) 0];
Dlin = @(V) [0; rho*V^2*Ad*c1/(2*m)];

% grid the state space
na = 31;
nV = 21;
A = zeros(2,2,na,nV);
B = zeros(2,1,na,nV);
C = zeros(2,2,na,nV);
D = zeros(2,1,na,nV);
aVec = linspace(-50*pi/180, 50*pi/180,na);
vVec = linspace(1000,3200,nV);

% populate matrices 
xeq = zeros(2,na,nV);   % Eq state
ueq = zeros(1,na,nV);   % Eq input 
yeq = zeros(2,na,nV);   % Eq output
for i = 1:na
    a = aVec(i);
    for j = 1:nV
        V = vVec(j);
        % to find d_eq we solve a quadratic equation; here we pick the
        % appropriate root
        if a>=0
            d_eq = min(r(a,V));
        else
            d_eq = max(r(a,V));
        end
        q_eq = q(a,V,d_eq); % not used
        A(:,:,i,j) = Alin(a,V,d_eq);
        B(:,:,i,j) = Blin(a,V,d_eq);
        C(:,:,i,j) = Clin(a,V);
        D(:,:,i,j) = Dlin(V);
        
        % Store Equilibrium Info
        xeq(:,i,j) = [a; q_eq];
        ueq(1,i,j) = d_eq;
        
        CNzeq = -21*a - 24.5*a^3 + -2*6/pi*d_eq; % XXX Check this
        Fzeq = 0.5*rho*V^2*Ad*CNzeq;
        yeq(:,i,j) = [q_eq; Fzeq/m];
    end
end

% transfer functions
sys = ss(A,B,C,D);
transfer_functions =tf(sys);
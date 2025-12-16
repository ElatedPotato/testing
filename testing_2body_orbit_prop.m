clear all

close all

global mu 
options     = odeset('RelTol',1e-13,'AbsTol',1e-13);
mu = 398600;

earth_radius = 6378.1;
alt = 500;

r = alt + earth_radius;

tf = 2 * pi * sqrt(r^3/mu);
t0 = 0;

x = [r;0;0;0;sqrt(mu/r);0];
phi0 = eye(6);


[tt,XX] = ode45(@(tt,XX) EOM_STM45(tt,XX), [t0 tf], [x; reshape(phi0,36,1)], options);


X           = XX(end,1:6).';
PHI        = reshape(XX(end,7:end),6,6);



function YY = EOM_STM45(t,XX)
% Integrate Equations of Motion and STM
global mu

x   = XX(1);
y   = XX(2);
z   = XX(3);
vx  = XX(4);
vy  = XX(5);
vz  = XX(6);

%uncomment for STM
PHI = reshape(XX(7:end),6,6);




r       = sqrt(x^2 + y^2 + z^2);




% % Integrate States
Xdot    = [vx vy vz (-mu*x/r^3) (-mu*y/r^3) (-mu*z/r^3)];



%TWO BODY ANALYTICAL PARTIALS
G(1,1)  = 3*mu*x^2/r^5 - mu/r^3 ;
G(1,2)  = 3*mu*x*y/r^5;
G(1,3)  = 3*mu*x*z/r^5 ;
G(2,1)  = 3*mu*x*y/r^5 ;
G(2,2)  = 3*mu*y^2/r^5 - mu/r^3 ;
G(2,3)  = 3*mu*y*z/r^5 ;
G(3,1)  = 3*mu*x*z/r^5 ;
G(3,2)  = 3*mu*y*z/r^5 ;
G(3,3)  = 3*mu*z^2/r^5 - mu/r^3 ;

AA11    = zeros(3,3);
AA12    = eye(3,3);
AA21    = G;
AA22    = zeros(3,3);
A       = [AA11 AA12; AA21 AA22];

PHIdot  = A*PHI;



YY      = [Xdot reshape(PHIdot,1,size(PHIdot,1)*size(PHIdot,2))]';
end

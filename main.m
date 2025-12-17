clear all
close all

startup %add paths kalman filter


global mu Re omega etrot grav_vec
%% Load IOD Data


% load marktestdatajuly128hrs
% load markIODTBwnoise30arc7point5sec2hrs
% load markHDdata71323 %loading the high fidelity data
load testdata17July2013LOS


const_moon
% mu = 1;








% omega = 0;


% cspice_furnsh( 'naif0012.tls.pc' )
% cspice_furnsh( 'de440.bsp' )

% et = cspice_str2et('2021 Aug 12 12:00');

%% Load Initial and Final Conditions
spice_path = 'C:\Users\Josh\Desktop\Orbital\SPICE\';
cspice_furnsh([spice_path,'receding_horiz_3189_1burnApo_DiffCorr_15yr.bsp'])
% cspice_furnsh('naif0012.tls')
cspice_furnsh([spice_path, 'naif0012.tls.pc'])
cspice_furnsh([spice_path, 'de440.bsp']) %use 440
cspice_furnsh([spice_path, 'pck00011.tpc'])



disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Parameters

%     Re      = 6378.1;       % km
%     Re      = 1738.1;
%     Day     = 86400;        % s
%     mu      = 398600.44;    % km^3 / s^2
%     mu      = 4902.799;
%     omega   = 2*pi/Day;     % Earth rotation
ObsErr  = 2;            % arcsec (white noise 1-sigma)

flag =0;                % 0 means that we're using the range and range rate


%     omega       = (2*pi)/(655.720*3600);% Moon Rotation Rate
% %     omega = 0;
%     Re      = 1738.0;
%     mu      = 4902.799056; %from lunar prospector data
%     et = 812545406.534024;
grav_vec = 20;
sv = states;
struc_srp.mass = 4.669874020241132e+02; %kg
struc_srp.C_r = 1.5; % reflectance coefficient
struc_srp.A_s = 10; % assume a sphere with a projected area of 10m^2



%% Change RA and Dec measurements to LOS vectors
%
% for i = 1:length(Dec)
%     L(:,i)  = [cos(Dec(i))*cos(RA(i));cos(Dec(i))*sin(RA(i));sin(Dec(i))];
%
%     L2(:,i)  = [cos(Decsv(i))*cos(RAsv(i));cos(Decsv(i))*sin(RAsv(i));sin(Decsv(i))];
% end
%

%% KM IOD
gm =  4902.799056;

meas = L(:,25:28);
meas = [RA(42:45),Dec(42:45)];
meas = [G(:,313),G(:,(313+20)),G(:,(313+40)),G(:,(313+60))];
dec =0;
%     [r1,v1] = KMIOD(x1(:,42),x1(:,43),x1(:,44),x1(:,45),meas,gm,t(26)-t(25),dec);
[r1,v1] = KMIOD(x1(:,313),x1(:,(313+20)),x1(:,(313+40)),x1(:,(313+60)),meas,gm,t(20)-t(1),dec);%need to space out measurements enough can perform analysis on how long I need to do this for


%     r2 = states(1:3,42) + 10*randn(3,1);
%     v2 = states(4:6,42) + 0.01*randn(3,1);
% 313:471





%     RA = [RA(329:494)];
%     Dec = [Dec(329:494)];
%     measurements = [RA,Dec];
measurements = -G(:,313:471).';
%     sigmaA = [1e-3 1e-3 1e-3].';%from Burt
%     sigmaA = [1e-1 1e-1 1e-1].';%from Burt
sigmaA = 1e-16*[1 1 1].'*TU^2/DU;
tspan = [tspan(313:471)];

obs = x1(:,313:471);
v_obs = vecGS(4:6,313:471);
%     obs = x1(:,:)/DU;
%     v_obs = vecGS(4:6,:)/VU;


%     for i=1:25
%     r2 = states(1:3,313) + 10*randn(3,1);
%     v2 = states(4:6,313) + 0.01*randn(3,1);
% %         r2 = states(1:3,1) + 10*randn(3,1);
% %     v2 = states(4:6,1) + 0.01*randn(3,1);
% %     disp((abs([r1;v1])-abs([r2;v2])))
%     r1=r2;
%     v1=v2;



%% Kalman Filter
% Measurement Error & Covariance
sigmaang = 30/3600*pi/180;
%     tic
R       = diag([sigunitvector sigunitvector sigunitvector]).^2;
P0      = blkdiag(10^1/DU*eye(3),10^-2/VU*eye(3)).^2;
%     P0      = blkdiag(10^1/*eye(3),10^-2*eye(3)).^2;
A = ones(6,6);
B = triu(A,1)*1e-12;
C = tril(A,-1)*1e-12;
P0 = P0+B+C;





%     measurements = [RA(1:end),Dec(1:end)]*pi/180;
%     sigmaA = [1e-15 1e-15 1e-15].';%from Burt
%     tspan = time(1:36);
%
%     obs = x1(:,1:end)/DU;
%     v_obs = vecGS(4:6,1:end)/VU;

% [Xout,P,Del_y,itr] = BatchLSQtest([r1/DU;v1/VU],P0,tspan/TU,measurements,Rsite,Vsite)
%     [X,Pk,Kk,res,strc] = kalmanfilter([r1/DU;v1/VU],tspan/TU,P0,measurements,R,obs,v_obs,sigmaA);
LOS = true;
tic
[X,Pk,Kk,res,strc] = kalmanfilter([r1/DU;v1*TU/DU],tspan/TU,P0,measurements,R,obs/DU,v_obs/DU*TU,sigmaA,1,et,struc_srp,DU,TU,LOS);
toc
%     strc.X(1:3,:) = strc.X(1:3,:)*DU;
%     strc.X(4:6,:) = strc.X(4:6,:)*VU;
XX = strc.X(1:6,:);

%     strc.ypost

yy2 =yy;
% % load radardata
% % R       = diag(sigmarho*10).^2;
% % measurements = rho2.';
% % 
% % LOS = false;
% % P0      = blkdiag(10^1/DU*eye(3),10^-2/VU*eye(3)).^2;
% % tic
% % [X2,Pk2,Kk2,res2,strc2] = kalmanfilter(XX(1:6,end),(tspan2+tspan(end))/TU,Pk,measurements/DU,R,obs/DU,v_obs/DU*TU,sigmaA,1,et,struc_srp,DU,TU,LOS);
% % toc


%     toc

% disp(Pk)
% disp(res)
%
% disp((abs(X(1:6))-abs([r2;v2])))



%         if itr < 600;
%             [elem] = OrbitElem(X(1:3),X(4:6),mu,1e-10);


%                 figure(1);
for i = 1:length(tspan)
    eigs(i,:) = (sqrt(eig(strc.S(:,:,i))));

    Xeigs(i) = eigs(i,1);
    Yeigs(i) =-eigs(i,2);
    Zeigs(i) =-eigs(i,3);

end


hold on
subplot(3,1,1)
%                 plot(res(1,:),'b.','LineWidth',2.0)
hold on
plot(strc.ypost(1,:),'r.-')
plot(3*Xeigs,'k--')
plot(-3*Xeigs,'k--')
%                 ay1 = yline(3*sigunitvector,'-','3-\sigma');
%                 ay2 = yline(-3*sigunitvector,'-','3-\sigma');
legend('Post-Fit Residuals','3\sigma Bounds')
ylim([-10*sigunitvector 10*sigunitvector])
grid on
grid minor
xlabel('Measurements');
ylabel('Residual (Non Dim)');
title('X-Component of LOS Unit Vector')


%                 figure(2);
hold on
subplot(3,1,2)
%                 plot(res(2,:),'b.','LineWidth',2.0)
hold on
plot(strc.ypost(2,:),'r.-')
plot(3*Yeigs,'k--')
plot(-3*Yeigs,'k--')
%                 ay1 = yline(3*sigunitvector,'-','3-\sigma');
%                 ay2 = yline(-3*sigunitvector,'-','3-\sigma');

ylim([-10*sigunitvector 10*sigunitvector])
grid on
grid minor
legend('Post-Fit Residuals','3\sigma Bounds')
xlabel('Measurements');
ylabel('Residual (Non Dim)');
title('Y-Component of LOS Unit Vector')



subplot(3,1,3)
plot(strc.ypost(3,:),'r.-')
hold on
plot(3*Zeigs,'k--')
plot(-3*Zeigs,'k--')
%                 ay1 = yline(3*sigunitvector,'-','3-\sigma');
%                 ay2 = yline(-3*sigunitvector,'-','3-\sigma');
ylim([-10*sigunitvector 10*sigunitvector])
legend('Post-Fit Residuals','3\sigma Bounds')
xlabel('Measurements');
ylabel('Residual (Non Dim)');
title('Z-Component of LOS Unit Vector')
grid on
grid minor
sgtitle('Line of Sight Residuals for 20 min of Lunar Orbiter Data (7.875 sec Measurement Cadence)')
% f = figure;
% movegui(f,'west');
set(gcf,'Position',[100 100 800 800]);
% end
hold off

%propagating the trajectory
%     end
mu = gm;
options     = odeset('RelTol',1e-10,'AbsTol',1e-12);
%     [tt,YY] = ode45(@(tt,YY) EOM_STM_OD_fast_state(tt,YY,1,et,struc_srp), tspan, [strc.X(1:6,1); reshape(eye(6),36,1)], options);

figure()
%%create sphere adjust the number of points to get a smoother sphere
[X1,Y1,Z1] = sphere(500);



X2 = X1 * Re;
Y2 = Y1 * Re;
Z2 = Z1 * Re;
surf(X2,Y2,-Z2)
%% combines the second sphere with Moon topography
moon = imread('moon.jpg');
h = findobj('Type','surface');
set(h,'CData',moon,'FaceColor','texturemap','edgecolor','none')
xlabel('X_{Inertial} (km)')
ylabel('Y_{Inertial} (km)')
zlabel('Z_{Inertial} (km)')
axis equal

hold on
plot3(XX(1,:)*DU,XX(2,:)*DU,XX(3,:)*DU,'k')
%     plot3(YY(:,1),YY(:,2),YY(:,3),'r')
plot3(yy2(313:471,1),yy2(313:471,2),yy2(313:471,3),	'b')
% % plot3(yy(:,1),yy(:,2),yy(:,3),	'b')
% % plot3(strc2.X(1,:)*DU,strc2.X(2,:)*DU,strc2.X(3,:)*DU,'r','LineWidth',1.25)
% load ODbatchdatacomparetoOpticaldata
% plot3(XX(:,1),XX(:,2),XX(:,3),'r')
title('Estimated Orbit Calculated from 20 min Line of sight data Observing Shackleton Crater')

set(gcf,'Position',[900 100 800 800]);




tesvar = yy2(313:471,1:6);
for  i =1:length(XX(1,:))

    positionerr(i) = norm(XX(1:3,i))*DU-norm(tesvar(i,1:3));
    velocityerr(i) = norm(XX(4:6,i))*VU-norm(tesvar(i,4:6));
   

end

% 
% load 
% for  i =1:length(XX(1,:))
% 
%     positionerr(i) = norm(XX(1:3,i))*DU-norm(tesvar(i,1:3));
%     velocityerr(i) = norm(XX(4:6,i))*VU-norm(tesvar(i,4:6));
%    
% 
% end
figure()
subplot(1,2,1)
plot((tspan-tspan(1))/60,positionerr)
xlabel('Time (min)')
ylabel('Position Error (km)')
grid on
grid minor
title('Position')
hold on
subplot(1,2,2)
plot((tspan-tspan(1))/60,velocityerr)
xlabel('Time (min)')
ylabel('Velocity Error (km/s)')
grid on
grid minor
title('Velocity')

sgtitle('Position and Velocity Error when Compared to Truth')




load ODbatchdatacomparetoOpticaldata
tesvar = yy(:,1:6);
for  i =1:length(XX(:,1))

    positionerr2(i) = norm(XX(i,1:3))-norm(tesvar(i,1:3));
    velocityerr2(i) = norm(XX(i,4:6))-norm(tesvar(i,4:6));
   

end
hold on

subplot(1,2,1)
plot((tt-tt(1))/60,positionerr2)
legend('EKF Optical Measurements','Batch LSQ Range & Range-Rate Measurements')
ylim([-3 55])

subplot(1,2,2)
plot((tt-tt(1))/60,velocityerr2)
legend('EKF Optical Measurements','Batch LSQ Range & Range-Rate Measurements')

%% range plots uncomment to make the plots
% % figure()
% % plot(res2*DU,'.')
% % %                 ay1 = yline(3*sigunitvector,'-','3-\sigma');
% % %                 ay2 = yline(-3*sigunitvector,'-','3-\sigma');
% % 
% % 
% % xlabel('Measurements');
% % ylabel('Residual (km)');
% % title('Range Residuals')
% % grid on
% % grid minor


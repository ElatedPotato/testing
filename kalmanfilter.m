function [X,Pk,Kk,res,strc] = kalmanfilter(X,time,P0,Measure,R,observer,Vsite,sigmaA,body,et,struc_srp,LU,TU,LOS)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
II = eye(6);
W = inv(R);
options     = odeset('RelTol',1e-12,'AbsTol',1e-13);
PHI0 = eye(6);
errTol = 1e-8;
IC = X;
err = 1;

strc.PK(:,:,1) = P0;
strc.X(:,1) = X;

if LOS
    first = 1;

else
    first = 2;
end

for i =first:length(time)
    if i ==1
        X = IC;
        PHI = PHI0;
        PHI0 = PHI;

        Pk_bar = PHI*P0*(PHI.') ;%+ blkdiag(10^40*eye(3),10^-40*eye(3));
    else
        dt = time(i) - time(i-1);
        Q = Qdyn(sigmaA.^2,dt,IC);
        % Integration Step
        %     [tt,XX] = ode45(@(tt,XX) EOM_STM_HW(tt,XX), [time(i-1) time(i)], [IC; reshape(PHI0,36,1)], options);
        [tt,XX] = ode45(@(tt,XX) EOM_STM_OD_fast(tt,XX,body,et,struc_srp,LU,TU), [time(i-1) time(i)], [IC; reshape(PHI0,36,1)], options);


        X           = XX(end,1:6).';
        PHI        = reshape(XX(end,7:end),6,6);
        PHI0 = PHI;

        Pk_bar = PHI*P0*(PHI.') + Q;%+blkdiag(1e-10*eye(3),1e-8*eye(3));%;%
    end

    %Measurement Calculations and Measurement Sensitivty matrix calc M is
    %sensitivity matrix and Yes is the measurement calc
    if LOS
        [M,Yest] = MatPopLOS(X,[observer(1:3,i);Vsite(1:3,i)]);

    else
        [M,Yest] = MatPopRadar(XX(1,1:3).',X(1:3));

    end


    yk = Measure(i,:).'-Yest;
    res(:,i)=yk;

    %kalman gain calculation
    var = (M*Pk_bar*(M.')+R);
    %     var(2,2) = var(2,2)*0.999;
    %     var(1,1) = var(1,1)*0.999;
    %     var(3,3) = var(3,3)*0.999;
    Kk  = Pk_bar*(M.')/var;
    X   = X + Kk*(yk);
    Pk = (II-Kk*M)*Pk_bar*(II-Kk*M).'+(Kk*R*Kk.') ;
    %         var=(M*Pk*(M.')+R);

    if LOS
        [~,Ypost] = MatPopLOS(X,[observer(1:3,i);Vsite(1:3,i)]);
        strc.ypost(:,i) = Measure(i,:).'-Ypost;

    else

        [~,Ypost] = MatPopRadar(XX(1,1:3).',X(1:3));
        strc.ypost(:,i) = Measure(i,:).'-Ypost;

    end

    strc.PK(:,:,i) = Pk;
    strc.X(:,i) = X;
    strc.S(:,:,i)=var;
    strc.K(:,:,i)=Kk;
    if i >1
        strc.Q(:,:,i)=Q;
    end
    P0=Pk;%*1.0595;
    IC=X;

end


%end


end

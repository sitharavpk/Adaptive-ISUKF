function  [xa , Pa, K] = ukfmrfnli(xa, Pa, R, y, alpha, beta, kappa,index)

global E; global V; global Q;

xa = [xa, 1]';
Pa = [Pa 0 ; 0 R];
na = size(xa,1);
 
%% Prediction Step

%Compute Sigma Points and Weights
[xSigmaPts, wSigmaPts, nsp] = scaledSymmetricSigmaPoints(xa, Pa, alpha, beta, kappa);
 
 Xa = xSigmaPts(1:na-1,:);
 xn_dep = xSigmaPts(na,:);
 
 gf = Xa;
 Ya = gf .* xn_dep;
 
 Wm = wSigmaPts(1:nsp);
 Wc(1) = wSigmaPts(nsp+1);
 Wc(2:nsp) = wSigmaPts(2:nsp);
 
 
 %Combine Prior intensity estimate and update error covariance
 xmn = Xa*Wm';
 Xmn = xmn*ones(1, nsp);
 dx = (Xa-Xmn)*sqrt(Wc(1));
 
 Ppred = abs((Wm.*dx)*dx') +Q;

 %% Updation Step
 [xSigmaPts1, wSigmaPts1, nsp1] = scaledSymmetricSigmaPoints([xmn;1], [Ppred 0;0 R], alpha, beta, kappa);
 
 Wm1 = wSigmaPts1(1:nsp);
 Ya = xSigmaPts1(1:na-1,:);
 Wc1(1) = wSigmaPts1(nsp+1);
 Wc1(2:nsp) = wSigmaPts1(2:nsp);
 
 %Update the intensity estimate, innovation covariance and cross covariance 
 ymn = Ya*Wm1';
 Ymn = ymn*ones(1, nsp);
 dy = (Ya-Ymn)*sqrt(Wc1(1));
 
 S1 = (Wc1.*dy)*dy'+ V;

 Pxz = dx*(Wc1.*dy)';



 K = Pxz * inv(S1);         %Kalman Gain
 xa = xmn + K*(y - ymn);    %Intensity Estimate
 Pa = abs(Ppred - K*S1*K'); %Error Covariance
 

E(index+14) = abs(y -(imnoise(xa/255,'speckle',R)*255))*10;  % Intensity Residual
c=1;F =0;

for k=index-15+1:index
    F = F + (E(c)^2);
    c=c+1;
end

temp = Ya + (-y+E(index+14))*ones(1,nsp);

V =((F/15) + Wc1*temp');%/1000;
 
Q = K*F*K';%/10000;

 

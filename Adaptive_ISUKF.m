

% =====================================================================================
% Adaptive Importance Sampling Unscented Kalman Filter based SAR Image Super Resolution
                                   
%                                  - Sithara Kanakaraj, Madhu S. Nair and Saidalavi Kalady

% =====================================================================================

% 
% This code is the Matlab implementation of the paper:
% ----------------------------------------------------
% Sithara Kanakaraj, Madhu S. Nair and Saidalavi Kalady, “Adaptive Importance 
% Sampling Unscented Kalman Filter based SAR Image Super Resolution”, Computers and 
% Geosciences, Elsevier, Vol. 133, Article No. 104310, December 2019. 
% DOI: 10.1016/j.cageo.2019.104310
%
% If you find the code useful please consider citing us.
%
% Input:
%     - filePath: the file name of the input containing 16 low-resolution
%     images with translational shift and speckle noise
%     
% Output:
%     - HR: the super-resolved image
%
% Example usage:
%   filePath = 'data\synthetic\synthetic1.mat';
%
%   HR = Adaptive_ISUKF(filePath);
% 
% Disclaimer: The assessment metric values in the paper are the best results for the ideal cases.
%
%   It is provided for educational/research purpose only.
%   

%function HR = Adaptive_ISUKF(filePath)
global E; global V; global Q;

filePath = 'data\synthetic\synthetic2.mat';
% =========================================================================
% Read the input image file
% =========================================================================

LR = cell2mat(struct2cell(load(filePath)));
figure,imshow(LR(:,:,1)),title('Low-Resolution Input');

E= ones(1,35)*10;   % Initialize the intensity residual
Q = 0.01;           % Initialize the Process Noise Covariance
V = 2;              % Initialize the Measurement Noise Covariance
d=2;                % Decimation Factor



% =========================================================================
% Estimating the Noise Level 
%                               Publically available Code by Xinhao Liu, et. al
% =========================================================================
[pvar i1 i2]=NoiseLevel(LR(:,:,1));         
fprintf('Noise Level of LR image: %.4f\n',pvar);

% =========================================================================
% Image Registration to alignment all LR images to the base image 
%                               Publically available Code by Manuel Guizar-Sicairos, et. al
% =========================================================================
REG_LR=[];
REG_LR(:,:,1)=LR(:,:,1);
usfac = 1;
for i=2:16      %Considering the first LR image to be the base image
       [out, Greg] = dftregistration(fft2(REG_LR(:,:,1)),fft2(LR(:,:,i)),usfac);
       REG_LR(:,:,i)=abs(ifft2(Greg));
end



% =========================================================================
% Initializing size and values of the final HR image
% =========================================================================
[M N] = size(LR(:,:,1));
M=d*M; N=d*N;
HR = ones(M+(2*d-d),N+(2*d-d))*median(median(REG_LR(1:5,:,1))).*255;

fprintf('Adaptive ISUKF....\n');
% =========================================================================
% Adaptive Importance Sampling Unscented Kalman Filtering
% =========================================================================
alpha = 0.1; beta = 2; kappa = 1;
gma = 5.8;
Hk = 1; kv = 0.00001;

for m = d:M+(d-1)
    fprintf('%d\t',m);
    for n = d:N+(d-1)
        %% Determining the estimates of the moments using Discontinuity Adaptive Markov random field prior
        px = [HR(m,n-1), HR(m-1,n), HR(m-1,n-1), HR(m-1,n+1)];
        mp = mean(px);
        Qis = kv*(mp^1.5);
        [phi, vr] = dmrfchym4(px, Qis, gma);    %finding mean and variance
        dxa = phi;
        dPa = vr;
        
        for i=1:16
            y=REG_LR(:,:,i);y=y(floor(m/d),floor(n/d)).*255; 
            [dxa,dPa,K] = ukfmrfnli(dxa, dPa, pvar, y, alpha, beta, kappa,i);
            if dxa(1) < 1
                dxa(1) = 1;
            elseif  dxa(1) > 255
                dxa(1) = 255;
            end
        end
        HR(m,n) = dxa(1);
        sigma(m,n)=sqrt(dPa);
    end
end
damrf = HR(d:M+(d-1), d:N+(d-1));

fprintf('\nPost-Processing....\n');

% ==============================================================================
% Post-processing using Discontinuity Adaptive Non-Local Mean (DA-NLM) Filtering
% ==============================================================================
HR = DANLMF(damrf,sigma);

 
 
% =========================================================================
% Output quality assessment metrics
% =========================================================================
O = HR./255;
figure,imshow(O),title('High-Resolution Output');

%       ###############################
%       ###############################
%           Only for Synthetic Images
%       ###############################
%       ###############################

% I = im2double(rgb2gray(imread('data\synthetic\synthetic1.png')));        %For Synthetic Image 1
I = im2double(imread('data\synthetic\synthetic2.png')); O = O(2:256,1:255);I = I(1:255,1:255);%For Synthetic Image 2

ps1=psnr(O,I);
[ps2, smap]=ssim(O,I,[0.01 0.03],fspecial('gaussian', 11, 1.5),1);
ps3=FeatureSIM(O,I);
ps4=epf(O,I);
fprintf('\nPSNR: %.4f\nSSIM: %.4f\nFSIM: %.4f\nEPF: %.4f\n',ps1,ps2,ps3,ps4);




%       ###############################
%       ###############################
%           Only for Real SAR Images
%       ###############################
%       ###############################

% reg = O(240:300,350:415);          % Real SAR Image 1
% reg = O(330:380,210:260);          % Real SAR Image 2
% reg = O(256:282,288:316);          % Real SAR Image 3

% ENL = (mean2(reg)/std2(reg))^2;
% fprintf('\nENL: %.4f\n',ENL);




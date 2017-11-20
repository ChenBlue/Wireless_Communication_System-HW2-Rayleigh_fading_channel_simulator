% Implement a Rayleigh fading channel simulator based on the Filtered 
% Gaussian Noise method, and plot the channel output for fmT = 0.01, 0.1 
% and 0.5. 
fmT = [0.01 0.1 0.5]; % 3 different fm*T
[row, num] = size(fmT);
Omgp = 1; % Set average power as 1
sample_num = 300; % channel output data point
T = 1; % simulation step size

% Filter two independent white Gaussian noise sources with low-pass
% filters. First-order lowpass filter is shown as following:
%
%       [Gi(k+1),Gq(k+1)] = sigma*[Gi(k),Gq(k)] + (1-sigma)*[W1(k),W2(k)]
%
% Gi is in-phase part of output, Gq is quadrature part of output, sigma is 
% the coefficient of first order low-pass filter.
sigma = 2-cos(pi.*fmT./2)-sqrt((2-cos(pi.*fmT./2)).^2-1);
var = (1+sigma)./(1-sigma).*Omgp./2; % variance of Gaussian noise source

w1 = zeros(num, sample_num); % Gaussian noise source 1 for Gi
w2 = zeros(num, sample_num); % Gaussian noise source 2 for Gq
for i = 1:num
    w1(i,:) = normrnd(0,sqrt(var(i)),1,sample_num);
    w2(i,:) = normrnd(0,sqrt(var(i)),1,sample_num);
end

gI = ones(num, sample_num); % In-phase part of output
gQ = ones(num, sample_num); % Quadrature part of output
% Derive the value by low-pass filter
sigma = sigma';
for j = 1:sample_num-1
    gI(:,j+1) = sigma.*gI(:,j)+(1-sigma).*w1(:,j);
    gQ(:,j+1) = sigma.*gQ(:,j)+(1-sigma).*w2(:,j);
end

envelope = sqrt(gI.^2+gQ.^2); % Power of output
envelope_dB = 10*log10(envelope); 

x_axis = (1:sample_num)./T; % time axis: t/T (0~300)
% Plot the channel output with 3 different fm*T
figure,plot(x_axis, envelope_dB(1,:),'g',x_axis, envelope_dB(2,:),'b',x_axis, envelope_dB(3,:),'r')
title('Filtered Gaussian Noise Method');
xlabel('Time, t/T');
ylabel('Envelope Level (dB)');
legend('fmT=0.01','fmT=0.1','fmT=0.5');

% Implement a Rayleigh fading channel simulator based on the Filtered 
% Gaussian Noise method, and plot the channel output autocorrelation
% for fmT = 0.01, 0.1 and 0.5. 
fmT = [0.01 0.1 0.5]; % 3 different fm*T
[row, num] = size(fmT);
Omgp = 1; % Set average power as 1
sample_num = 30000;  % channel output data point
T = 1;  % simulation step size

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
g = sqrt(2)*(gI+1i*gQ);

fm = fmT./ T;
tau = 10./fm; % largest tau, Requirement: fm*tau=0~10
% Autocorrelation for three different fm*T, the size of autocorrelation is
% different.
phi1 = zeros(1, tau(1)+1);
phi2 = zeros(1, tau(2)+1);
phi3 = zeros(1, tau(3)+1);

% fm*T = 0.01
for i=0:tau(1)
    g_shift = zeros(1,sample_num); 
    g_shift(1, i+1:end) = g(1,1:end-i); % gI shift right i index and put it in gI_shift
    phi1(1,i+1) = mean(conj(g(1,:)).*g_shift); % Multiply each other and average them
end

% fm*T = 0.1
for i=0:tau(2)
    g_shift = zeros(1,sample_num);
    g_shift(1, i+1:end) = gI(2,1:end-i);
    phi2(1,i+1) = mean(conj(g(2,:)).*g_shift);
end

% fm*T = 0.5
for i=0:tau(3)
    g_shift = zeros(1,sample_num);
    g_shift(1, i+1:end) = gI(3,1:end-i);
    phi3(1,i+1) = mean(conj(g(3,:)).*g_shift);
end

% Plot the channel autocorrelation with 3 different fm*T
figure,plot(fm(1).*(0:tau(1)), phi1./abs(phi1(1)),'r',fm(2).*(0:tau(2)), phi2./abs(phi2(1)),'g',fm(3).*(0:tau(3)), phi3./abs(phi3(1)),'b');
title('Autocorrelation of Filtered Gaussian Method');
xlabel('f_m\tau');
ylabel('Autocorrelation');
legend('fmT=0.01','fmT=0.1','fmT=1');
grid on

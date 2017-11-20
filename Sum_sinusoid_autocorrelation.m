% Implement a Rayleigh fading channel simulator based on the Sum of
% Sinusoids method, and plot the channel output for fmT = 0.01, 0.1 
% and 0.5, and for M = 8, 16.
fmT = [0.01;0.1;0.5]; 
T = 1; % simulation step size
fm = fmT./T; % fm
num = 3;
Omgp = 1; % Set average power as 1
sample_num = 30000; % channel output data point
M = 8;
m = (1:M);
N = 4*M+2;
n = (1:N);
theta_n = 2*pi*n/N;  % theta_n is uniformly distributed
theta_m = theta_n(1:M);
beta_m = repmat(pi*m/M,num,1);
alpha = 0;
fn = repmat(fm,1,M).*repmat(cos(theta_m),num,1);

cosfn = cos(2*pi*fm(1).*1.*cos(theta_m));
cosfm = cos(2*pi*fm(1)*1);
gI = zeros(num,sample_num+1);
gQ = zeros(num,sample_num+1);

% Use sum of sinusoids to derive gI and gQ
for t = 0:sample_num
    gI(:,t+1) = 2*sum(cos(beta_m).*cos(2*pi*t.*fn),2)+sqrt(2)*cos(alpha).*cos(2*pi.*fm*t);
    gQ(:,t+1) = 2*sum(sin(beta_m).*cos(2*pi*t.*fn),2)+sqrt(2)*sin(alpha).*cos(2*pi.*fm*t);
end
g = sqrt(2)*(gI+i*gQ);

tau = 10./fm; % largest tau. Requirement: fm*tau=0~10
% 3 autocorrelation with 3 different fm
phi1 = zeros(1, tau(1)+1);
phi2 = zeros(1, tau(2)+1);
phi3 = zeros(1, tau(3)+1);

% fm*T = 0.01
for i=0:tau(1)
    gI_shift = zeros(1,sample_num+1);
    gI_shift(1, i+1:end) = gI(1,1:end-i);
    phi1(1,i+1) = mean(gI(1,:).*gI_shift);
end

% fm*T = 0.1
for i=0:tau(2)
    gI_shift = zeros(1,sample_num+1);
    gI_shift(1, i+1:end) = gI(2,1:end-i);
    phi2(1,i+1) = mean(gI(2,:).*gI_shift);
end

% fm*T = 0.5
for i=0:tau(3)
    gI_shift = zeros(1,sample_num+1);
    gI_shift(1, i+1:end) = gI(3,1:end-i);
    phi3(1,i+1) = mean(gI(3,:).*gI_shift);
end

S=besselj(0,2*pi*fm(1).*(0:tau(1))*T); % Ideal autocorrelation

% Plot channel autocorrelation
figure,plot(fm(1).*(0:tau(1)), phi1./abs(phi1(1)),'r',fm(2).*(0:tau(2)), phi2./abs(phi2(1)),'k',fm(3).*(0:tau(3)), phi3./abs(phi3(1)),'b',fm(1).*(0:tau(1)), S,'m--');
title('Autocorrelation of Sum of Sinusoids Method for M=8');
xlabel('f_m\tau');
ylabel('Autocorrelation');
legend('fmT=0.01','fmT=0.1','fmT=1','Ideal');
grid on

% For M = 16; basically the code below is the same as above
M = 16;
m = (1:M);
N = 4*M+2;
n = (1:N);
theta_n = 2*pi*n/N;  % theta_n is uniformly distributed
theta_m = theta_n(1:M);
beta_m = repmat(pi*m/M,num,1);
alpha = 0;
fn = repmat(fm,1,M).*repmat(cos(theta_m),num,1);

cosfn = cos(2*pi*fm(1).*1.*cos(theta_m));
cosfm = cos(2*pi*fm(1)*1);
gI = zeros(num,sample_num+1);
gQ = zeros(num,sample_num+1);

% Use sum of sinusoids to derive gI and gQ
for t = 0:sample_num
    gI(:,t+1) = 2*sum(cos(beta_m).*cos(2*pi*t.*fn),2)+sqrt(2)*cos(alpha).*cos(2*pi.*fm*t);
    gQ(:,t+1) = 2*sum(sin(beta_m).*cos(2*pi*t.*fn),2)+sqrt(2)*sin(alpha).*cos(2*pi.*fm*t);
end
g = sqrt(2)*(gI+i*gQ);

tau = 10./fm; % largest tau. Requirement: fm*tau=0~10
% 3 autocorrelation with 3 different fm
phi1 = zeros(1, tau(1)+1);
phi2 = zeros(1, tau(2)+1);
phi3 = zeros(1, tau(3)+1);

% fm*T = 0.01
for i=0:tau(1)
    gI_shift = zeros(1,sample_num+1);
    gI_shift(1, i+1:end) = gI(1,1:end-i);
    phi1(1,i+1) = mean(gI(1,:).*gI_shift);
end

% fm*T = 0.1
for i=0:tau(2)
    gI_shift = zeros(1,sample_num+1);
    gI_shift(1, i+1:end) = gI(2,1:end-i);
    phi2(1,i+1) = mean(gI(2,:).*gI_shift);
end

% fm*T = 0.5
for i=0:tau(3)
    gI_shift = zeros(1,sample_num+1);
    gI_shift(1, i+1:end) = gI(3,1:end-i);
    phi3(1,i+1) = mean(gI(3,:).*gI_shift);
end

S=besselj(0,2*pi*fm(1).*(0:tau(1))*T); % Ideal autocorrelation

% Plot channel autocorrelation
figure,plot(fm(1).*(0:tau(1)), phi1./abs(phi1(1)),'r',fm(2).*(0:tau(2)), phi2./abs(phi2(1)),'k',fm(3).*(0:tau(3)), phi3./abs(phi3(1)),'b',fm(1).*(0:tau(1)), S,'m--');
title('Autocorrelation of Sum of Sinusoids Method for M=16');
xlabel('f_m\tau');
ylabel('Autocorrelation');
legend('fmT=0.01','fmT=0.1','fmT=1','Ideal');
grid on

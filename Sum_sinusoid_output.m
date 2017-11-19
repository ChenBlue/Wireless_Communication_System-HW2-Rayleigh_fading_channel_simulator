fmT = [0.01;0.1;0.5];
T = 1;
fm = fmT./T;
num = 3;
Omgp = 1;
sample_num = 300;
M = 8;
m = (1:M);
N = 4*M+2;
n = (1:N);
theta_n = 2*pi*n/N;
theta_m = theta_n(1:M);
beta_m = repmat(pi*m/M,num,1);
alpha = 0;
fn = repmat(fm,1,M).*repmat(cos(theta_m),num,1);

cosfn = cos(2*pi*fm(1).*1.*cos(theta_m));
cosfm = cos(2*pi*fm(1)*1);
gI = zeros(num,sample_num+1);
gQ = zeros(num,sample_num+1);
g = zeros(num,sample_num+1);

for t = 0:sample_num
    gI(:,t+1) = 2*sum(cos(beta_m).*cos(2*pi*t.*fn),2)+sqrt(2)*cos(alpha).*cos(2*pi.*fm*t);
    gQ(:,t+1) = 2*sum(sin(beta_m).*cos(2*pi*t.*fn),2)+sqrt(2)*sin(alpha).*cos(2*pi.*fm*t);
end
g = sqrt(2)*(gI+i*gQ);
envelope_dB = 10*log(abs(g)/mean2(abs(g)))
x_axis = (0:sample_num)
figure,plot(x_axis,envelope_dB(1,:),'k',x_axis,envelope_dB(2,:),'b',x_axis,envelope_dB(3,:),'r' )
title('Sum of Sinusoids Method for M=8');
xlabel('Time, t/T');
ylabel('Envelope Level (dB)');
legend('fmT=0.01','fmT=0.1','fmT=0.5');
grid on

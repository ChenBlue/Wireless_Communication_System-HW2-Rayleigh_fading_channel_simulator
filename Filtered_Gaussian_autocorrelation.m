fmT = [0.01 0.1 0.5];
[row, num] = size(fmT);
Omgp = 1;
sample_num = 30000;
T = 1;

sigma = 2-cos(pi.*fmT./2)-sqrt((2-cos(pi.*fmT./2)).^2-1);
var = (1+sigma)./(1-sigma).*Omgp./2;

w1 = zeros(num, sample_num);
w2 = zeros(num, sample_num);
for i = 1:num
    w1(i,:) = normrnd(0,sqrt(var(i)),1,sample_num);
    w2(i,:) = normrnd(0,sqrt(var(i)),1,sample_num);
end

gI = ones(num, sample_num);
gQ = ones(num, sample_num);

for i = 1:num
    for j = 1:sample_num-1
        gI(i,j+1) = sigma(i)*gI(i,j)+(1-sigma(i))*w1(i,j);
        gQ(i,j+1) = sigma(i)*gQ(i,j)+(1-sigma(i))*w2(i,j);
    end
end

fm = fmT./ T;
tau = 10./fm;
phi1 = zeros(1, tau(1)+1);
phi2 = zeros(1, tau(2)+1);
phi3 = zeros(1, tau(3)+1);

gI_shift = zeros(1,sample_num);

for i=0:tau(1)
    gI_shift = zeros(1,sample_num);
    gI_shift(1, i+1:end) = gI(1,1:end-i);
    phi1(1,i+1) = mean(gI(1,:).*gI_shift);
end

for i=0:tau(2)
    gI_shift = zeros(1,sample_num);
    gI_shift(1, i+1:end) = gI(2,1:end-i);
    phi2(1,i+1) = mean(gI(2,:).*gI_shift);
end

for i=0:tau(3)
    gI_shift = zeros(1,sample_num);
    gI_shift(1, i+1:end) = gI(3,1:end-i);
    phi3(1,i+1) = mean(gI(3,:).*gI_shift);
end

figure,plot(fm(1).*(0:tau(1)), phi1./abs(phi1(1)),'r',fm(2).*(0:tau(2)), phi2./abs(phi2(1)),'g',fm(3).*(0:tau(3)), phi3./abs(phi3(1)),'b');
title('Autocorrelation');
xlabel('f_m\tau');
ylabel('Autocorrelation');
legend('fmT=0.01','fmT=0.1','fmT=1');
grid on

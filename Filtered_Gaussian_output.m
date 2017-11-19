fmT = [0.01 0.1 0.5];
[row, num] = size(fmT);
Omgp = 1;
sample_num = 300;
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

envelope = sqrt(gI.^2+gQ.^2);
envelope_dB = 10*log10(envelope);

x_axis = (1:sample_num)./T;
figure,plot(x_axis, envelope_dB(1,:),'g',x_axis, envelope_dB(2,:),'b',x_axis, envelope_dB(3,:),'r')
title('Filtered Gaussian Noise method');
xlabel('Time, t/T');
ylabel('Envelope Level (dB)');
legend('fmT=0.01','fmT=0.1','fmT=0.5');


function quiz3()
close all

% Q 1
x = 0:.01:10;
thresholds = [5.667 5.978 5.830 2.69];
delta = 0:.01:1;
s1 = normpdf(x, 5, 0.5);
s2 = normpdf(x, 7, 1);
figure;
hold on;
plot(x, s1, 'r-', x, s2, 'b-');
for threshold = thresholds
  plot(threshold, delta, 'k-');
end % for


% Q 5
p = 0.1;
HF = -p * log2(p) - (1 - p) * log2(1 - p);


% Q 6
ps = 0.1;
qs = 1/2;    % P(firing | stimulus)
qns = 1/18;  % P(firing | no stimulus)

MI = HF - ps * (-qs * log2(qs) - (1 - qs) * log2(1 - qs)) - (1 - ps) * (-qns * log2(qns) - (1 - qns) * log2(1 - qns));


% Q 11
load tuning
figure;
hold on;
plot(stim, mean(neuron1), 'r-');
plot(stim, mean(neuron2), 'g-');
plot(stim, mean(neuron3), 'b-');
plot(stim, mean(neuron4), 'm-');
%legend('\mu_1', '\sigma_1', '\mu_2', '\sigma_2', '\mu_3', '\sigma_3', '\mu_4', '\sigma_4');
legend('\mu_1', '\mu_2', '\mu_3', '\mu_4');


% Q 12
figure;
hold on;
T = 10; % s
plot(mean(neuron1)*T, var(neuron1)*T^2, 'r.');
plot(mean(neuron2)*T, var(neuron2)*T^2, 'g.');
plot(mean(neuron3)*T, var(neuron3)*T^2, 'b.');
plot(mean(neuron4)*T, var(neuron4)*T^2, 'm.');
legend('F_1', 'F_2', 'F_3', 'F_4');
axis equal


% Q 13 (Note 0 is +y and 90 is +x)
rmean = mean([r1' r2' r3' r4']);
rmax = max([mean(neuron1)' mean(neuron2)' mean(neuron3)' mean(neuron4)']);
c = [c1; c2; c3; c4];
vpop = rmean ./ rmax * c;
degrees = mod(90 - atand(vpop(2)/vpop(1)), 360);
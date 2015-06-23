function quiz5

close all

% problem 1
pkg load image
img = imread('q5_pave.jpg');
rimg = imread('q5_paveBlurred.jpg');
f1 = [1 0 0 0 0; 1 0 0 0 0; 1 0 0 0 0; 0 0 0 1 0; 0 0 0 1 0];
f2 = [1 0 0 0 -1; -1 1 0 0 0; 0 -1 1 0 0; 0 0 -1 1 0; 0 0 0 -1 1]*3/4;
f3 = eye(5)/2;
f4 = [1 1 0 0 1; 1 1 1 0 0; 0 1 1 1 0 ; 0 0 1 1 1; 1 0 0 1 1]/3; <== blurring

% problem 3
W = 0.1 * ones(5) + 0.5*eye(5);
u = [0.6 .5 .6 .2 .1]';
M = 1/4 * ([0 0 1 1 0; 0 0 0 1 1; 1 0 0 0 1; 1 1 0 0 0; 0 1 1 0 0] - eye(5));
I = eye(5);

vss = -(M - I) \ W*u

h = W*u;
[e lambdas] = eig(M);
lambdas = diag(lambdas);
vss = zeros(5, 1);
for i = 1:5
  vss = vss + dot(e(:, i), W*u) * e(:, i) / (1 - lambdas(i));
end % for
vss

% problem 4
Q = [0.15 0.1; 0.1 0.12];
[e lambdas] = eig(Q);
w = 2*e(:, 2);


% problem 8
% t_peaks = 0.5:0.5:10; % ms
% no_peaks = length(t_peaks);
% output_spike_counts = zeros(1, no_peaks);
% for i = 1:no_peaks
  % output_spike_counts(i) = alpha_neuron(t_peaks(i));
% end % for t_peak
% figure;
% plot(t_peaks, output_spike_counts, 'r-')

% problem 13
% eta = 1;
% alpha = 1;
% dt = 0.01;

% data = load('c10p1.mat');
% X = bsxfun(@minus, data, mean(data)); % center about (0,0)
% N = size(X, 1);
% R = X'*X/N;
% [e l] = eig(R);

% figure;
% plot(X(:, 1), X(:, 2), 'b.');
% axis equal
% hold on;
% plot(e(1, 1), e(2, 1), 'mx', e(1, 2), e(2, 2), 'mx');
% drawnow;
% for i = 1:5
  % beep
  % fprintf(1, 'Trial %d\n', i);
  % fflush(1);
  % w = 2*rand(size(X, 2), 1) - 1;
  % for j = 1:100000
    % for k = 1:N
      % u = X(k, :)';
      % v = dot(u, w);
      % w = w + dt * eta * (v * u - alpha * v * v * w);
    % end % for k
    % if mod(j, 1000) == 0
      % fprintf(1, '.');
      % if mod(j, 10000) == 0
        % fprintf(1, '%d\n', j);
      % end % if
     % fflush(1);
    % end % if
  % end % for j
  % fprintf(1, '\n');
  % plot(w(1), w(2), 'r+');
  % drawnow;
% end % for i

% problem 14
data = load('c10p1.mat');
eta = 1;
alpha = 1;
dt = 0.01;

for h = 1:3
  figure;
  axis equal
  hold on;
  beep
  x = 2*rand - 1;
  y = 2*rand - 1;
  X = bsxfun(@minus, data, [x y] + mean(data)); % center about (x, y)
  N = size(X, 1);
  R = X'*X/N;
  [e l] = eig(R);
  mu = mean(X);
  plot(X(:, 1), X(:, 2), 'b.', mu(1), mu(2), 'k+');
  plot(e(1, 1), e(2, 1), 'mx', e(1, 2), e(2, 2), 'mx');
  drawnow;
  for i = 1:5
    fprintf(1, 'Trial %d,%d\n', h, i);
    fflush(1);
    w = 2*rand(size(X, 2), 1) - 1 - [x y]';
    ws = w';
    for j = 1:100000
      for k = 1:N
        u = X(k, :)';
        v = dot(u, w);
        w = w + dt * eta * (v * u - alpha * v * v * w);
      end % for k
      ws = [ws; w'];
      if mod(j, 1000) == 0
        fprintf(1, '.');
        if mod(j, 10000) == 0
          fprintf(1, '%d\n', j);
        end % if
       fflush(1);
       end % if
    end % for j
    fprintf(1, '\n');
    plot(ws(:, 1), ws(:, 2), 'r-');
    plot(w(1), w(2), 'r+');
    drawnow;
  end % for i
end % for h
figure;
hold on;
X = bsxfun(@minus, data, mean(data)); % center about (0,0)
N = size(X, 1);
R = X'*X/N;
[e l] = eig(R);
plot(X(:, 1), X(:, 2), 'b.');
plot(e(1, 1), e(2, 1), 'mx', e(1, 2), e(2, 2), 'mx');
drawnow;

% problem 15
% data = load('c10p1.mat');
% eta = 1;
% alpha = 1;
% dt = 0.01;

% X = bsxfun(@minus, data, mean(data)); % center about (0,0)
% N = size(X, 1);
% R = X'*X/N;
% [e l] = eig(R);

% figure;
% plot(X(:, 1), X(:, 2), 'b.');
% a=2; axis([-a a -a a])
% hold on;
% plot(e(1, 1), e(2, 1), 'mx', e(1, 2), e(2, 2), 'mx');
% drawnow;
% for i = 1:5
  % beep
  % fprintf(1, 'Trial %d\n', i);
  % fflush(1);
  % w = 2*rand(size(X, 2), 1) - 1;
  % ws = w';
  % for j = 1:100000
    % for k = 1:N
      % u = X(k, :)';
      % v = dot(u, w);
      % w = w + dt * eta * v * u;
    % end % for k
    % ws = [ws; w'];
    % if mod(j, 1000) == 0
      % fprintf(1, '.');
      % if mod(j, 10000) == 0
        % fprintf(1, '%d\n', j);
      % end % if
      % fflush(1);
    % end % if
  % end % for j
  % fprintf(1, '\n');
  % plot(ws(:, 1), ws(:, 2), 'r-');
  % drawnow;
% end % for i

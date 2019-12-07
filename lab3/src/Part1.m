clear all, close all

%% Part 1.4 Delay-and-sum beamforming for uniform linear arrays
theta_s = pi/2;           % source angle is 90 degrees
f=2000;                   % frequancy is 2kHz
w=2*pi*f;
c=340;                    % speed of sound in air is ~ 340m/s
l=c/f;                    % wave length
theta=0:0.001:pi;         % theta is a parameter [0,180] degrees

%% Part 1.4.1 
d=0.04;

% We calculate the delay-and-sum beam pattern  in respect to the angle theta
N=4;
B4 = zeros(1,length(theta));
for i=1:length(theta);
    B4(i)=(sin((cos(theta(i))-cos(theta_s))*N*d*w/(2*c))/sin((cos(theta(i))-cos(theta_s))*d*w/(2*c)))/N;
end

N=8;
B8 = zeros(1,length(theta));
for i=1:length(theta);
    B8(i)=(sin((cos(theta(i))-cos(theta_s))*N*d*w/(2*c))/sin((cos(theta(i))-cos(theta_s))*d*w/(2*c)))/N;
end

N=16;
B16 = zeros(1,length(theta));
for i=1:length(theta);
    B16(i)=(sin((cos(theta(i))-cos(theta_s))*N*d*w/(2*c))/sin((cos(theta(i))-cos(theta_s))*d*w/(2*c)))/N;
end

% We plot the absolute values of the results in logarithmic scale using the
% function semilogy
figure;
semilogy(theta,abs(B4),theta,abs(B8),theta,abs(B16));
legend('N=4','N=8','N=16');
title('Delay and Sum beam pattern for N=4, N=8, N=16 (d=0.04)');
xlabel('è (rad)');
ylabel('(dB)');

figure;
semilogy(theta,abs(B4));
title('Delay and Sum beam pattern for N=4 (d=0.04)');
xlabel('è (rad)');
ylabel('(dB)');

figure;
semilogy(theta,abs(B8));
title('Delay and Sum beam pattern for N=8 (d=0.04)');
xlabel('è (rad)');
ylabel('(dB)');

figure;
semilogy(theta,abs(B16));
title('Delay and Sum beam pattern for N=16 (d=0.04)');
xlabel('è (rad)');
ylabel('(dB)');

%% Part 1.4.2 
N=8;

d=0.04;
for i=1:length(theta);
    B4(i)=(sin((cos(theta(i))-cos(theta_s))*N*d*w/(2*c))/sin((cos(theta(i))-cos(theta_s))*d*w/(2*c)))/N;
end

d=0.08;
for i=1:length(theta);
    B8(i)=(sin((cos(theta(i))-cos(theta_s))*N*d*w/(2*c))/sin((cos(theta(i))-cos(theta_s))*d*w/(2*c)))/N;
end

d=0.16;
for i=1:length(theta);
    B16(i)=(sin((cos(theta(i))-cos(theta_s))*N*d*w/(2*c))/sin((cos(theta(i))-cos(theta_s))*d*w/(2*c)))/N;
end

figure;
semilogy(theta,abs(B4),theta,abs(B8),theta,abs(B16));
legend('d=4cm','d=8cm','d=16cm');
title('Delay and Sum beam pattern for d=4cm, d=8cm, d=16cm (N=8)');
xlabel('è (rad)');
ylabel('(dB)');

figure;
semilogy(theta,abs(B4));
title('Delay and Sum beam pattern for d=4cm (N=8)');
xlabel('è (rad)');
ylabel('(dB)');

figure;
semilogy(theta,abs(B8));
title('Delay and Sum beam pattern for d=8cm (N=8)');
xlabel('è (rad)');
ylabel('(dB)');

figure;
semilogy(theta,abs(B16));
title('Delay and Sum beam pattern for d=16cm (N=8)');
xlabel('è (rad)');
ylabel('(dB)');

%% Part 1.4.3
N=8;
d=0.04;
theta_s=0;
theta=-pi:0.001:pi;         %theta is a parameter [-180,180] degrees

% We calculate the absolute value of the delay-and-sum beam pattern 
% in respect to the angle theta (theta_s = 0 degrees)
B = zeros(1,length(theta));
for i=1:length(theta);
    B(i)=(sin((cos(theta(i))-cos(theta_s))*N*d*w/(2*c))/sin((cos(theta(i))-cos(theta_s))*d*w/(2*c)))/N;
end
B1=abs(B);

% Plot the results in logarithmic scale (polar diagram)
figure;
semilogr_polar(theta,B1);
title('Delay and Sum beam pattern for N=8, d=4cm, ès=0\circ');

% We calculate the absolute value of the delay-and-sum beam pattern 
% in respect to the angle theta (theta_s = 45 degrees)
theta_s=pi/4;
theta=-pi:0.001:pi;
for i=1:length(theta);
    B(i)=(sin((cos(theta(i))-cos(theta_s))*N*d*w/(2*c))/sin((cos(theta(i))-cos(theta_s))*d*w/(2*c)))/N;
end
B2=abs(B);

figure;
semilogr_polar(theta,B2);
title('Delay and Sum beam pattern for N=8, d=4cm, ès=45\circ');

% We calculate the absolute value of the delay-and-sum beam pattern 
% in respect to the angle theta (theta_s = 90 degrees)
theta_s=pi/2;
theta=-pi:0.001:pi;
for i=1:length(theta);
    B(i)=(sin((cos(theta(i))-cos(theta_s))*N*d*w/(2*c))/sin((cos(theta(i))-cos(theta_s))*d*w/(2*c)))/N;
end
B3=abs(B);
figure;
semilogr_polar(theta,B3);
title('Delay and Sum beam pattern for N=8, d=4cm, ès=90\circ');

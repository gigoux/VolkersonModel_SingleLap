close all;
clear;
clc;

%Data in order of crosshead displacement [in], loading [lbf], time [sec]
specA_data = importdata('Lab1A.dat');
specB_data = importdata('Lab1B.dat');

specA_data(:, 1) = 0.0254*specA_data(:, 1); %converting in to m
specB_data(:, 1) = 0.0254*specB_data(:, 1); %converting in to m
specA_data(:, 2) = 4.448*specA_data(:, 2); %converting lbf to N
specB_data(:, 2) = 4.448*specB_data(:, 2); %converting lbf to N

specA_loadfail = max(specA_data(:, 2)); %failure loading for specimen A [N]
specB_loadfail = max(specB_data(:, 2)); %failure laoding for specimen B [N]

figure(1)
hold on;
grid on;
plot(specA_data(:, 3), specA_data(:, 2));
plot(specB_data(:, 3), specB_data(:, 2));
xlabel('Time (s)');
ylabel('Loading (N)');
title('Specimen Loading vs Time');
legend('Specimen A', 'Specimen B', 'Location', 'northwest');
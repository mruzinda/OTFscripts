clearvars;
close all;

tstamps = {'2017_05_20_01:46:47',...
           '2017_05_21_03:10:22',...
           '2017_05_21_03:12:45',...
           '2017_05_21_03:24:53',...
           '2017_05_22_01:06:18',...
           '2017_05_22_01:13:27'};

for idx = 1:length(tstamps)
    tstamp = tstamps{idx};
    filename = sprintf('%s_tsys.mat', tstamp);
    A = load(filename);
    if idx == 1
        Tsys = A.Tsys;
    else
        Tsys = Tsys + A.Tsys;
    end
end
Tsys = Tsys/length(tstamps);
figure();
plot(A.freqs, Tsys.');
xlabel('Frequency (MHz)');
ylabel('T_s_y_s (K)');
title('Average Tsys in OTF - Poor Quantization');
ylim([0,100]);
xlim([freqs(1), freqs(end)]);
grid on;
grid minor;
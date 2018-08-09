% Simple script that extracts the covariance matrices for the entire band
% and computes Y-factors + Tsys

close all;
clearvars;

tic;
% System parameters
Dir = '/lustre/gbtdata/TGBT16A_508_01/TMP/BF';

% May 19th, 2017 - 01
% Quant gain = 10
% hot_tstamp = '2017_05_20_01:39:52';
% cold_tstamp = '2017_05_20_01:46:47';
% Thot = 290;
% Tcold = 7.5;

% May 20th, 2017 - 01
% Quant gain = 10
% hot_tstamp = '2017_05_21_03:06:05';
% cold_tstamp = '2017_05_21_03:10:22';
% Thot = 290;
% Tcold = 7.5;

% May 20th, 2017 - 02
% Quant gain = 10
% hot_tstamp = '2017_05_21_03:16:27';
% cold_tstamp = '2017_05_21_03:12:45';
% Thot = 290;
% Tcold = 7.5;

% May 20th, 2017 - 03
% Quant gain = 10
% hot_tstamp = '2017_05_21_03:19:47';
% cold_tstamp = '2017_05_21_03:24:53';
% Thot = 290;
% Tcold = 7.5;

% May 20th, 2017 - 04
% Quant gain = 20
% hot_tstamp = '2017_05_21_03:32:31';
% cold_tstamp = '2017_05_21_03:28:16';
% Thot = 290;
% Tcold = 7.5;

% May 21st, 2017 - 01
% Quant gain = 10
% hot_tstamp = '2017_05_22_01:02:42';
% cold_tstamp = '2017_05_22_01:06:18';
% Thot = 290;
% Tcold = 7.5;

% May 21st, 2017 - 02
% Quant gain = 10
% hot_tstamp = '2017_05_22_01:09:49';
% cold_tstamp = '2017_05_22_01:13:27';
% Thot = 290;
% Tcold = 7.5;

% May 21st, 2017 - 03
% Quant gain = 20
hot_tstamp = '2017_05_22_01:40:27';
cold_tstamp = '2017_05_22_01:36:57';
Thot = 290;
Tcold = 7.5;

% May 21st, 2017 - 03
% Quant gain = 40
% hot_tstamp = '2017_05_22_01:43:10';
% cold_tstamp = '2017_05_22_01:46:57';
% Thot = 290;
% Tcold = 7.5;

% May 22nd, 2017 - 03
% Quant gain = 10
% hot_tstamp = '2017_05_23_02:57:22';
% cold_tstamp = '2017_05_23_03:05:01';
% Thot = 290;
% Tcold = 7.5;

% % May 22nd, 2017 - 04  %%%%%% Not accurate %%%%%%
% % Quant gain = 10
% hot_tstamp = '2017_05_24_20:19:16';
% cold_tstamp = '2017_05_23_03:55:51';
% Thot = 290;
% Tcold = 7.5;

% May 24th, 2017 - 05  %%%%%% GBT test %%%%%%
% Quant gain = 10
% on_tstamp = {'2017_05_25_01:31:35', '2017_05_25_01:32:19', '2017_05_25_01:33:03', '2017_05_25_01:33:47', '2017_05_25_01:34:31', '2017_05_25_01:35:15', '2017_05_25_01:37:26'}; % {'2017_05_26_03:41:28'};
% off_tstamp = {'2017_05_25_01:30:31', '2017_05_25_01:36:22'}; % {'2017_05_26_03:42:05'};
% % beam_num = 1;
% Thot = 290;
% Tcold = 7.5;

LO_freq = 1450;
freqs = ((-249:250)*303.24e-3)+LO_freq;

banks = {'A', 'B', 'C', 'D',...
         'E', 'F', 'G', 'H',...
         'I', 'J', 'K', 'L',...
         'M', 'N', 'O', 'P',...
         'Q', 'R', 'S', 'T'};

bad_freqs = [46:50, 81:100, 146:150, 181:200, 246:250, 281:300, 346:350, 381:400, 446:450, 481:500];
good_idx = [9:19,20,21,23:28,30:33,35:39];
Y = zeros(length(good_idx), 500);

% Get Y factors for each bank
for idx = 1:length(banks)
    bank = banks{idx};
    fprintf('%s', bank);
    cold_filename = sprintf('%s/%s%s.fits', Dir, cold_tstamp, bank);
    hot_filename  = sprintf('%s/%s%s.fits', Dir, hot_tstamp, bank);

    % Get COLD data
    if exist(cold_filename)
        [Rcold, ~, xid, ~] = extract_covariances(cold_filename);
    else
        continue;
    end
    
    % Get HOT data
    if exist(hot_filename)
        [Rhot,  ~, xid, ~] = extract_covariances(hot_filename);
    else
        continue;
    end

    chans = [1:5, 101:105, 201:205, 301:305, 401:405] + 5*xid;
    for i = 1:length(good_idx)
        Phot  = squeeze(mean(Rhot (good_idx(i),good_idx(i),:,:),4)).';
        Pcold = squeeze(mean(Rcold(good_idx(i),good_idx(i),:,:),4)).';
        Y(i,chans) = Phot./Pcold;
    end
end
fprintf('\n');
Tsys = real((Thot - Y*Tcold)./(Y - 1));
Tsys(:,bad_freqs) = NaN;

% Plot results
figure();
plot(freqs, Tsys.');
xlabel('Frequency (MHz)');
ylabel('T_s_y_s (K)');
ylim([0,100]);
xlim([freqs(1), freqs(end)]);
grid on;
grid minor;

save(sprintf('%s_tsys.mat', cold_tstamp), 'freqs', 'Tsys');
% Aggregate weights and plot sensitvity map
close all;
clearvars;

tic;
% System parameters
data_dir = '/lustre/projects/flag';
meta_dir = '/home/gbtdata';


% May 25, 2017 Calibration Grid
% on_tstamp = {'2017_05_26_03:44:16', '2017_05_26_03:45:02', '2017_05_26_03:45:50', ...
%              '2017_05_26_03:47:49', '2017_05_26_03:48:35', '2017_05_26_03:49:21', ...
%               '2017_05_26_03:51:35', '2017_05_26_03:52:21', '2017_05_26_03:53:07', ...
%               '2017_05_26_03:55:06', '2017_05_26_03:55:52', '2017_05_26_03:56:37', ...
%                '2017_05_26_03:58:50', '2017_05_26_03:59:35', '2017_05_26_04:00:21', ...
%                '2017_05_26_04:02:19', '2017_05_26_04:03:05', '2017_05_26_04:03:50', ...
%                '2017_05_26_04:06:01', '2017_05_26_04:06:46', '2017_05_26_04:07:32', ...
%                '2017_05_26_04:09:28', '2017_05_26_04:10:13', '2017_05_26_04:10:59', ...
%                 '2017_05_26_04:13:08', '2017_05_26_04:13:54', '2017_05_26_04:14:39', ...
%                 '2017_05_26_04:16:34', '2017_05_26_04:17:19', '2017_05_26_04:18:05', ...
%                 '2017_05_26_04:20:14', '2017_05_26_04:21:00', '2017_05_26_04:21:45'};
% off_tstamp = {'2017_05_26_03:47:03', '2017_05_26_03:50:40', '2017_05_26_03:54:20', ...
%               '2017_05_26_03:57:55', '2017_05_26_04:01:33', '2017_05_26_04:05:07', ...
%               '2017_05_26_04:08:43', '2017_05_26_04:12:15', '2017_05_26_04:15:49', ...
%               '2017_05_26_04:19:21', '2017_05_26_04:22:55'};
% 
% bad_freqs = [81:100,181:200,281:300,381:400,481:500];
% proj_id = 'AGBT16B_400_02';
% save_dir = sprintf('%s/%s/BF', data_dir, proj_id);
% out_dir = sprintf('%s/mat', save_dir);

% May 24, 2017 Calibration Grid
%on_tstamp = {'2017_05_24_23:50:18','2017_05_24_23:51:00', ...
%             '2017_05_24_23:53:50', '2017_05_25_00:47:53', '2017_05_25_00:48:38', ...
%             '2017_05_25_00:50:52', '2017_05_25_00:51:37', '2017_05_25_00:53:51', ...
%             '2017_05_25_00:54:36', '2017_05_25_00:56:50', '2017_05_25_00:57:35', ...
%             '2017_05_25_00:59:49', '2017_05_25_01:00:34', '2017_05_25_01:02:50', ...
%             '2017_05_25_01:03:33', '2017_05_25_01:31:35', '2017_05_25_01:32:19', ...
%             '2017_05_25_01:33:03', '2017_05_25_01:33:47', '2017_05_25_01:34:31', ...
%             '2017_05_25_01:35:15', '2017_05_25_01:37:26', '2017_05_25_01:38:10', ...
%             '2017_05_25_01:43:04', '2017_05_25_01:43:48', '2017_05_25_01:44:32', ...
%             '2017_05_25_01:45:16', '2017_05_25_01:46:00', '2017_05_25_01:48:12', ...
%             '2017_05_25_01:48:58', '2017_05_25_01:49:40', '2017_05_25_01:50:24', ...
%             '2017_05_25_01:51:08', '2017_05_25_01:51:52', '2017_05_25_01:54:05', ...
%             '2017_05_25_01:54:49', '2017_05_25_01:55:35', '2017_05_25_01:56:19'};
%
%off_tstamp = {'2017_05_24_23:49:15', '2017_05_24_23:52:04', '2017_05_25_00:46:46', ...
%              '2017_05_25_00:49:45', '2017_05_25_00:52:44', '2017_05_25_00:55:43', ...
%              '2017_05_25_00:55:43', '2017_05_25_00:58:42', '2017_05_25_01:01:42', ...
%              '2017_05_25_01:30:31', '2017_05_25_01:36:22', '2017_05_25_01:47:07', ...
%               '2017_05_25_01:53:00'};
          
% bad_freqs = [91:100,191:200,291:300,391:400,491:500];
% good_idx = [20,21,23:33,35:39];
% proj_id = 'AGBT16B_400_01';
% save_dir = sprintf('%s/%s/BF', data_dir, proj_id);
% out_dir = sprintf('%s/mat', save_dir);
          

% May 26, 2017 Calibration Grid
on_tstamp = {'2017_05_27_04:34:08', '2017_05_27_04:35:01', ...
             '2017_05_27_04:36:51','2017_05_27_04:37:44', '2017_05_27_04:38:37', ...
             '2017_05_27_04:40:38','2017_05_27_04:41:31', '2017_05_27_04:42:24', ...
             '2017_05_27_04:44:12','2017_05_27_04:45:05', '2017_05_27_04:45:57', ...
             '2017_05_27_04:47:57','2017_05_27_04:48:49', '2017_05_27_04:49:41', ...
             '2017_05_27_04:51:29','2017_05_27_04:52:21', '2017_05_27_04:53:13', ...
             '2017_05_27_04:55:12','2017_05_27_04:56:03', '2017_05_27_04:56:54', ...
             '2017_05_27_04:58:41','2017_05_27_04:59:32', '2017_05_27_05:00:23', ...
             '2017_05_27_05:02:20','2017_05_27_05:03:11', '2017_05_27_05:04:02', ...
             '2017_05_27_05:05:47','2017_05_27_05:06:38', '2017_05_27_05:07:29', ...
             '2017_05_27_05:09:26','2017_05_27_05:10:16', '2017_05_27_05:11:06'};
         

off_tstamp = {'2017_05_27_04:36:09', '2017_05_27_04:39:49', '2017_05_27_04:43:31', ...
              '2017_05_27_04:47:09', '2017_05_27_04:50:48', '2017_05_27_04:54:24', ...
              '2017_05_27_04:58:00', '2017_05_27_05:01:33', '2017_05_27_05:05:07', ...
              '2017_05_27_05:08:39'};

bad_freqs = [21:25, 41:45, 81:100, 121:125, 141:145, 181:200, 221:225, 241:245, 281:300, 321:325, 341:345, 381:400, 421:425, 441:445, 481:500];
proj_id = 'AGBT16B_400_03';
save_dir = sprintf('%s/%s/BF', data_dir, proj_id);
out_dir = sprintf('%s/mat', save_dir);

% May 26, 2017 Calibration Grid
% on_tstamp = {'2017_05_28_19:53:44'};
% off_tstamp = {'2017_05_28_19:52:05', '2017_05_28_20:16:21'};
% 
% bad_freqs = [81:100, 181:200, 281:300, 381:400, 481:500];
% proj_id = 'AGBT16B_400_05';
% save_dir = sprintf('%s/%s/BF', data_dir, proj_id);
% out_dir = sprintf('%s/mat', save_dir);


% Constants
overwrite = 0;
k = 0;
good_idx = [20,21,23:33,35:39];
% good_idx = [1:19];
kB = 1.38*1e-23;

rad = 50;
Ap = (rad^2)*pi;

LO_freq = 1450;
freqs = ((-249:250)*303.75e-3) + LO_freq;
a = 1.4701;
b = -0.7658;
c = -0.2780;
d = -0.0347;
e = 0.0399;
x = a + b*log10(freqs./1e3) + c*log10(freqs./1e3).^2 + d*log10(freqs./1e3).^3 + e*log10(freqs./1e3).^4;
flux_density = 10.^x;

ant_dir = sprintf('%s/%s/Antenna', meta_dir, proj_id);

% Iterate over off pointings just to have them ready
% Iterate over on pointings and look for closest off pointing

% Iterate over off pointings
AZoff = zeros(length(off_tstamp), 1);
ELoff = zeros(length(off_tstamp), 1);
fprintf('Processing off pointings...\n');
for j = 1:length(off_tstamp)
    tmp_stmp = off_tstamp{j};
    fprintf('    Time stamp: %s - %d/%d\n', tmp_stmp, j, length(off_tstamp));
    
    % Generate filename
    filename = sprintf('%s/%s.mat', out_dir, tmp_stmp);
    
    % Extract data and save
    if ~exist(filename, 'file') || overwrite == 1
        [Rtmp, az_tmp, el_tmp, ~] = aggregate_banks(save_dir, ant_dir, tmp_stmp);
        
        % Off pointings are dwell scans; need single R, az, and el
        R = mean(Rtmp,4);
        az = mean(az_tmp);
        el = mean(el_tmp);
        save(filename, 'R', 'az', 'el');
    else
        load(filename);
    end
    
    % Create entry in position table
    AZoff(j) = az;
    ELoff(j) = el;
end
    
figure(1);
plot(AZoff, ELoff, 'x');

% Iterate over on pointings
fprintf('Processing on pointings...\n');
Sens = [];
AZ = [];
EL = [];
for i = 1:length(on_tstamp)
    tmp_stmp = on_tstamp{i};
    fprintf('    Time stamp: %s - %d/%d\n', tmp_stmp, i, length(on_tstamp));
    
    % Generate filename
    filename = sprintf('%s/%s.mat', out_dir, tmp_stmp);
    
    % Extract data and save
    if ~exist(filename, 'file') || overwrite == 1
        [R, az, el, ~] = aggregate_banks(save_dir, ant_dir, tmp_stmp);
        save(filename, 'R', 'az', 'el');
    else
        load(filename);
    end
    
    hold on;
    plot(az, el, '-b');
    hold off;
    drawnow;
    
    % Find nearest off poinitng
    for j = 1:length(off_tstamp)
        az_dist = az - AZoff(j);
        el_dist = el - ELoff(j);
        
        vector_distance(j) = mean(sqrt(az_dist.^2 + el_dist.^2));
    end
    
    [~, idx] = min(vector_distance);
    OFF = load(sprintf('%s/%s', out_dir, off_tstamp{idx}));
    
    % Compute max-SNR weights and compute sensitivity
    Senstmp = zeros(size(R,4),size(R,3));
    for t = 1:size(R,4)
        for b = 1:size(R,3)
            if sum(bad_freqs == b) == 0
                [a, ~] = eigs(R(good_idx, good_idx, b, t), OFF.R(good_idx, good_idx, b), 1);
                w = OFF.R(good_idx, good_idx, b)\a;
                w = w./(w'*a);
                Pon = w'*R(good_idx, good_idx, b, t)*w;
                Poff = w'*OFF.R(good_idx, good_idx, b)*w;
                SNR(t,b) = (Pon - Poff)/Poff;
                Senstmp(t,b) = 2*kB*SNR(t,b)./(flux_density(b)*1e-26);
            else
                Senstmp(t,b) = 0;
            end
        end
        AZ = [AZ; az(t)];
        EL = [EL; el(t)];
    end
    Sens = [Sens; Senstmp];
end

% Plot map

% Interpolated map
Npoints = 80;
minX = 212.3;
maxX = 213.4;
minY = 51.85;
maxY = 52.5;
xval = linspace(minX, maxX, Npoints);
yval = linspace(minY, maxY, Npoints);
[X,Y] = meshgrid(linspace(minX,maxX,Npoints), linspace(minY,maxY,Npoints));
map_fig = figure();
for b = 101:101 %Nbins
    fprintf('Bin %d/%d\n', b, 500);
    Sq = griddata(AZ, EL, real(squeeze(Sens(:,b))), X, Y);
    imagesc(xval, yval, Sq);
    colorbar;
    set(gca, 'ydir', 'normal');
    colormap('jet');
    xlabel('Azimuth Offset (degrees)');
    ylabel('Elevation Offset (degrees)');
    title('Formed Beam Sensitivity Map');
end

Tsys = Ap/max(max(Sens));

% Tsys_eta = Ap./Sens;

% Get the angle of arrival for max sensitivity beam
[s_max,max_idx] = max(Sens(:,101));
Tsys_eta = Ap./(Sens(max_idx,:)*22.4./flux_density); % Remember to change when running again.


figure(4);
plot(freqs,real(Tsys_eta).');
title('T_s_y_s/\eta_a_p vs Frequency');
xlabel('Frequency (MHz)');
ylabel('T_s_y_s/\eta_a_p (K)');
grid on;

% keyboard;

% % Iterate over off pointings
% for j = 1:length(off_tstamp)
%     for i = 1:length(on_tstamp)
%         k = k+1;
%         tmp_stmp = on_tstamp{i};
%         tmp_off_stmp = off_tstamp{j};
%         % % Generate filenames
%         Hot_file = sprintf('/home/groups/flag/mat/Hot_cov_%s.mat', tmp_stmp);
%         Cold_file = sprintf('/home/groups/flag/mat/Cold_cov_%s.mat', tmp_off_stmp); 
%         
%         % Get HOT data
%         if ~exist(Hot_file, 'file') || overwrite == 1
%             [R_hot, az, el, dmjd_hot] = aggregate_banks(save_dir, tmp_stmp);
%             save(Hot_file, 'R_hot');
%         else
%             load(Hot_file);
%         end
%         
%         % Get COLD data
%         if ~exist(Cold_file, 'file') || overwrite == 1
%             [R_cold, az_off, el_off, dmjd_cold] = aggregate_banks(save_dir, tmp_off_stmp);
%             save(Cold_file, 'R_cold');
%         else
%             load(Cold_file);
%         end
%         
%         Az(:,k) = az;
%         El(:,k) = el;
%         Az_off(:,k) = az_off;
%         El_off(:,k) = el_off;
%         
%         good_idx = [20,21,23:33,35:39];
%         for m = 1:size(R_cold,4)
%             for h = 1:size(R_hot,4)
%                 Ron = R_hot(good_idx, good_idx, :,h);
%                 Roff = R_cold(good_idx, good_idx, :,m);
%                 
%                 
%                 Nele = size(Ron, 1);
%                 Nele_act = size(R_hot,1);
%                 Nbins = size(Ron, 3);
%                 
%                 [w, Pon, Poff, SNR] = Multi_beam(Ron, Roff, h, Nele_act, Nbins, good_idx, bad_freqs, tmp_off_stmp, overwrite);
%                 
%                 Sens(k,h,:) = 2*kB*real(SNR)./(flux_density*1e-26);
%                 
%                 Tsys_eta(k,h,:) = Ap./Sens(k,h,:);
%             end
%         end
%     end
% end
% 
% 
% % figure(3);
% % imagesc(10*log10(squeeze(abs(Ron(:,:,5)))));
% % title('On covariance bin5');
% % xlabel('No. of elements');
% % ylabel('No. of elements');
% % 
% % figure(4);
% % imagesc(10*log10(squeeze(abs(Roff(:,:,5)))));
% %  title('Off covariance bin 5');
% % xlabel('No. of elements');
% % ylabel('No. of elements');
% % 
% % figure(5);
% % imagesc(10*log10(squeeze(abs(Ron(:,:,6)))));
% % title('On covariance bin 6');
% % xlabel('No. of elements');
% % ylabel('No. of elements');
% % 
% % figure(6);
% % imagesc(10*log10(squeeze(abs(Roff(:,:,6)))));
% % title('Off covariance bin 6');
% % xlabel('No. of elements');
% % ylabel('No. of elements');
% 
% LO_freq = 0;
% freqs = ((-249:250)*303.24e-3)+LO_freq;
% figure(7);
% % plot(1:Nbins,10*log10(Pon-Poff).');
% plot(freqs(1:Nbins),(Pon-Poff).');
% title('Beamformed Power');
% xlabel('Relative Frequency to LO (MHz)');
% % ylabel('Power (dB)');
% ylabel('Power');
% grid on;
% 
% figure(8);
% plot(freqs(1:Nbins),SNR);
% title('SNR vs Frequency');
% xlabel('Frequency');
% ylabel('SNR');
% grid on;
% 
% figure(9);
% plot(1:Nbins, squeeze(Sens(1,1,:)));
% % plot(freqs(1:Nbins), Sens(1,1,:));
% title('Sensitivity vs Frequency bins');
% xlabel('Frequency bins');
% ylabel('Sensitivity');
% grid on;
% 
% figure(10);
% plot(freqs(1:Nbins), squeeze(Tsys_eta(1,1,:)));
% title('Tsys/\eta_a_p vs Frequency bins');
% xlabel('Frequency bins');
% ylabel('Tsys/\eta_a_p');
% grid on;
% 
% figure(11);
% plot(1:Nbins, squeeze(Sens(1,:,:)));
% % plot(freqs(1:Nbins), Sens);
% title('Sensitivity map');
% xlabel('Frequency bins');
% ylabel('Sensitivity per beam');
% grid on;
% 
% figure(12);
% plot(1:Nbins, squeeze(Tsys_eta(1,:,:)));
% % plot(freqs(1:Nbins), Tsys_eta);
% title('Tsys/\eta_a_p map (K)');
% xlabel('Frequency (MHz)');
% ylabel('Tsys/\eta_a_p per beam');
% grid on;
% 
% 
% % Interpolated map
% Npoints = 80;
% use_radec = 0;
% if use_radec
%     minX = max(min(Az(:,1)), -0.45);
%     maxX = min(max(Az(:,1)), 0.45);
%     minY = max(min(El(:,1)), -0.35);
%     maxY = min(max(El(:,1)), 0.35);
% else
%     minX = min(Az(:,1));
%     maxX = max(Az(:,1));
%     minY = min(El(:,1));
%     maxY = max(El(:,1));
% end
% xval = linspace(minX, maxX, Npoints);
% yval = linspace(minY, maxY, Npoints);
% [X,Y] = meshgrid(linspace(minX,maxX,Npoints), linspace(minY,maxY,Npoints));
% for b = 101:101 %Nbins
%     fprintf('Bin %d/%d\n', b, Nbins);
%     map_fig = figure(13);
%     Sq = griddata(Az(:,1), El(:,1), real(squeeze(Sens(1,:,b))), X, Y);
%     imagesc(xval, yval, Sq);
%     colorbar;
%     set(gca, 'ydir', 'normal');
%     colormap('jet');
%     xlabel('Azimuth Offset (degrees)');
%     ylabel('Elevation Offset (degrees)');
%     title('Formed Beam Sensitivity Map');
% end
% 
% figure(14);
% for k = 1:size(Az,1)
%     plot(Az(k,:),El(k,:),'k-'); hold on;
% end
% title('Azimuth vs Elevation');
% xlabel('Azimuth Offset (degrees)');
% ylabel('Elevation Offset (degrees)');
% grid on;
% hold off;
% 

toc; 

% print(gcf, '-dpng', sprintf('fig/OTF_Tsys_%s.png', on_tstamp));


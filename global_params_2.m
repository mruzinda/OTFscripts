% Specifies common-to-all parameters

% data_set selects a data set to analyze
% 1 => May 20, 2017, OTF Scan, Quant gain = 10
% 2 => May 21, 2017, OTF Scan, Quant gain = 10
% 3 => May 21, 2017, OTF Scan, Quant gain = 10
% 4 => May 21, 2017, OTF Scan, Quant gain = 10
% 5 => May 21, 2017, OTF Scan, Quant gain = 20
% 6 => May 21, 2017, OTF Scan, Quant gain = 10
% 7 => May 21, 2017, OTF Scan, Quant gain = 10
% 8 => May 21, 2017, OTF Scan, Quant gain = 20
% 9 => May 21, 2017, OTF Scan, Quant gain = 40
data_set = 1;

% Specify if output files should be overwritten
overwrite = 1;

% Specify bank ID (only single-bank scripts use this)
bank = 'D';

% Create FITS file names
switch data_set
    case 1
        projID = 'TGBT16A_508_01';
        lustre_dir = '/lustre/gbtdata';
        sub_dir = 'TMP/BF';
        stamp = '2017_05_20_01:39:52';
        off_stamp = '2017_05_20_01:46:47';
        Thot = 290;
        Tcold = 7.5;
%         meta_projID = 'TGBT16A_508_02';
%         meta_stamp = '2016_07_25_04:32:33';
%         meta_off_projID = 'TGBT16A_508_02';
%         meta_off_stamp = '2016_07_25_04:32:33';
%         flux_density = 14.90; % 3C286 @ 1465 MHz
%         x_idx = [1:19];
%         y_idx = [21:39];
%         badX = [];
%         badY = [2,9,12,17];
%         idxs = [13:105, 209:1263];
%         off_idxs = [5130:5140];
%         center_freq = 1445.565e6;
%         dmjd_delay = 10.8/(24*3600);
%         dmjd_delay_off = dmjd_delay;
%         az_vals = [-0.08496, 0.09283, -0.1496,   0.003936,  0.1494,  -0.09304, 0.09283];
%         el_vals = [ 0.1639,  0.1639, 0.002164,  0.002164, 0.002164,  -0.1014, -0.1014];
%         use_radec = 0;
    case 2
        projID = 'TGBT16A_508_01';
        lustre_dir = '/lustre/gbtdata';
        sub_dir = 'TMP/BF';
        stamp = '2017_05_21_03:06:05';
        off_stamp = '2017_05_21_03:10:22';
        Thot = 290;
        Tcold = 7.5;
    case 3
        projID = 'TGBT16A_508_01';
        lustre_dir = '/lustre/gbtdata';
        sub_dir = 'TMP/BF';
        stamp = '2017_05_21_03:16:27';
        off_stamp = '2017_05_21_03:12:45';
        Thot = 290;
        Tcold = 7.5;
    case 4
        projID = 'TGBT16A_508_01';
        lustre_dir = '/lustre/gbtdata';
        sub_dir = 'TMP/BF';
        stamp = '2017_05_21_03:19:47';
        off_stamp = '2017_05_21_03:24:53';
        Thot = 290;
        Tcold = 7.5;
    case 5
        projID = 'TGBT16A_508_01';
        lustre_dir = '/lustre/gbtdata';
        sub_dir = 'TMP/BF';
        stamp = '2017_05_21_03:32:31';
        off_stamp = '2017_05_21_03:28:16';
        Thot = 290;
        Tcold = 7.5;
    case 6
        projID = 'TGBT16A_508_01';
        lustre_dir = '/lustre/gbtdata';
        sub_dir = 'TMP/BF';
        stamp = '2017_05_22_01:02:42';
        off_stamp = '2017_05_22_01:06:18';
        Thot = 290;
        Tcold = 7.5;
    case 7
        projID = 'TGBT16A_508_01';
        lustre_dir = '/lustre/gbtdata';
        sub_dir = 'TMP/BF';
        stamp = '2017_05_22_01:09:49';
        off_stamp = '2017_05_22_01:13:27';
        Thot = 290;
        Tcold = 7.5;
    case 8
        projID = 'TGBT16A_508_01';
        lustre_dir = '/lustre/gbtdata';
        sub_dir = 'TMP/BF';
        stamp = '2017_05_22_01:40:27';
        off_stamp = '2017_05_22_01:36:57';
        Thot = 290;
        Tcold = 7.5;
    case 9
        projID = 'TGBT16A_508_01';
        lustre_dir = '/lustre/gbtdata';
        sub_dir = 'TMP/BF';
        stamp = '2017_05_22_01:43:10';
        off_stamp = '2017_05_22_01:46:57';
        Thot = 290;
        Tcold = 7.5;
end

meta_dir = '/home/gbtdata';
ant_sub_dir = 'Antenna';

dir = sprintf('%s/%s/%s', lustre_dir, projID, sub_dir);
% ant_dir = sprintf('%s/%s/%s', meta_dir, meta_projID, ant_sub_dir);
% ant_off_dir = sprintf('%s/%s/%s', meta_dir, meta_off_projID, ant_sub_dir);
% ant_fits_file = sprintf('%s/%s.fits', ant_dir, meta_stamp);
% ant_off_fits_file = sprintf('%s/%s.fits', ant_dir, meta_off_stamp);

% % Specify which antenna elements will be used in beamformer and sensitivity
% % measurements (data_set dependent)
% all_idx = 1:40;
% goodX = all_idx(x_idx);
% goodX(badX) = [];
% goodY = all_idx(y_idx);
% goodY(badY) = [];

% % Select y-pol for now
% idx = goodY;

% Boltzmann's Constant
kb = 1.38e-23;

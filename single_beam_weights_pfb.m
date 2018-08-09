% Simple script that extracts the covariance matrices for the entire band
% and computes Y-factors

close all;
clearvars;

tic;
% System parameters
Dir = '/lustre/projects/flag/'; % '/lustre/projects/flag/TMP/BF'; %'/lustre/gbtdata/TGBT16A_508_01/TMP/BF';
projID = '/AGBT16B_400_09';
sub_dir = '/BF';
save_dir = sprintf('%s/%s/%s', Dir, projID, sub_dir);

% % May 24th, 2017 - 05  %%%%%% GBT test %%%%%%
% % Quant gain = 10
% on_tstamp = {'2017_05_28_05:39:03'};
% off_tstamp = {'2017_05_28_05:39:52'};

% July, 28th 2017 - 05  %%%%%% GBT test %%%%%%
% Quant gain = 
on_tstamp = {'2017_07_28_05:49:52'};
off_tstamp = {'2017_07_28_05:50:50'};

overwrite = 1;

AZoff = zeros(length(off_tstamp),1);
ELoff = zeros(length(off_tstamp),1);
good_idx = [1:7, 35, 9:19]; % [20, 21, 23:34, 8, 36:39];
bad_freqs = []; %[81:100,181:200,281:300,381:400,481:500];

% Off Pointing
fprintf('Getting OFF pointing from %s...\n', off_tstamp{1});
tmp_stmp = off_tstamp{1};
filename = sprintf('%s/mat/%s.mat', save_dir, tmp_stmp);

if ~exist(filename, 'file') || overwrite == 1
    [R, az_off, el_off] = aggregate_banks1_pfb(save_dir, projID, tmp_stmp);
    save(filename, 'R', 'az_off', 'el_off');
else
    load(filename);
end
Roff = R(good_idx, good_idx, :);
    

% ON pointing
fprintf('Getting ON pointing from %s...\n', on_tstamp{1});
tmp_stmp = on_tstamp{1};

filname = sprintf('%s/mat/%s.mat', save_dir, tmp_stmp);
if ~exist(filename, 'file') || overwrite == 1
    [R, az, el] = aggregate_banks1_pfb(save_dir, projID, tmp_stmp);
    save(filename, 'R', 'az', 'el');
else
    load(filename);
end

Ron = R(good_idx, good_idx, :);

Nele = size(Ron, 1);
Nele_act = size(R,1);
Nbins = size(Ron, 3);

[w, Pon, Poff, snr] = Single_beam_pfb(Ron, Roff, Nbins, good_idx, bad_freqs);

weights = zeros(64, 3200, 14);
for i = 1:14
    weights(:,:,i) = w;
end

% Reshape weights
N_beam = size(weights,3)/2;
N_ele= 64;
N_bin = 160;
N_pol = 2;

% Save data into weight file formatted for RTBF code
banks = {'A', 'B', 'C', 'D',...
    'E', 'F', 'G', 'H',...
    'I', 'J', 'K', 'L',...
    'M', 'N', 'O', 'P',...
    'Q', 'R', 'S', 'T'};

interleaved_w = zeros(2*N_ele*N_bin*N_beam*N_pol,1);
weight_dir = sprintf('%s/weight_files', Dir);
steer_dir = sprintf('%s/steer_files', Dir);
% Frequency index stitching
idx = 1:160;
idx1 = reshape(idx, [5,32]);
idx2 = idx1';
chan_idx=reshape(idx2,[160,1]);

for k = 1:5
    chan_idx((k-1)*32+1:k*32) = flip(fftshift(chan_idx((k-1)*32+1:k*32)));
end

for b = 1:length(banks)
    % Get bank name
    bank_name = banks{b};

    % Extract channels for bank
    w1 = weights(:,N_bin*(b-1)+1:N_bin*b,:);
    
    % Reshape for file format
    w2 = reshape(w1, N_ele*N_bin, N_beam*N_pol);
    w_real = real(w2(:));
    w_imag = imag(w2(:));
    interleaved_w(1:2:end) = w_real(:);
    interleaved_w(2:2:end) = w_imag(:);
    
    % Get filename
    weight_file = sprintf('%s/w_%s_%s.bin', weight_dir, tmp_stmp, bank_name);
    
    % Create metadata for weight file
    offsets_el = el;
    offsets_az = az;
    offsets = [offsets_el; offsets_az; offsets_el; offsets_az; offsets_el; offsets_az; offsets_el; offsets_az; offsets_el; offsets_az; offsets_el; offsets_az; offsets_el; offsets_az];
    offsets = offsets(:);
    cal_filename = sprintf('%s%s.fits',on_tstamp{1}, banks{b});
    to_skip1 = 64 - length(cal_filename);
    algorithm_name = 'Max Signal-to-Noise Ratio';
    to_skip2 = 64 - length(algorithm_name);
    xid = b-1;
    
    % Write to binary file
    WID = fopen(weight_file,'wb');
    if WID == -1
        error('Author:Function:OpenFile', 'Cannot open file: %s', weight_file);
    end
    
    % Write payload
    fwrite(WID, single(interleaved_w), 'single');
    
    % Write metadata
    fwrite(WID,single(offsets),'float');
    fwrite(WID,cal_filename, 'char*1');
    if to_skip1 > 0
        fwrite(WID, char(zeros(1,to_skip1)));
    end
    fwrite(WID,algorithm_name, 'char*1');
    if to_skip2 > 0
        fwrite(WID, char(zeros(1,to_skip2)));
    end
    fwrite(WID, uint64(xid), 'uint64');
    fclose(WID);
    
    fprintf('Saved to %s\n', weight_file);
end

toc; 


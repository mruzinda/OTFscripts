function [ R, my_az, my_el, info ] = aggregate_banks1_pfb( save_dir, projID, tstamp)
%AGGREGATE_BANKS Collects all the BANK data into a single covariance matrix
%variable
%   Takes the data from dir/tstamp and finds banks A-T. Opens them,
%   reconstructs the covariance matrices and aggregates them into a
%   40x40x500 matrix.

    % Constants
    banks = {'A', 'B', 'C', 'D',...
        'E', 'F', 'G', 'H',...
        'I', 'J', 'K', 'L',...
        'M', 'N', 'O', 'P',...
        'Q', 'R', 'S', 'T'};
    
    % Frequency index stitching
    idx = 1:160;
    idx1 = reshape(idx, [5,32]);
    idx2 = idx1';
    chan_idx=reshape(idx2,[160,1]);

    for k = 1:5
        chan_idx((k-1)*32+1:k*32) = flip(fftshift(chan_idx((k-1)*32+1:k*32)));
    end
    
    Nele = 40;
    Nchan = 160;
    Nchan_tot = 3200;
    R = zeros(Nele,Nele,Nchan_tot);
    fprintf('     ');
    for i = 1:length(banks)
        filename = sprintf('%s/%s%s.fits', save_dir, tstamp, banks{i});
        % Check to see that file exists
        if exist(filename, 'file') ~= 0
            [Rtmp, dmjd, xid, info] = extract_pfb_covariances(filename);
            Ravg = mean(Rtmp, 4);
            R(:,:,chan_idx + Nchan*xid) = Ravg(1:Nele, 1:Nele, :);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract antenna positions for scan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Extracting antenna positions and computing offsets');

% Extract offsets
ant_dir = sprintf('/home/gbtdata/%s/Antenna', projID);
use_radec = 1;
ant_fits_file = sprintf('%s/%s.fits', ant_dir, tstamp);
fprintf('Reading from %s\n', ant_fits_file);
[ant_dmjd, az_off, el_off] = get_antenna_positions(ant_fits_file, use_radec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Associate offsets with correlation matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Associating pointing offsets with correlation matrices');
Ntime = 1;
my_el  = zeros(Ntime, 1);
my_az  = zeros(Ntime, 1);
my_ra  = zeros(Ntime, 1);
my_dec = zeros(Ntime, 1);
for t = 1:Ntime
    cur_dmjd = dmjd(t);
    tmp_dmjd = ant_dmjd;
    [~, idx1] = min(abs(tmp_dmjd - cur_dmjd));
    
    tmp_dmjd(idx1) = NaN;
    [~, idx2] = min(abs(tmp_dmjd - cur_dmjd));
    
    x1 = min(ant_dmjd(idx1), ant_dmjd(idx2));
    x2 = max(ant_dmjd(idx1), ant_dmjd(idx2));
    az1 = min(az_off([idx1,idx2]));
    az2 = max(az_off([idx1,idx2]));
    el1 = min(el_off([idx1,idx2]));
    el2 = max(el_off([idx1,idx2]));
    
    my_el(t) = (el2 - el1)/(x2 - x1)*(cur_dmjd - x1) + el2;
    my_az(t) = (az2 - az1)/(x2 - x1)*(cur_dmjd - x1) + az2;
    my_ra(t) = (el2 - el1)/(x2 - x1)*(cur_dmjd - x1) + el2;
    my_dec(t) = (az2 - az1)/(x2 - x1)*(cur_dmjd - x1) + az2;
end

end

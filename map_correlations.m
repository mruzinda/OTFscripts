function [R, az_off, el_off] = map_correlations(projID, stamp, ant_fits_file, use_radec)

lustre_dir = '/lustre/gbtdata';
sub_dir = 'TMP/BF';
dir = sprintf('%s/%s/%s', lustre_dir, projID, sub_dir);

fprintf('Aggregating and extracting covariances for %s\n', stamp);
[Rtmp, dmjd] = aggegrate_banks(dir, stamp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ON Pointing Aggregation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply additional integration
Ntime = size(Rtmp,4);
R = zeros(size(Rtmp,1), size(Rtmp,2), size(Rtmp,3), floor(Ntime/4)-1);
dmjd = dmjd + dmjd_delay;
% Correct for dmjd skew from time
corrected_dmjd = (15150/15187.5)*(dmjd - dmjd(1)) + dmjd(1);

data_dmjd = zeros(floor(Ntime/4)-1,1);
for j = 1:floor(Ntime/4)-1
    R(:,:,:,j) = sum(Rtmp(:,:,:,4*(j-1)+1:4*j),4)/4;
    data_dmjd(j) = sum(corrected_dmjd(4*(j-1)+1:4*j))/4;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract antenna positions for scan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Extracting antenna positions and computing offsets');

% Extract offsets
% use_radec = 0; % Flag for azimuth and elevation calculation using ra and dec
[ant_dmjd, az_off, el_off, ra, dec] = get_antenna_positions(ant_fits_file, use_radec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Associate offsets with correlation matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Associating pointing offsets with correlation matrices');
Ntime = size(R, 4);
my_el  = zeros(Ntime, 1);
my_az  = zeros(Ntime, 1);
my_ra  = zeros(Ntime, 1);
my_dec = zeros(Ntime, 1);
for t = 1:Ntime
    cur_dmjd = data_dmjd(t);
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
    ra1 = min(ra([idx1,idx2]));
    ra2 = max(ra([idx1,idx2]));
    dec1 = min(dec([idx1,idx2]));
    dec2 = max(dec([idx1,idx2]));
    
    my_el(t) = (el2 - el1)/(x2 - x1)*(cur_dmjd - x1) + el2;
    my_az(t) = (az2 - az1)/(x2 - x1)*(cur_dmjd - x1) + az2;
    my_ra(t) = (el2 - el1)/(x2 - x1)*(cur_dmjd - x1) + el2;
    my_dec(t) = (az2 - az1)/(x2 - x1)*(cur_dmjd - x1) + az2;
end

ANT.az_off = my_az;
ANT.el_off = my_el;
ANT.ra = my_ra;
ANT.dec = my_dec;
ANT.dmjd = data_dmjd;

end
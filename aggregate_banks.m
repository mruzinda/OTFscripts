function [ R, my_az, my_el, info ] = aggregate_banks( save_dir, ant_dir, tstamp, Nint)
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
    
    chan_idx = [1:5, 101:105, 201:205, 301:305, 401:405];
    % chan_idx = [1:25];

    Nele = 40;
    Nchan = 500;
    Rtmp1 = zeros(Nele,Nele,Nchan);
    fprintf('     ');
    for i = 1:length(banks)
        filename = sprintf('%s/%s%s.fits', save_dir, tstamp, banks{i});
        fprintf('%c', banks{i});
        % Check to see that file exists
        if exist(filename, 'file') ~= 0
            [Rtmp, dmjd, xid, info] = extract_covariances(filename);
    %         Ravg = mean(Rtmp, 4);
    %         R(:,:,chan_idx + 5*xid) = Ravg(1:Nele, 1:Nele, :);
            if i == 1
                Rtmp1 = zeros(Nele, Nele, Nchan, size(Rtmp,4));
            end
            if size(Rtmp,4) ~= 0
                Rtmp1(:,:,chan_idx + 5*xid,:) = Rtmp;
            else
                fprintf('!');
            end

        end
    end

    N_time = size(Rtmp1,4);
    if Nint == -1
       Nint = N_time; 
    end
    R = zeros(size(Rtmp1,1), size(Rtmp1,2), size(Rtmp1,3), floor(N_time/Nint)-1);
    dec_dmjd = zeros(floor(N_time/Nint), 1);
    for j = 1:floor(N_time/Nint)-1
        R(:,:,:,j) = sum(Rtmp1(:,:,:,Nint*(j-1)+1:Nint*j),4)/Nint;
        dec_dmjd(j) = sum(dmjd(Nint*(j-1)+1:Nint*j))/Nint;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract antenna positions for scan
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Extracting antenna positions and computing offsets');

    % Extract offsets
    % ant_dir = '/home/gbtdata/AGBT16B_400_05/Antenna';
    use_radec = 1;
    ant_fits_file = sprintf('%s/%s.fits', ant_dir, tstamp);
    [ant_dmjd, az_off, el_off] = get_antenna_positions(ant_fits_file, use_radec);

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
        cur_dmjd = dec_dmjd(t);
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


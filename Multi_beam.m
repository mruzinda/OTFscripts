function [w, Pon, Poff, SNR] = Multi_beam(Ron, Roff, h, Nele, Nbins, good_idx, bad_freqs, off_tstamp, overwrite)
v = zeros(length(good_idx), Nbins);
w = zeros(length(good_idx), Nbins);
wei = zeros(Nele, Nbins);

for b = 1:size(Ron,3)
    bad_flag = sum(bad_freqs == b);
    if bad_flag
        w(:,b) = zeros(length(good_idx), 1);
    else
        [v(:,b),lambda] = eigs(Ron(:,:,b), Roff(:,:,b), 1);
        w(:,b) = Roff(:,:,b)\v(:,b);
        w(:,b) = w(:,b)./(w(:,b)'*v(:,b));
    end
    Pon(b) = w(:,b)'*Ron(:,:,b)*w(:,b);
    Poff(b) = w(:,b)'*Roff(:,:,b)*w(:,b);
    SNR(b) = (Pon(b) - Poff(b))/Poff(b);
    wei(good_idx,b) = w(:,b);
end

weight_dir = '/home/groups/flag/weight_files';
tmp_w = wei.';
w_real = real(tmp_w(:));
w_imag = imag(tmp_w(:));
interleaved_w = zeros(2*Nele*Nbins,1);
interleaved_w(1:2:end) = w_real;
interleaved_w(2:2:end) = w_imag;
weight_file = sprintf('%s/weights_%dbeam_%s.in', weight_dir, h, off_tstamp);
if ~exist(weight_file) || overwrite == 1
    save(weight_file,'interleaved_w');
end

end
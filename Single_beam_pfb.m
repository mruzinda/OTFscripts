function [wei, Pon, Poff, SNR] = Single_beam_pfb(Ron, Roff, Nbins, good_idx, bad_freqs)
v = zeros(length(good_idx), Nbins);
w = zeros(length(good_idx), Nbins);
Nbin_Act = 3200;
Nele_Act = 64;
wei = zeros(Nele_Act, Nbin_Act);
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
end





w_file = '/lustre/projects/flag/weight_files/wA.bin';
w = fopen(w_file,'rb');
wei = fread(w, 2*7*2*25*64, 'single');
weicomp = wei(1:2:end) + 1j*wei(2:2:end);
weights = reshape(weicomp, 64, 25, 14);

azoffset1 = fread(w, 1, 'single');
eloffset1 = fread(w, 1, 'single');
disp([azoffset1, eloffset1]);

azoffset1 = fread(w, 1, 'single');
eloffset1 = fread(w, 1, 'single');
disp([azoffset1, eloffset1]);

azoffset1 = fread(w, 1, 'single');
eloffset1 = fread(w, 1, 'single');
disp([azoffset1, eloffset1]);

azoffset1 = fread(w, 1, 'single');
eloffset1 = fread(w, 1, 'single');
disp([azoffset1, eloffset1]);

azoffset1 = fread(w, 1, 'single');
eloffset1 = fread(w, 1, 'single');
disp([azoffset1, eloffset1]);

azoffset1 = fread(w, 1, 'single');
eloffset1 = fread(w, 1, 'single');
disp([azoffset1, eloffset1]);

azoffset1 = fread(w, 1, 'single');
eloffset1 = fread(w, 1, 'single');
disp([azoffset1, eloffset1]);


figure(1);
plot(squeeze(abs(weights(:,2,6))));

fclose(w);

function shFPM_getSubimage()
%{
    FPM simulation 
    generate low resolution image
%}


%% read object imaging
pix = 1540;
I = imread('resource\\intensity.png');
I = double(I(:,:,1));I = I-min(I(:));I = I/max(max(I));I = imresize(I,[pix,pix]); %intensity
P = imread('resource\\phase.jpg');
P = double(P(:,:,1));P = P-min(P(:));P = P/max(max(P));P = imresize(P,[pix,pix]);
O = sqrt(I).*exp(1i*pi*P); % complex amplitude of object

%% initializing enviroment
[lambda,n_LED,~,CTF_object,N_lens,~,~,...
                     pix_CCD,sub_pix,plane_wave,df]=ini_enviroment();

%% Simulation of the process of snapshot
F = fftshift(fft2(O));
I_camera = zeros(sub_pix*N_lens,sub_pix*N_lens,n_LED^2);


for con = 1:n_LED^2
    fxc = round((pix+1)/2+(plane_wave(1,con)/lambda)/df);
    fyc = round((pix+1)/2+(plane_wave(2,con)/lambda)/df);
    
    fxl=round(fxc-(pix_CCD-1)/2);fxh=round(fxc+(pix_CCD-1)/2);
    fyl=round(fyc-(pix_CCD-1)/2);fyh=round(fyc+(pix_CCD-1)/2);
    
    F_sub = F(fyl:fyh,fxl:fxh) .* CTF_object; 
    I_camera(:,:,con) = abs(arrayfft2(F_sub,N_lens,N_lens,sub_pix,'ifft')).^2;
    con
end

%% intensity in CCD
I_camera = I_camera-min(I_camera(:));
I_camera = I_camera/max(I_camera(:));
save('resource//Lseq.mat','I_camera')
close all
end
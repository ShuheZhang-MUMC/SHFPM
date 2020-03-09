%% Initialize experimental parameter
%{
    parameters

    CTF_object0: pupil function without aberration
    CTF_object : pupil function with aberration
    
    N_lens     : Row number of microlens array
    CTF_subLA  : sub-pupil of lenslet
    NA         : numerical aperture
    pix_CCD    : pixel number of camera sensor
    sub_pix    : pixel number of sub-images
    plane_wave : illumination plane_wave
    df         : step length of frequency for each pixel
%}


function [lambda,n_LED,CTF_object0,CTF_object,N_lens,CTF_subLA,NA,...
                     pix_CCD,sub_pix,plane_wave,df]=ini_enviroment()
  clc
  clear
%% Objective properties
lambda = 0.532; % wavelength um
k = 2*pi/lambda;
NA = 0.15; %nurmical aperture


%% Lens array properties
N_lens = 11;     % row number of Lens array
sub_pix = 35;    % corresponding pixel on CCD of each lens
[tempfx,tempfy] = meshgrid(linspace(-1,1,sub_pix));
CUT_OFF_LA = 0.99; % ratio of lenslet aperture in each square grid
CTF_subLA = sqrt(tempfx.^2 + tempfy.^2) <= CUT_OFF_LA;

%% CCD properties
Mag = 4;
pix_CCD = sub_pix * N_lens;
pix_SIZ = 6.4;
sample_size = pix_SIZ/Mag * pix_CCD;

fx_CCD = (-pix_CCD/2:pix_CCD/2-1)/(pix_CCD*pix_SIZ/Mag);
df = fx_CCD(2)-fx_CCD(1);
[fx_CCD,fy_CCD] = meshgrid(fx_CCD);
CTF_CCD = (fx_CCD.^2+fy_CCD.^2)<(NA/lambda).^2;


CTF_LA  = ones(pix_CCD);
for x_con = 1:N_lens
    for y_con = 1:N_lens
        CTF_LA(sub_pix*(y_con-1)+1:sub_pix*y_con,...
            sub_pix*(x_con-1)+1:sub_pix*x_con)= CTF_subLA;
    end
end
CTF_object0 = CTF_LA.*CTF_CCD;

%% aberration 
CTF_CCD = CTF_CCD.*exp(1i*12*atan2(fy_CCD,fx_CCD)).*exp(1i*200*sqrt(1-fy_CCD.^2-fx_CCD.^2));
CTF_object = CTF_LA.*CTF_CCD;

figure();
imshow(abs(CTF_object),[]); %modulus of pupil function 
figure();
imshow(mod(angle(CTF_object),2*pi),[]) %angle of pupil function 

P = mod(angle(CTF_CCD),2*pi);
P = P - min(P(:));
P = P / max(P(:));
imwrite(P,'output//or_mod_20200225.png')

df_LA = df * 64;
Ratio_LED_LA = 0.6;

%% LED properties
d_LED = 4    *1000; % distance between adjust LED dot
h_LED = sqrt((N_lens*pix_SIZ*d_LED/(lambda*Mag*Ratio_LED_LA)).^2-d_LED^2)  % distance between LED matrix and sample

n_LED = 13;     % number of LED
x_LED = -(n_LED/2*d_LED-d_LED/2):d_LED:(n_LED/2*d_LED-d_LED/2);
[x_LED,y_LED] = meshgrid(x_LED);
[xpos,ypos] = getoddseq(n_LED);

R0 = lambda*(d_LED/(lambda*sqrt(d_LED^2+h_LED^2)))/NA;
R_overlap = 1/pi*(2*acos(R0/2)-R0*sqrt(1-(R0/2)^2)) %R_overlap


plane_wave = zeros(2,n_LED^2); %tilted plane wave

for con = 1:n_LED^2
    v = [0,0,h_LED]-[x_LED(ypos(con),xpos(con)),y_LED(ypos(con),xpos(con)),0];
    v = v/norm(v);
    plane_wave(1,con)=v(1);
    plane_wave(2,con)=v(2);
end

end

function [xseq,yseq] = getoddseq(N)
if mod(N,2)==0
    error('N must be an odd integer');
end

num = N^2;
count = 1;
xseq = zeros(1,num);
yseq = zeros(1,num);
xseq(1,1) = (N+1)/2;
yseq(1,1) = (N+1)/2;
for n = 1:num
    if mod(n,2)==0
    for m = 1:n
        count = count + 1;
        xseq(1,count) = xseq(1,count-1);
        yseq(1,count) = yseq(1,count-1) -1;
        if count == num
            return;
        end
    end
    for m = 1:n
        count = count + 1;
        xseq(1,count) = xseq(1,count-1) -1;
        yseq(1,count) = yseq(1,count-1);
        if count == num
            return;
        end
    end
    else
        for m = 1:n
            count = count + 1;
            xseq(1,count) = xseq(1,count-1);
            yseq(1,count) = yseq(1,count-1) +1;
            if count == num
                return;
            end
        end
        for m = 1:n
            count = count + 1;
            xseq(1,count) = xseq(1,count-1) +1;
            yseq(1,count) = yseq(1,count-1);
            if count == num
                return;
            end
        end
    end
end

end
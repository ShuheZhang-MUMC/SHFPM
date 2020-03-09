function [lambda,n_LED,CTF_object0,CTF_object,...
                     pix_CCD,plane_wave,df]=ini_enviroment()
  clc
  clear
%% Objective properties
lambda = 0.532; % wavelength um
k = 2*pi/lambda;
NA = 0.15; %nurmical aperture




%% CCD properties
Mag = 4;
pix_CCD = 385;
pix_SIZ = 6.4;

fx_CCD = (-pix_CCD/2:pix_CCD/2-1)/(pix_CCD*pix_SIZ/Mag);
df = fx_CCD(2)-fx_CCD(1);
[fx_CCD,fy_CCD] = meshgrid(fx_CCD);
CTF_CCD = (fx_CCD.^2+fy_CCD.^2)<(NA/lambda).^2;

%% aberration 
CTF_object0 = CTF_CCD;

CTF_CCD = CTF_object0.*exp(1i*12*(atan2(fy_CCD,fx_CCD))+pi/3).*exp(1i*0*sqrt(1-fy_CCD.^2-fx_CCD.^2));%.*exp(-fr.^2/0.8^2).*exp(1i*1.3*pi*exp(-fr.^2/0.4^2));%
CTF_object = CTF_CCD;

figure();
imshow(abs(CTF_object),[])
figure();
imshow(mod(angle(CTF_object),2*pi),[])


%% LED properties
d_LED = 4    *1000;%sqrt(df_LED^2*lambda^2*h_LED^2/(1-df_LED^2*lambda^2)); % distance between adjust LED dot
h_LED = 250  *1000;%sqrt((N_lens*pix_SIZ*d_LED/(lambda*Mag*Ratio_LED_LA)).^2-d_LED^2) % distance between LED matrix and sample
n_LED = 13;     % number of LED
x_LED = -(n_LED/2*d_LED-d_LED/2):d_LED:(n_LED/2*d_LED-d_LED/2);
[x_LED,y_LED] = meshgrid(x_LED);
[xpos,ypos] = getoddseq(n_LED);
% R0 = lambda*(d_LED/(lambda*sqrt(d_LED^2+h_LED^2)))/NA;
R0 =d_LED/sqrt(d_LED^2+h_LED^2)/NA;
R_overlap = 1/pi*(2*acos(R0/2)-R0*sqrt(1-(R0/2)^2))

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

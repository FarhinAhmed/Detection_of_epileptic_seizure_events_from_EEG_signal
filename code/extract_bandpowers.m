clc
clear all
close all
load filt_chb01_03.mat
load chb01_03.mat
load datacut_01_03_ch1.mat
[r,c]=size(cut_data);
fs = 256;
%% entropy
% [se,te] = pentropy(record(1,:),fs);
% se2=[];
% for i=1:r
%     
%    [se2] =[se2 pentropy(cut_data(i,:),1,'Instantaneous',false)];
% end
%     
% subplot(2,1,1);plot(te,se,'linewidth',2)
% title('entire record')
% ylabel('Spectral Entropy')
% xlabel('Time')
% subplot(2,1,2); plot(se2,'linewidth',2)

%% Correlation coefficients
% p=[];
% for i=1:r-1
%     [p]=[p corr(cut_data(i,:)',cut_data(i+1,:)')]; % not good
% end
% plot(p(2990:3100))

%% FFT
% Ts=1/fs;
% filtered_eeg=cut_data(3036,:);
% t=0:Ts:length(filtered_eeg);
% X=fft(filtered_eeg);
% L=length(filtered_eeg);
% Y=fft(filtered_eeg);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = fs*(0:(L/2))/L;
% plot(f,P1) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
%% PSD
% pxx1=pwelch(cut_data(1000,:));
% pxx2=pwelch(cut_data(2998,:));
% pxx3=pwelch(cut_data(3020,:));

% composite_non=channel_record(1:40);
% composite_seizures=channel_record(2996:3036); % this part not imp
% figure
% cn=pwelch(composite_non);
% cs=pwelch(composite_seizures);
% plot(10*log10(cn));
% hold on
% plot(10*log10(cs));
% figure

% plot(10*log10(pxx1),'linewidth',2);
% xlim([0 50])
% hold on
% plot(10*log10(pxx2),'linewidth',2);
% xlim([0 50])
% hold on
% plot(10*log10(pxx3),'linewidth',2);
% xlim([0 50])
% legend('no seizure','seizure 1', 'seizure 2')
%% Energy
% e1=[];
% e2=[];
% for i=1:40
%   [e1]=[e1 sum(abs(cut_data(i,:).^2))];
% end
% for j=2996:3036
%   [e2]=[e2 sum(abs(cut_data(j,:).^2))];
% end
% plot(e1,'linewidth',2);
% hold on
% plot(e2,'linewidth',2);
% legend('normal','seizure')
% % e2=sum(abs(cut_data(2998,:).^2));
% % e3=sum(abs(cut_data(3020,:).^2));
%% band characteristics

% 
%delta power

Ts=1/fs;
delta=[];
indn=1:3600;
inds=2996:3036;
for i=indn(1):indn(end)
 t=0:Ts:length(cut_data(i,:));
 X=fft(cut_data(i,:));
 L=length(cut_data(i,:));
 Y=fft(cut_data(i,:));
 P2 = abs(Y/L);
 P1 = P2(1:L/2+1);
 P1(2:end-1) = 2*P1(2:end-1);
 f = fs*(0:(L/2))/L;

 f1=find(f>=0.5 & f<=4);
%  f2=find(f>=-4 & f<=-0.5);
 for ii = 1:length(P1)
        if ((ii>=f1(1) & ii<=f1(end)))
            t=1;
        end
        if (t==0)
            P1(ii)=0;

        end
        t=0;
  end
    P1=ifftshift(P1);
    s=real(ifft(P1));
    [delta]=[delta sum(s.^2)/length(s)];
end
subplot(4,1,1);plot(delta,'linewidth',2)
title("Delta Band Power")
hold on
Ts=1/fs;
deltas=[];
for i=2996:3036
 t=0:Ts:length(cut_data(i,:));
 X=fft(cut_data(i,:));
 L=length(cut_data(i,:));
 Y=fft(cut_data(i,:));
 P2 = abs(Y/L);
 P1 = P2(1:L/2+1);
 P1(2:end-1) = 2*P1(2:end-1);
 f = fs*(0:(L/2))/L;
% P1 is the fft

%delta energy
 
 f1=find(f>=0.5 & f<=4);
%  f2=find(f>=-4 & f<=-0.5);
 for ii = 1:length(P1)
        if ((ii>=f1(1) & ii<=f1(end)))
            t=1;
        end
        if (t==0)
            P1(ii)=0;

        end
        t=0;
  end
    P1=ifftshift(P1);
    s=real(ifft(P1));
    [deltas]=[deltas sum(s.^2)/length(s)];
end
plot(deltas,'linewidth',2)
legend('normal','seizure')
% 
 

Ts=1/fs;
theta=[];
for i=indn(1):indn(end)
 t=0:Ts:length(cut_data(i,:));
 X=fft(cut_data(i,:));
 L=length(cut_data(i,:));
 Y=fft(cut_data(i,:));
 P2 = abs(Y/L);
 P1 = P2(1:L/2+1);
 P1(2:end-1) = 2*P1(2:end-1);
 f = fs*(0:(L/2))/L;

%P1 is the fft

%theta energy
 
 f1=find(f>=4 & f<=8);
%  f2=find(f>=-4 & f<=-0.5);
 for ii = 1:length(P1)
        if ((ii>=f1(1) & ii<=f1(end)))
            t=1;
        end
        if (t==0)
            P1(ii)=0;

        end
        t=0;
  end
    P1=ifftshift(P1);
    s=real(ifft(P1));
    [theta]=[theta sum(s.^2)/length(s)];
end
subplot(4,1,2);plot(theta,'linewidth',2)
title("Theta Band Power")
hold on
Ts=1/fs;
thetas=[];
for i=2996:3036
 t=0:Ts:length(cut_data(i,:));
 X=fft(cut_data(i,:));
 L=length(cut_data(i,:));
 Y=fft(cut_data(i,:));
 P2 = abs(Y/L);
 P1 = P2(1:L/2+1);
 P1(2:end-1) = 2*P1(2:end-1);
 f = fs*(0:(L/2))/L;

%P1 is the fft

%theta energy
 
 f1=find(f>=4 & f<=8);
%  f2=find(f>=-4 & f<=-0.5);
 for ii = 1:length(P1)
        if ((ii>=f1(1) & ii<=f1(end)))
            t=1;
        end
        if (t==0)
            P1(ii)=0;

        end
        t=0;
  end
    P1=ifftshift(P1);
    s=real(ifft(P1));
    [thetas]=[thetas sum(s.^2)/length(s)];
end
plot(thetas,'linewidth',2)
legend('normal','seizure')
% 
Ts=1/fs;
alph=[];
for i=indn(1):indn(end)
 t=0:Ts:length(cut_data(i,:));
 X=fft(cut_data(i,:));
 L=length(cut_data(i,:));
 Y=fft(cut_data(i,:));
 P2 = abs(Y/L);
 P1 = P2(1:L/2+1);
 P1(2:end-1) = 2*P1(2:end-1);
 f = fs*(0:(L/2))/L;

% P1 is the fft

% alph power
 
 f1=find(f>=8 & f<=13);
 for ii = 1:length(P1)
        if ((ii>=f1(1) & ii<=f1(end)))
            t=1;
        end
        if (t==0)
            P1(ii)=0;

        end
        t=0;
  end
    P1=ifftshift(P1);
    s=real(ifft(P1));
    [alph]=[alph sum(s.^2)/length(s)];
end
subplot(4,1,3); plot(alph,'linewidth',2)
title("Alpha Band Power")

hold on
Ts=1/fs;
alphs=[];
for i=2996:3036
 t=0:Ts:length(cut_data(i,:));
 X=fft(cut_data(i,:));
 L=length(cut_data(i,:));
 Y=fft(cut_data(i,:));
 P2 = abs(Y/L);
 P1 = P2(1:L/2+1);
 P1(2:end-1) = 2*P1(2:end-1);
 f = fs*(0:(L/2))/L;

%P1 is the fft

%alpha energy
 
 f1=find(f>=8 & f<=13);
 f2=find(f>=-4 & f<=-0.5);
 for ii = 1:length(P1)
        if ((ii>=f1(1) & ii<=f1(end)))
            t=1;
        end
        if (t==0)
            P1(ii)=0;

        end
        t=0;
  end
    P1=ifftshift(P1);
    s=real(ifft(P1));
    [alphs]=[alphs sum(s.^2)/length(s)];
end
plot(alphs,'linewidth',2)
legend('normal','seizure')

%beta power
Ts=1/fs;
beta=[];
for i=indn(1):indn(end)
 t=0:Ts:length(cut_data(i,:));
 X=fft(cut_data(i,:));
 L=length(cut_data(i,:));
 Y=fft(cut_data(i,:));
 P2 = abs(Y/L);
 P1 = P2(1:L/2+1);
 P1(2:end-1) = 2*P1(2:end-1);
 f = fs*(0:(L/2))/L;

%P1 is the fft

%Beta energy
 
 f1=find(f>=13 & f<=30);
%  f2=find(f>=-4 & f<=-0.5);
 for ii = 1:length(P1)
        if ((ii>=f1(1) & ii<=f1(end)))
            t=1;
        end
        if (t==0)
            P1(ii)=0;

        end
        t=0;
  end
    P1=ifftshift(P1);
    s=real(ifft(P1));
    [beta]=[beta sum(s.^2)/length(s)];
end
subplot(4,1,4);plot(beta,'linewidth',2) 
title("Beta band Power")
hold on
Ts=1/fs;
betas=[];
for i=2996:3036
 t=0:Ts:length(cut_data(i,:));
 X=fft(cut_data(i,:));
 L=length(cut_data(i,:));
 Y=fft(cut_data(i,:));
 P2 = abs(Y/L);
 P1 = P2(1:L/2+1);
 P1(2:end-1) = 2*P1(2:end-1);
 f = fs*(0:(L/2))/L;

%P1 is the fft

%beta energy
 
 f1=find(f>=13 & f<=30);
%  f2=find(f>=-4 & f<=-0.5);
 for ii = 1:length(P1)
        if ((ii>=f1(1) & ii<=f1(end)))
            t=1;
        end
        if (t==0)
            P1(ii)=0;

        end
        t=0;
  end
    P1=ifftshift(P1);
    s=real(ifft(P1));
    [betas]=[betas sum(s.^2)/length(s)];
end
plot(betas,'linewidth',2)
legend('normal','seizure')


band_features_n=[delta' theta' alph' beta'];
band_features_s=[deltas' thetas' alphs' betas'];
save('band_features01_03_ch1','band_features_n','band_features_s')

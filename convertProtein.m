% Created 2025-03-15
% Authors: Jun Steed Huang et al
% Convert the Influenza A viral hemagglutinin H1, H3 and H5 proteins sequence for analysis  
  clc;
  clear;
  clf;
% read data from Excelsheet InfluenzaAhemagNewYorkH3
  [NUM,TXT,RAW]=xlsread('InfluenzaAhemagNewYorkH3.xls', 1, 'A1:BH12');
% convert table to string
  XTX=char(TXT');
% convert letter to number (low case or up case of your choice)
  letter2num=inline('x-''a''+1');
  NTX=letter2num(XTX);
% prepare charge acid mass hydro tables
  QT=[0,0,0,-1,-1,0,0,1,0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0];
  PT=[6.11,4.13,5.05,2.85,3.15,5.49,6.06,7.6,6.05,6.03,9.6,6.01,5.74,5.41,5.855,6.3,5.65,10.76,5.68,5.6,5.8,6,5.89,5.765,5.64,4.4];
  WT=[89.1,132.6,121.2,133.1,147.1,165.2,75.1,155.2,131.2,131.2,146.2,131.2,149.2,132.1,255.3,115.1,146.2,174.2,105.1,119.1,168.1,117.1,204.2,192.7,181.2,146.7];
  HT=[41,-41.5,49,-55,-31,100,0,8,99,98,-23,97,74,-28,-37,-46,-10,-14,-5,13,44.5,76,97,80,63,-20.5];
% look up tables
  [SZ,UN]=size(NTX);
for i=1:SZ
    t=NTX(i);
    if t<0
        break;
    else
    Q(i)=QT(t);
    P(i)=PT(t);
    W(i)=WT(t);
    H(i)=HT(t);
    end
end
% row to column transpose
QF=Q';
PF=P';
WF=W';
HF=H';
% plot raw data
figure(1);
plot(QF);
hold on;
plot(PF);
plot(WF);
plot(HF);
hold off;
grid on;
title('Flu Sequence Values');
legend('charge','isoelectric','mass','hydroph');
% calculate Langlands spectrum
[U,N]=size(Q);
i=1;
% z steps
for z = 0.1:0.01:10
    recordZ(i)= z;
    % sub 11th semi variance for acid by acid
for j=1:N    
    record(j)= ((W(j)+z*H(j))/(Q(j)+z*P(j)) - mean((W+z*H)/(Q+z*P)))^(1/11);
end
    record2(i)=mean(record);
    i=i+1;
end
figure(2);
loglog(recordZ, abs(record2));
title('Log Amplitude of Spectrum');
grid on;
figure(3);
loglog(recordZ, angle(record2));
title('Log Phase of Spectrum');
grid on;


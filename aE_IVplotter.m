%%
close all
clear all
setup='2P3DAO';
file='1310302rm';
gsc='g1_s1_c3';
sweeps=[1,15];
filter=[5000,15000];
locations=marcicucca_locations;
load([locations.tgtardir,'MATLABdata/IV/',setup,'/',file]);
iv=iv.(gsc);
figure(1)
clf
hold on;
for i=1:length(sweeps)
    si=mode(diff(iv.time));
    
    v=iv.(['v',num2str(sweeps(i))]);
    if  filter(min(length(filter),i))>0
        [b,a]=butter(1,filter(i)/(1/mode(diff(iv.time)))/2,'low');
        v=filtfilt(b,a,v);
    end
    plot(iv.time,v,'k-','LineWidth',2);
    axis tight
    xlim([0 1])
    axis off
end
axis tight
% plot2svg([figuresfolder,'IV/',num2str(cellstoplot(tempi)),'-',num2str(tempii),'.svg'],gcf);
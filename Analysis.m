%% Analyse V0.2.00
%           spike2 Data and Tracked Data (Manu) in 1WS.mat - 10WS.mat
%           Social SetUp Ben for Data from: Aron Isabella Inga Anni
%

%% The MIT License (MIT)
% 
% Copyright (c) 2016 Benjamin H Paffhausen
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%

%% FIGURE 1 A
% plotting the closest bees with neuronal activity as color on tracks relative to recBee rA1 
% angle and position of periferBee RELATIV to egocentric recBee

x_rel=NaN(size(x));
y_rel=NaN(size(x));
angle_rel=NaN(size(x));
for i=1:length(x)
    if isnan(angle(i,1))
    else
        x_rel(i,:)=x(i,:)-x(i,1);
        y_rel(i,:)=y(i,:)-y(i,1);
        angle_rel(i,:)=mod(angle(i,:)-angle(i,1),360);
        R=rotx(angle(i,1));
        for j=2:12
            temp_a=[1;x_rel(i,j);y_rel(i,j)];
            temp_b=R*temp_a;
            y_rel(i,j)=temp_b(3);
            x_rel(i,j)=temp_b(2);
        end
    end
end

combiTest_rel=cat(3,distance_to_main,x_rel,y_rel);
for i=1:length(x)
    temp(:,:)=combiTest_rel(i,2:12,:);
    temp= sortrows(temp); % all other then the recBee sorted
    combiTest_rel(i,2:12,:)=temp;
    clear temp;
end

figure
combiTest_rA1=[ combiTest_rel(:,12,2),combiTest_rel(:,12,3),...
    combiTest_rel(:,11,2),combiTest_rel(:,11,3),combiTest_rel(:,10,2),combiTest_rel(:,10,3),...
    combiTest_rel(:,9,2),combiTest_rel(:,9,3),combiTest_rel(:,8,2),combiTest_rel(:,8,3),...
    combiTest_rel(:,7,2),combiTest_rel(:,7,3),combiTest_rel(:,6,2),combiTest_rel(:,6,3),...
    combiTest_rel(:,5,2),combiTest_rel(:,5,3),combiTest_rel(:,4,2),combiTest_rel(:,4,3),...
    combiTest_rel(:,3,2),combiTest_rel(:,3,3),combiTest_rel(:,2,2),combiTest_rel(:,2,3),rAverage1];
combiTest_rA1=sortrows(combiTest_rA1,23);
scatter(combiTest_rA1(:,21),combiTest_rA1(:,22),[],combiTest_rA1(:,23),'filled')
xlim([-250 250])
ylim([-250 250])
title('Spike rate as color on relative coordinated to focal bee')


%% FIGURE 1 B
% distance to closest vs activity boxplots 
figure
if exist('rAverage2')==0
    rAverage2=rAverage1;
    rAverage3=rAverage1;
    rAverage4=rAverage1;
end    
distance_rAverage1234=[distance_to_main(:,2:12) rAverage1 rAverage2 rAverage3 rAverage4];
teiler=10;
k=1;
distance_bool=nan(length(x),teiler);
for i = 1:length(x)
    distance_rAverage1234(i,1:11)=sortrows(distance_rAverage1234(i,1:11)')';
end
distance_rAverage1234=sortrows(distance_rAverage1234,-(11+k));
maxNOW=max(distance_rAverage1234(:,11+k));
for j=1:teiler
    for i=1:length(x)
        if distance_rAverage1234(i,11+k)>(j-1)*(maxNOW/teiler) && distance_rAverage1234(i,11+k)<(j)*(maxNOW/teiler)
            distance_bool(i,j) =  distance_rAverage1234(i,1);
        end
    end
end
if nansum(distance_bool(:,10))>0 && nansum(distance_bool(:,9))>0 && ...
        nansum(distance_bool(:,8))>0 && nansum(distance_bool(:,7))>0 && ...
        nansum(distance_bool(:,6))>0 && nansum(distance_bool(:,5))>0 && ...
        nansum(distance_bool(:,4))>0 && nansum(distance_bool(:,3))>0 && ...
        nansum(distance_bool(:,2))>0 && nansum(distance_bool(:,1))>0
    [p,h,stats] = ranksum(distance_bool(:,9),distance_bool(:,10));p
    [p,h,stats] = ranksum(distance_bool(:,8),distance_bool(:,9));p
    [p,h,stats] = ranksum(distance_bool(:,7),distance_bool(:,8));p
    [p,h,stats] = ranksum(distance_bool(:,6),distance_bool(:,7));p
    [p,h,stats] = ranksum(distance_bool(:,5),distance_bool(:,6));p
    [p,h,stats] = ranksum(distance_bool(:,4),distance_bool(:,5));p
    [p,h,stats] = ranksum(distance_bool(:,3),distance_bool(:,4));p
    [p,h,stats] = ranksum(distance_bool(:,2),distance_bool(:,3));p
    [p,h,stats] = ranksum(distance_bool(:,1),distance_bool(:,2));p
end

boxplot(distance_bool,'Orientation','horizontal','Labels',{num2str(round(maxNOW/10,1)),num2str(round(2*maxNOW/10,1)),...
    num2str(round(3*maxNOW/10,1)),num2str(round(4*maxNOW/10,1)),num2str(round(5*maxNOW/10,1)),...
    num2str(round(6*maxNOW/10,1)),num2str(round(7*maxNOW/10,1)),...
    num2str(round(8*maxNOW/10,1)),num2str(round(9*maxNOW/10,1)),num2str(round(maxNOW,1))})
text(100,10,num2str(length(distance_bool)))
xlim([0 200])
title('distance of closest PeriferBee (x-axis) vs. spike rate in 10% bins (y-axis)')


%% FIGURE 2
% spike rate variance distribution for different behaviour

walkingspeed=smooth(double(distance(:,1)),10);
breite = 40;
clear pauses
n=1;
for i=100:length(walkingspeed)-100
    if sum(walkingspeed(i-10:i))==0 && walkingspeed(i+1)>0
        pauses(n)=i;
        n=n+1;
    end
end

walkSpike=zeros(breite*2,length(pauses));
o=1;
for i=1:length(pauses)
    walkSpike(:,i)=(rAverage1(pauses(i)-breite:pauses(i)+breite-1));
    walkPSTH(:,o)=distance_block(pauses(i)-breite*2:pauses(i)+breite*2,3)/nanmean(distance_block(pauses(i)-breite*2:pauses(i)+breite*2,3));
    o=o+1;
end
walkSpikeRAND=zeros(breite*2,length(pauses));
randP=1:round(length(rAverage1)/length(pauses)):length(rAverage1);
randP(1)=100;
for i=1:length(pauses)-100
    walkSpikeRAND(:,i)=(rAverage1(randP(i)-breite:randP(i)+breite-1));
end
varRand=var(walkSpikeRAND);
varwalk=var(walkSpike);
varW=[varwalk', varRand'];

% boxplots alone | contact | random
m=1;
max_abstand=70;
for j=1:length(x(:,1))
    if  distance_block(j,1)>max_abstand &&distance_block(j,2)>max_abstand
        m=m+1;
    else
        if m> 2*breite
            distance_block(j-(m-breite),8)=1;
        end
        m=1;
    end
end
g=1;
varianz_alone=zeros(sum(distance_block(:,8)),1);
for j=breite+1:length(x(:,1))-breite+1
    if distance_block(j,8)==1
        varianz_alone(g,1)=var(distance_block(j-breite:j,3));
        aloneSpikes(g,:)=distance_block(j-breite:j+breite,3);
        g=g+1;
    end
end

g=1;
varianz_contact=zeros(sum(distance_block(:,6)),1);
for j=breite+1:length(x(:,1))-breite+1
    if distance_block(j,6)==1
        varianz_contact(g,1)=var(distance_block(j-breite:j,3));
        g=g+1;
    end
end

g=1;
n=1;
hupiflupiTausend=length(x(:,1))/sum(distance_block(:,6));
hupiflupiTausend=int32(hupiflupiTausend);
varianz_rand=zeros(sum(distance_block(:,6)),1);
for j=breite+1:length(x(:,1))-(breite+100)
    n=n+1;
    if n==hupiflupiTausend
        varianz_rand(g,1)=var(distance_block(j-breite:j,3));
        g=g+1;
        n=1;
    end
end

% walking behavior in contact as divider
speed_grenze=.1;
m=1;
p=1;
g=1;
varianz_walkingPRE=nan(sum(distance_block(:,6)),2);
for j=breite*2+1:length(x(:,1))-breite*3
    if distance_block(j,6)==1
        if nanmean(distance_block(j-breite:j,11))*2<nanmean(distance_block(j:j+breite,11))
            varianz_walkingPRE(g,1)=var(distance_block(j-breite:j,3));
            ActivePSTH(:,m)=distance_block(j-breite*2:j+breite*2,3)/nanmean(distance_block(j-breite*2:j+breite*2,3));
            m=m+1;
        else
            varianz_walkingPRE(g,2)=var(distance_block(j-breite:j,3));
            passivePSTH(:,p)=distance_block(j-breite*2:j+breite*2,3)/nanmean(distance_block(j-breite*2:j+breite*2,3));
            p=p+1;
        end
        g=g+1;
    end
end

varianz_walking=varianz_walkingPRE(:,1);
varianz_walking(isnan(varianz_walking))=[];
varianz_NOTwalking=varianz_walkingPRE(:,2);
varianz_NOTwalking(isnan(varianz_NOTwalking))=[];
varBEHAV=nan(900,5);
varBEHAV(1:length(varianz_alone),2)=varianz_alone;
varBEHAV(1:length(varianz_NOTwalking),4)=varianz_NOTwalking; %% new conmtact ohne walking onset (nur passive)
varBEHAV(1:length(varianz_rand),1)=varianz_rand;
varBEHAV(1:length(varwalk),3)=varwalk;
varBEHAV(1:length(varianz_walking),5)=varianz_walking;

figure
boxplot(varBEHAV,'Labels',{['random ', num2str(length(varianz_rand))],...
    ['alone ', num2str(length(varianz_alone))],...
    ['walking onset ',num2str(length(varwalk))],...
    ['passive contact ', num2str(length(varianz_NOTwalking))],...
    ['active contact ', num2str(length(varianz_walking))]...
    })


%% FIGURE 3

figure
ActivePSTH(isnan(ActivePSTH))=0;
plot(sum(ActivePSTH')/(m-1),'b')
hold on
passivePSTH(isnan(passivePSTH))=0;
plot(sum(passivePSTH')/(p-1),'k')
walkPSTH(isnan(walkPSTH))=0;
plot(sum(walkPSTH')/(o-1),'r')

[p,h,stats] = ranksum(varBEHAV(:,2), varBEHAV(:,1));
statistik(3)=p;
[p,h,stats] = ranksum(varBEHAV(:,3), varBEHAV(:,1));
statistik(4)=p;
[p,h,stats] = ranksum(varBEHAV(:,4), varBEHAV(:,1));
statistik(5)=p;
[p,h,stats] = ranksum(varBEHAV(:,5), varBEHAV(:,1));
statistik(6)=p;
[p,h,stats] = ranksum(varBEHAV(:,5), varBEHAV(:,4));
statistik(7)=p;
[p,h,stats] = ranksum(varBEHAV(:,5), varBEHAV(:,3));
statistik(8)=p;
distance(distance>10)=10;
data=[distance(:,1), rAverage1];
dataSort=sortrows(data,2);
[rho, alpha]=corr(double(data),'Type','Pearson');
statistik(1)=alpha(1,2);
statistik(2)=rho(1,2);
disp(statistik')

%% FIGURE 4 A
clear peakabooSpeed peakabooSpeedVar VORpeakabooSpeedVar peakabooDist...
    peakabooDist VORpeakabooDistVar RANDpeakabooSpeed RANDpeakabooSpeedVar...
    RANDpeakabooDist RANDpeakabooDistVar
[pks,locs,widths,proms]=findpeaks(rAverage1,'MinPeakProminence',4,'Annotate','extents');
subs=round(sqrt(length(locs)))+1;
randLOCS=(1:round(length(rAverage1)/length(locs)):length(rAverage1));
breite=100;
moment_peak=zeros(1,length(angle));
for i=10:length(locs)-5
    if locs(i+1)+breite < size(distance_block,1)
        peakabooSpeed(:,i)=distance_block(locs(i+1)-breite:locs(i+1)+breite,11);
        peakabooSpeedVar(i)=var(distance_block(locs(i+1):locs(i+1)+breite,11));
        VORpeakabooSpeedVar(i)=var(distance_block(locs(i+1)-breite:locs(i+1),11));
        GESAMTpeakabooSpeedVar(i)=var(distance_block(locs(i+1)-breite:locs(i+1)+breite,11));
        peakabooDist(:,i)=distance_block(locs(i+1)-breite:locs(i+1)+breite,1);
        peakabooDistVar(i)=var(distance_block(locs(i+1):locs(i+1)+breite,1));
        VORpeakabooDistVar(i)=var(distance_block(locs(i+1)-breite:locs(i+1),1));
        GESAMTpeakabooDistVar(i)=var(distance_block(locs(i+1)-breite:locs(i+1)+breite,1));
        RANDpeakabooSpeed(:,i)=distance_block(randLOCS(i+1)-breite:randLOCS(i+1)+breite,11);
        RANDpeakabooSpeedVar(i)=var(distance_block(randLOCS(i+1):randLOCS(i+1)+breite,11));
        RANDpeakabooDist(:,i)=distance_block(randLOCS(i+1)-breite:randLOCS(i+1)+breite,1);
        RANDpeakabooDistVar(i)=var(distance_block(randLOCS(i+1):randLOCS(i+1)+breite,1));
        moment_peak(locs)=1;
    end
end

nanmean(peakabooSpeed);
bothSpeedVar=[peakabooSpeedVar', RANDpeakabooSpeedVar'];
figure
boxplot(bothSpeedVar,'Labels',{'PSRC','random'})
ranksum(bothSpeedVar(:,1),bothSpeedVar(:,2));
disp(['speed-trend    ',num2str(ans)])
title('Variance of walking speed distribution');

bothDistVar=[peakabooDistVar', RANDpeakabooDistVar'];
figure
boxplot(bothDistVar,'Labels',{'PSRC','random'})
ranksum(bothDistVar(:,1),bothDistVar(:,2));
disp(['dist-trend    ',num2str(ans)])
title('Variance of distance to closest bee distribution');

%% FIGURE 4 B

figure
plot(median((peakabooSpeed(80:120,:)')));
title('median of walking speed');

figure
plot(median((peakabooDist(80:120,:)')));
title('median of distance to closest bee');

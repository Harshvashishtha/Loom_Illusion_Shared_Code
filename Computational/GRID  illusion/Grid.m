
clear all
clc
%%
sizeX=7;
sizeY=7;
xx = [1:sizeX];
yy = [1:sizeY];

[X,Y] = meshgrid(xx,yy);

DC=0;

%% Create offsets

offset = (randperm(numel(X))-1)/(numel(X)-1);
offset = offset + 1E-10*(~offset);
offset = offset - 1E-10*(offset == 1);

%% Define parameters
P=0.5; %period in s
T=P*60; %Convert it to frames
del=1/T; %Offset Increament in each frame
frames=360;
if (DC~=0)|| (DC~=1)
S1=(1)/(DC);
S2=(-1/(1-DC));
In=(1-(DC*S2));
end

if DC==1
    fun = @(x) (x == 0).*(-1)+(0 < x).*(-1+2*x);
elseif DC==0
    fun = @(x) (x == 0).*(1)+(0 < x).*(1-2*x);
end


LinearMap = offset;
patt(1,:)=fun(LinearMap);
premap(:,:,1)=double(reshape(patt(1,:),[sizeY,sizeX]));
m(1,1)=mean(patt(1,:));
for i=2:frames
LinearMap=mod(LinearMap+del,1);
patt(i,:)=fun(LinearMap);
grid=reshape(patt(i,:),[sizeY,sizeX]);
premap(:,:,i)=grid;
m(i,1)=mean(grid(:));
end

%%
% sigma1=3;
% sigma2=1;
% 
% RF1=-(1/(2*pi*sigma1*sigma1))*exp(-((X-sizeX/2).^2+(Y-sizeY/2).^2)/(2*sigma1*sigma1));
% RF2=(1/(2*pi*sigma2*sigma2))*exp(-((X-sizeX/2).^2+(Y-sizeY/2).^2)/(2*sigma2*sigma2));
% 

% % RF=zeros(sizeX,sizeY);
% RFsurround=RF1;
% RFcenter=RF2;
%%

[X,Y] = meshgrid(xx,yy);

rr=0;
R=2;
RFsurround=double(((X-sizeX/2).^2 + ((Y-sizeY/2).^2) <= R*R)-((X-sizeX/2).^2 + ((Y-sizeY/2).^2) <= rr*rr));
RFcenter=double((X-sizeX/2).^2 + ((Y-sizeY/2).^2) <= rr*rr);


RF0=RFcenter+(-1*RFsurround);
RF=RF0;
RF(RF==0)=nan;

RFcenter(RFcenter==0)=nan;
RFsurround(RFsurround~=1)=nan;
RFsurround(RFsurround==1)=-2;



%%

for  i=1:frames

sC=RFcenter.*premap(:,:,i);
sS=RFsurround.*premap(:,:,i);
SRF=RF.*premap(:,:,i);
SigC(i)=nanmean(sC(:));
SigS(i)=nanmean(sS(:));
SigRF(i)=nanmean(SRF(:));

end


signal=SigC+SigS;
% signal(signal>0)=0;

%%
figure
plot(SigRF)
hold on
yline(0)

%%
clims = [-1 1];
figure
subplot(1,3,1)
imagesc(premap(:,:,1),clims)
colormap gray
subplot(1,3,2)
imagesc(RF0,clims)
colormap gray
subplot(1,3,3)
imagesc(RF0.*premap(:,:,1),clims)
colormap gray

%%
figure

for  i=1:frames
colormap gray
imagesc(premap(:,:,i))
drawnow


end


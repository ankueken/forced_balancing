clear
filesC = dir('Results_lethality/*cancer_balanced.mat');
filesH = dir('Results_lethality/*normal_balanced.mat');

for f=1:length(filesH)

    C=load(['Results_lethality/' filesC(f).name]);
    H=load(['Results_lethality/' filesH(f).name]);
    
    disp(filesC(f).name)
    disp(filesH(f).name)
    
    C.complexes = cell(length(C.model.complexes),1);
        for i=1:length(C.model.complexes)
            C.complexes{i,1} = strjoin(strcat(num2str(C.model.Y(C.model.Y(:,i)~=0,i)), '*', C.model.mets(C.model.Y(:,i)~=0) ),'+');
        end
    
    H.complexes = cell(length(H.model.complexes),1);
        for i=1:length(H.model.complexes)
            H.complexes{i,1} = strjoin(strcat(num2str(H.model.Y(H.model.Y(:,i)~=0,i)), '*', H.model.mets(H.model.Y(:,i)~=0) ),'+');
        end
    
    [joint_complexes, inx_joint_C, inx_joint_H] = intersect(C.complexes,H.complexes);
    length(joint_complexes)
    
    
    lethal_candidates = find(C.Bio_after_balancing(inx_joint_C)==0 & H.Bio_after_balancing(inx_joint_H)>0.9*H.Bio_opt);
        
    CANDIDATE_COMPLEXES{f,1} = C.complexes(inx_joint_C(lethal_candidates));
    [~,candidate_rxns{f,1}] = find(C.model.A(inx_joint_C(lethal_candidates),:)<0);
    CANDIDATE_EC{f,1} = C.model.rxnECNumbers(candidate_rxns{f,1});
    CANDIDATE_SYSTEM{f,1} = C.model.subSystems(candidate_rxns{f,1});

end

% which candidate complexes are joint across tissue types
All_complexes_joint = cell(0,0);
for f=1:9
        All_complexes_joint = union(All_complexes_joint,CANDIDATE_COMPLEXES{f});
end

Data = zeros(length(All_complexes_joint),9);
for f=1:9
    [~,inx_all] = intersect(All_complexes_joint,CANDIDATE_COMPLEXES{f});
    Data(inx_all,f) = 1;
end

% upSetMDemo
% @author : slandarer
% Zhaoxu Liu / slandarer (2023). upset plot 
% (https://www.mathworks.com/matlabcentral/fileexchange/123695-upset-plot), 
% MATLAB Central File Exchange. 检索来源 2023/1/22.
setName={'breast' 'gastric' 'kidney' 'blood' 'liver' 'lung (SCC)' 'lung (AC)' 'ovarian' 'pancreas'};
% setName={'A','B','C','D','E'};
% Data=rand([200,7])>.9;

% bar1Color=[61,58,61]./255;
bar2Color=[161,158,161]./255;
% lineColor=[61,58,61]./255;
bar1Color=[195,66,66]./255;% vertical bars
lineColor=[61,58,61]./255;
hbarColor=[158,158,158]./255;
% bar1Color=[0,0,245;245,0,0]./255;
% bar2Color=cool;
% lineColor=[61,58,61]./255;
%% =========================================================================

pBool=abs(dec2bin((1:(2^size(Data,2)-1))'))-48; % combination of classes
% pBool(sum(pBool,2)==1,:)=[];
[pPos,~]=find(((pBool*(1-Data'))|((1-pBool)*Data'))==0); 
sPPos=sort(pPos);dPPos=find([diff(sPPos);1]);
pType=sPPos(dPPos);pCount=diff([0;dPPos]);
[pCount,pInd]=sort(pCount,'descend'); % count intersecting candidates

%% added to cut the plot at certain threshold of intersection
% small_intersection=find(pCount<5);
% pInd(small_intersection)=[];
% pCount(small_intersection)=[];

pType=pType(pInd);
sCount=sum(Data,1); % count of candidates per group
[sCount,sInd]=sort(sCount,'descend');
sType=1:size(Data,2);
sType=sType(sInd);
%% ========================================================================
% figure axes
fig=figure('Units','normalized','Position',[.3,.2,.5,.63],'Color',[1,1,1]);
axI=axes('Parent',fig);hold on;
set(axI,'Position',[.33,.35,.655,.61],'LineWidth',1.2,'Box','off','TickDir','out',...
    'FontName','Helvetica','FontSize',12,'XTick',[],'XLim',[0,length(pType)+1])
axI.YLabel.String='Number of shared candidate complexes';
axI.YLabel.FontSize=14;
%
axS=axes('Parent',fig);hold on;
set(axS,'Position',[.03,.08,.245,.26],'LineWidth',1.2,'Box','off','TickDir','out',...
    'FontName','Helvetica','FontSize',12,'YColor','none','YLim',[.5,size(Data,2)+.5],...
    'YAxisLocation','right','XDir','reverse','YTick',[])
axS.XLabel.String='Number of candidate complexes';
axS.XLabel.FontSize=14;
%
axL=axes('Parent',fig);hold on;
set(axL,'Position',[.33,.08,.655,.26],'YColor','none','YLim',[.5,size(Data,2)+.5],'XColor','none','XLim',axI.XLim)
%% ========================================================================
%  -----------------------------------------------------------
% colorMapLength = length(unique(pCount))-8;
% pink = [153, 51, 0]./255;
% red = [255, 145, 0]./255;
% vbar2Color = [linspace(red(1),pink(1),colorMapLength)', linspace(red(2),pink(2),colorMapLength)', linspace(red(3),pink(3),colorMapLength)'];
vbar2Color = autumn(50);
vbar2Color = [repmat([1 0 0],7,1); vbar2Color];
[upCount,inxupCount] = unique(pCount);
barHdlI=bar(axI,pCount,'FaceColor','flat');
barHdlI.EdgeColor='none';

for i=1:length(unique(pCount))-8
    barHdlI.CData(pCount==upCount(i),:) = repmat(vbar2Color(size(vbar2Color,1)-i,:),sum(pCount==upCount(i)),1);
end
barHdlI.CData([1:8,19],:) = repmat([1 0 0],9,1);

% if size(bar1Color,1)==1
%     bar1Color=[bar1Color;bar1Color];
% end
% tx=linspace(0,1,size(bar1Color,1))';
% ty1=bar1Color(:,1);ty2=bar1Color(:,2);ty3=bar1Color(:,3);
% tX=linspace(0,1,length(pType))';
% bar1Color=[interp1(tx,ty1,tX,'pchip'),interp1(tx,ty2,tX,'pchip'),interp1(tx,ty3,tX,'pchip')];
% barHdlI.FaceColor='flat';
% for i=1:length(pType)
%     barHdlI.CData(i,:)=bar1Color(i,:);
% end
text(axI,1:length(pCount),pCount,string(pCount),'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontName','Helvetica','FontSize',10,'Color',[61,58,61]./255)
% 集合统计图 ---------------------------------------------------------------
barHdlS=barh(axS,sCount,'BarWidth',.6);
barHdlS.EdgeColor='none';
barHdlS.BaseLine.Color='none';
for i=1:size(Data,2)
    annotation('textbox',[(axS.Position(1)+axS.Position(3)+axI.Position(1))/2-.02,...
        axS.Position(2)+axS.Position(4)./size(Data,2).*(i-.5)-.02,.04,.04],...
        'String',setName{sInd(i)},'HorizontalAlignment','center','VerticalAlignment','middle',...
        'FitBoxToText','on','LineStyle','none','FontName','Helvetica','FontSize',10)
end
if size(bar2Color,1)==1
    bar2Color=[bar2Color;bar2Color];
end
tx=linspace(0,1,size(bar2Color,1))';
ty1=bar2Color(:,1);ty2=bar2Color(:,2);ty3=bar2Color(:,3);
tX=linspace(0,1,size(Data,2))';
bar2Color=[interp1(tx,ty1,tX,'pchip'),interp1(tx,ty2,tX,'pchip'),interp1(tx,ty3,tX,'pchip')];
barHdlS.FaceColor='flat';
sstr{size(Data,2)}='';
for i=1:size(Data,2)
    barHdlS.CData(i,:)=bar2Color(i,:);
    sstr{i}=[num2str(sCount(i)),' '];
end
text(axS,sCount,1:size(Data,2),sstr,'HorizontalAlignment','right',...
    'VerticalAlignment','middle','FontName','Helvetica','FontSize',10,'Color',[61,58,61]./255)
% 绘制关系连线 ---------------------------------------------------------------
patchColor=[248,246,249;255,254,255]./255;
for i=1:size(Data,2)
    fill(axL,axI.XLim([1,2,2,1]),[-.5,-.5,.5,.5]+i,patchColor(mod(i+1,2)+1,:),'EdgeColor','none')
end
[tX,tY]=meshgrid(1:length(pType),1:size(Data,2));
plot(axL,tX(:),tY(:),'o','Color',[233,233,233]./255,...
    'MarkerFaceColor',[233,233,233]./255,'MarkerSize',8);
for i=1:length(pType)
    tY=find(pBool(pType(i),:));
    oY=zeros(size(tY));
    for j=1:length(tY)
        oY(j)=find(sType==tY(j));
    end
    tX=i.*ones(size(tY));
    plot(axL,tX(:),oY(:),'-o','Color',lineColor(1,:),'MarkerEdgeColor','none',...
        'MarkerFaceColor',lineColor(1,:),'MarkerSize',8,'LineWidth',2);
end

% Zhaoxu Liu / slandarer (2023). upset plot 
% (https://www.mathworks.com/matlabcentral/fileexchange/123695-upset-plot), 
% MATLAB Central File Exchange. 检索来源 2023/1/22.

function results = msd_analysis(tracks, shortestTrack, pointsToFit, maxOffset, msdOrLinear)

%
%Analyze the mean squared displacements (msd) of a set of single-molecule tracks
%by fitting either a power-law or a linear function. When msd are fitted
%with a power-law, it is additionally fitted with a confined motion model i
%to calculate the confinement radius. See TrackIt manual for more
%information.
%
%results = msd_analysis(tracks, shortestTrack, pointsToFit, maxOffset, msdOrLinear)
%
%Input:
%   tracks          -   Cell array containing one track per cell. Each cell
%                       contains an array with at least 3 columns: [frame,xpos,ypos].                      
%   shortestTrack   -   Tracks with less detections than shortestTrack are
%                       ignored in the msd analysis
%   pointsToFit     -   The first n points of the msd that should be fitted.
%                       Must be a string containing either a numeric value 
%                       or percentage (eg. '10' or '90%': means fitting is 
%                       performed on the first 10 points or on 90% of the
%                       points of the msd, respectively.
%   maxOffset       -   Maximum allowed offset value for the fitting
%                       function
%   msdOrLinear     -   Use 'MSD' for fitting a power-law or 'Linear' for
%                       fitting a linear function
%
%Output:
%   results          -  Struct with following fields:
%                           msdDiffConst: diffusion constants
%                           alphaValues: alpha values
%                           confRad: confinement radii (if fit with power-law)
%                           meanJumpDist: average jump distance
%                           msd: mean squared displacement
%                           offset: offset
%                           nPointsFitted: amount of points that have been fitted
%                       
                              

%Default parameters
para=struct('shortest_track',shortestTrack,...
    'max_offset', maxOffset,...
    'points_to_fit', pointsToFit,...
    'msd_or_linear', msdOrLinear);

%Initialize
nTracks = length(tracks);
msdDiffConst = NaN(nTracks,1);
alphaValues = NaN(nTracks,1);
confRad = NaN(nTracks,1);
meanJumpDist = NaN(nTracks,1);
msd = cell(nTracks,1);
nPointsFitted = NaN(nTracks,1);
offset = NaN(nTracks,1);

%Correct for dark frames
[tr_temp,~]=correctdark(tracks);

%Iterate through tracks
if nTracks == 1
    [res_t, curMsd] = alphaana(tr_temp{1},para);
    
    if ~isempty(res_t)
        %Fit was succesfull so save results in structure array
        msdDiffConst(1) = res_t.msdDiffConst;
        alphaValues(1) = res_t.alphaValues;
        confRad(1) = res_t.confRad;
        meanJumpDist(1) = res_t.meanJumpDist;
        offset(1) = res_t.offset;
        nPointsFitted(1) = res_t.nPointsFitted;
        msd{1} = curMsd;
    end
else
    parfor trackIdx = 1:nTracks
        [res_t, curMsd] = alphaana(tr_temp{trackIdx},para);
        
        if ~isempty(res_t)
            %Fit was succesfull so save results in structure array
            msdDiffConst(trackIdx) = res_t.msdDiffConst;
            alphaValues(trackIdx) = res_t.alphaValues;
            confRad(trackIdx) = res_t.confRad;
            meanJumpDist(trackIdx) = res_t.meanJumpDist;
            offset(trackIdx) = res_t.offset;
            nPointsFitted(trackIdx) = res_t.nPointsFitted;
            msd{trackIdx} = curMsd;
        end
    end
end

results.msdDiffConst = msdDiffConst;
results.alphaValues = alphaValues;
results.confRad = confRad;
results.meanJumpDist = meanJumpDist;
results.msd = msd;
results.offset = offset;
results.nPointsFitted = nPointsFitted;
end

%--------------------------------------------------------------------------
function [result,msd]=alphaana(track,pa)


pointsToFit = pa.points_to_fit; %Amount of points in the msd to fit (either numeric or percentage)

%Calculate msd
msd = calcmsd(track);
trackLength = size(track,1);

if contains(pointsToFit,'%')
    %pointsToFit is a percentage value so calcuelate the corresponding
    %amount of points in the msd
    percentage = str2double(regexp(pointsToFit,'\d*','Match'));
    pointsToFit = round(numel(msd)*percentage/100);
else    
    %pointsToFit is not percentage value so conert to numeric value
    pointsToFit = str2double(pointsToFit);
end

%Make sure pointsToFit is not higher than the number of elements in msd
pointsToFit = min(numel(msd), pointsToFit);

%Cut off msd at the number of points to fit
msdToFit = msd(1:pointsToFit);

%Create time vector
t = (1:numel(msdToFit))';

%Get maximum allowed offset value
maxOffset = pa.max_offset;

%Fit msd
para=[NaN,NaN,NaN];
options=optimoptions('lsqnonlin','Display','off');

if trackLength>=pa.shortest_track
    if strcmp(pa.msd_or_linear,'Linear')   
        %Fit msd using a linear fit
        msd_power=@(para,x,y)(y-(4*para(1)*x.^para(2)+para(3)));
        para=lsqnonlin(@(para)msd_power(para,t,msdToFit),[1,1,1],[0,1,0],[inf,1,maxOffset],options);
    else
        %Fit msd using a powerlaw
        msd_power=@(para,x,y)(y-(4*para(1)*x.^para(2)+para(3)));
        para=lsqnonlin(@(para)msd_power(para,t,msdToFit),[1,1,1],[0,0,0],[inf,inf,maxOffset],options);
    end        
end

%Check if fit was succesful
if not(isnan(para(1)))     %&& para(2) > 0.001
    if strcmp(pa.msd_or_linear,'MSD')
        %Fit Radius of confinement
        para1=lsqnonlin(@(para)conf_obj(para,t,msdToFit),[1,1],[0,0],[],options);
    else
        para1 = NaN;
    end
    %y1=-conf_obj(para1,t,msd*0);
    %distances=jumpdistance(track);
    %if sum(distances>(msd+1.5))>0
%       disp('butterfly!')
    %end
    result.msdDiffConst = para(1);
    result.alphaValues = para(2);
    result.offset = para(3);       
    result.nPointsFitted = pointsToFit;
else
    result=[];
    return
end

%Save Results
if ~any(isnan(para1)) && para1(1) < 20
    result.confRad = para1(1);    
    result.meanJumpDist = meanJumpDistance(track); 
else
    result.confRad = NaN;
    result.meanJumpDist = NaN;
end

end

%--------------------------------------------------------------------------
function msd = calcmsd(positions)
%Calculate mean squared displacement

trackLength = size(positions,1);

msd=zeros(trackLength-1,1);

for deltaT=1:trackLength-1
    nSegments = trackLength-deltaT;
    for w=1:nSegments
        msd(deltaT)=msd(deltaT)+(positions(w,1)-positions(w+deltaT,1))^2+(positions(w,2)-positions(w+deltaT,2))^2;
    end
    msd(deltaT)=msd(deltaT)/nSegments;
end

end

%--------------------------------------------------------------------------
function [d]=conf_obj(para,x,y)
%Fit function for confinement radius fit
R=para(1);
sigma=para(2);
offset=0;

d=y-(R^2*(1-exp(-sigma*x/R^2))+offset);

end

%--------------------------------------------------------------------------
function [tracks,addedspots]=correctdark(tracks)
%Correct for bridged dark frames
    addedspots=0;

    for j=1:length(tracks)
            if not(isempty((tracks{j})))
                timescale=tracks{j}(:,1);
                timescale=timescale-timescale(1)+1;
                X = abs(tracks{j}(:,2));
                Y = abs(tracks{j}(:,3));
                                
                initiallength=numel(X);

                %Search for dark frames
                idx=1;
                while idx<=numel(X)
                   if not(idx==timescale(idx))
                       jumpto=timescale(idx)-1;
                       unitvector=[X(idx)-X(idx-1),Y(idx)-Y(idx-1)]/numel((1:(jumpto-idx+2)));
                       %fillX=X(idx-1)+(X(idx)-X(idx-1))*(1:(jumpto+1-idx))/numel((1:(jumpto+1-idx)));
                       %fillY=Y(idx-1)+(Y(idx)-Y(idx-1))*(1:(jumpto+1-idx))/numel((1:(jumpto+1-idx)));
                       fillX=X(idx-1)+(1:(jumpto-idx+1))*unitvector(1);
                       fillY=Y(idx-1)+(1:(jumpto-idx+1))*unitvector(2);
                       X=[X(1:idx-1);fillX';X(idx:end)];
                       Y=[Y(1:idx-1);fillY';Y(idx:end)];
                       timescale=[timescale(1:idx-1);(idx:jumpto)';timescale(idx:end)];                                   
                   end
                   idx=idx+1;
                end                    
            end
            %Save into variable
            tracks{j}=zeros(numel(X),2);
            tracks{j}(:,[1,2])=[X,Y];   
            addedspots=addedspots+numel(X)-initiallength;
    end
end

function meanJumpDist = meanJumpDistance(track)
%Calculate mean jump distance

%Extact variables
x=track(:,1);
y=track(:,2);

%Calculate distances between two points
distances=zeros(1,numel(x));
for q=1:numel(x)-1
    distances(q)=sqrt((x(q)-x(q+1))^2+(y(q)-y(q+1))^2);
end

meanJumpDist = mean(distances);

end
%% Preapare weather dataset (lon x lat x year x  day of year)
clear;
clc;

% download data from: https://drive.google.com/file/d/1-FTJFXhTOpQQhfgTgq6s3mMHo4nKYI1O/view?usp=sharing

data=ncread('df_2004_2021.nc', 'TMAX');
allDates = datetime(2004, 11, 16): datetime(2020, 11, 15);

startDate = datetime(2005, 01, 01);
endDate = datetime(2019, 12, 31);

sel_idx = (allDates >= startDate) & (allDates<=endDate);
selDates = allDates(sel_idx);
selData = data(:,:, sel_idx);

lons = linspace(4.0, 50.5, 94)';
lats = linspace(30.0, 54.5, 50)';

[reshaped, years, daysCount] = yearlyReshapeFunc(selData, startDate, endDate);

%% Fill missing entries with moving avg.
mask = isnan(reshaped);
dim = 4; % dim. along which to compute moving AVEG. 
ma = movmean(reshaped, [5 0], dim, 'omitnan');
filledData = reshaped;
filledData(mask) = ma(mask);

% If there are still NaNs, issue a warning
if any(isnan(filledData(:)))
    warning('Some NaN values could not be filled.');
end

% add small noise
filledData  = noisy(filledData, 99);

%%=============================
%% SAVE YEARLY DATA TO NETCDF
%%=============================

outFile = 'YearlyDataf.nc';
[nLon,nLat] = deal(length(lons), length(lats));
nYears = length(years);
maxDays = 366;  % fixed day-of-year dimension

nccreate(outFile,'lon','Dimensions',{'lon',nLon});
nccreate(outFile,'lat','Dimensions',{'lat',nLat});
nccreate(outFile,'year','Dimensions',{'year',nYears});
nccreate(outFile,'dayOfYear','Dimensions',{'dayOfYear',maxDays});

nccreate(outFile,'TMAX', ...
    'Dimensions',{'lon',nLon,'lat',nLat,'year',nYears,'dayOfYear',maxDays}, ...
    'Datatype','double');

ncwrite(outFile,'lon',lons);
ncwrite(outFile,'lat',lats);
ncwrite(outFile,'year',years);
ncwrite(outFile,'dayOfYear',1:maxDays);
ncwrite(outFile,'TMAX',filledData);

ncwriteatt(outFile,'lon','units','degrees_east');
ncwriteatt(outFile,'lat','units','degrees_north');
ncwriteatt(outFile,'year','description','Calendar year');
ncwriteatt(outFile,'dayOfYear','description','Day of year (1=Jan1,...,366)');
ncwriteatt(outFile,'TMAX','description','Yearly reshaped data array');

function [reshaped, years, daysCount] = yearlyReshapeFunc(data, startDate, endDate)
%YEARLYRESHAPEFUNC   Time-aware annual reshape
%   [RESHAPED, YEARS, DAYSCOUNT] = yearlyReshapeFunc(DATA, STARTDATE, ENDDATE)
    timeVec = startDate : endDate;
    [nLon,nLat,nTime] = size(data);
    assert(nTime==numel(timeVec),'Data time dim must match date span');

    years = unique(year(timeVec));
    nYears = numel(years);

    % Days per year (for reference)
    daysCount = arrayfun(@(y) sum(year(timeVec)==y), years);

    % Preallocate for full Julian year (max 366 days)
    maxDays = 366;
    reshaped = nan(nLon, nLat, nYears, maxDays);

    % Fill by day-of-year index
    for t = 1:nTime
        Y = year(timeVec(t));
        doy = day(timeVec(t), 'dayofyear');
        yi = find(years==Y, 1);
        reshaped(:,:,yi,doy) = data(:,:,t);
    end
end




function reconstructed = yearlyReconstruct(reshaped, years, timeVec)
%YEARLYRECONSTRUCT   Restore original 3D from annual array
    [nLon,nLat,nYears,maxDays] = size(reshaped);
    nTime = numel(timeVec);
    reconstructed = nan(nLon,nLat,nTime);

    for t = 1:nTime
        Y = year(timeVec(t));
        doy = day(timeVec(t), 'dayofyear');
        yi = find(years==Y, 1);
        reconstructed(:,:,t) = reshaped(:,:,yi,doy);
    end
end


%% Preapare weather dataset (lon x lat x year xseason x day of season)
clear;
clc;
% download data from: https://drive.google.com/file/d/1-FTJFXhTOpQQhfgTgq6s3mMHo4nKYI1O/view?usp=sharing
data=ncread('df_2004_2021.nc', 'TMAX');
startDate = datetime(2004, 11, 16);
endDate = datetime(2020, 11, 15);
lons = linspace(4.0, 50.5, 94)';
lats = linspace(30.0, 54.5, 50)';

[reshaped, years, defSeasons, daysCount] = seasonalReshapefunc(data, startDate, endDate);

clear data

%% Fill missing entries with kNN
filledData = reshaped;
dims = ndims(filledData);

% Iterate over each dimension
for dim = 1:dims
    filledData = fillmissing(filledData, 'nearest', dim);
end

% After one pass, some NaNs might still remain if they were isolated.
% Repeat the process until all NaNs are filled or a maximum number of iterations is reached.
maxIterations = 5;
iteration = 1;
while any(isnan(filledData(:))) && iteration <= maxIterations
    for dim = 1:dims
        filledData = fillmissing(filledData, 'nearest', dim);
    end
    iteration = iteration + 1;
end

% If there are still NaNs, issue a warning
if any(isnan(filledData(:)))
    warning('Some NaN values could not be filled after %d iterations.', maxIterations);
end


%% save data
% Define file name
outFile = 'seasonalDataOrd5.nc';
timeVec = startDate:endDate;
seasonNames = {'DJF', 'MAM', 'JJA', 'SON'};

nSeasons = numel(defSeasons);
nLon = length(lons);
nLat = length(lats);
nYears = length(years);
maxDays = max(daysCount(:));

% Create dimensions & variables
nccreate(outFile, 'lon',    'Dimensions', {'lon', nLon});
nccreate(outFile, 'lat',    'Dimensions', {'lat', nLat});
nccreate(outFile, 'year',   'Dimensions', {'year', nYears});
nccreate(outFile, 'season', 'Dimensions', {'season', nSeasons});
nccreate(outFile, 'day',    'Dimensions', {'day', maxDays});
nccreate(outFile, 'TMAX', ...
    'Dimensions', {'lon', nLon, 'lat', nLat, 'year', nYears, ...
                   'season', nSeasons, 'day', maxDays}, ...
    'Datatype', 'double');

%  Write the coordinate & data arrays
ncwrite(outFile, 'lon',    lons);
ncwrite(outFile, 'lat',    lats);
ncwrite(outFile, 'year',   years);
ncwrite(outFile, 'season', 1:nSeasons);
ncwrite(outFile, 'day',    1:maxDays);
ncwrite(outFile, 'TMAX',   filledData);

%  Add metadata
ncwriteatt(outFile, 'lon',    'units', 'degrees_east');
ncwriteatt(outFile, 'lat',    'units', 'degrees_north');
ncwriteatt(outFile, 'year',   'description', 'Calendar year');
ncwriteatt(outFile, 'season', 'description', '1=DJF,2=MAM,3=JJA,4=SON');
ncwriteatt(outFile, 'day',    'description', 'Day index within season');
ncwriteatt(outFile, 'TMAX',   'description', 'Seasonally reshaped TMAX data');


%% Functions
function [reshaped, years, defSeasons, daysCount] = seasonalReshapefunc(data, startDate, endDate)
% SEASONALRESHAPE   Time-aware seasonal reshape 
%
%   [RESHAPED, YEARS, DEFSEASONS, DAYSCOUNT] = seasonalReshape(DATA, STARTDATE, ENDDATE, LONS, LATS, OUTFILE)
%
%   Inputs:
%     DATA      - [nLon x nLat x nTime] spatiotemporal array
%     STARTDATE - datetime scalar (e.g., datetime(2004,11,16))
%     ENDDATE   - datetime scalar (e.g., datetime(2020,11,15))
%   Outputs:
%     RESHAPED   - [nLon x nLat x nYears x 4 x maxDays] seasonal data
%     YEARS      - vector of full calendar years used
%     DEFSEASONS - cell array of month vectors for each season
%     DAYSCOUNT  - [nYears x 4] number of days in each season per year

    %% 1. Build time vector and detect full years
    timeVec = startDate : endDate;
    nTime   = numel(timeVec);
    [nLon, nLat, tLen] = size(data);
    assert(tLen == nTime, 'Third dimension of data must match date range length');

    % Identify full calendar years
    allYears = unique(year(timeVec));
    years = allYears(allYears > year(startDate) & allYears < year(endDate));
    nYears = numel(years);

    %% 2. Define seasons
    defSeasons = {
        [12, 1, 2];  % DJF
        [3, 4, 5];   % MAM
        [6, 7, 8];   % JJA
        [9,10,11];   % SON
    };
    nSeasons = numel(defSeasons);
    seasonNames = {'DJF', 'MAM', 'JJA', 'SON'};

    %% 3. Count days per season & determine max
    daysCount = zeros(nYears, nSeasons);
    for yi = 1:nYears
        for si = 1:nSeasons
            mask = year(timeVec) == years(yi) & ismember(month(timeVec), defSeasons{si});
            daysCount(yi, si) = sum(mask);
        end
    end
    maxDays = max(daysCount(:));

    %% 4. Preallocate and fill reshaped array
    reshaped = nan(nLon, nLat, nYears, nSeasons, maxDays);
    for yi = 1:nYears
        for si = 1:nSeasons
            sel = find(year(timeVec) == years(yi) & ismember(month(timeVec), defSeasons{si}));
            for di = 1:numel(sel)
                reshaped(:, :, yi, si, di) = data(:, :, sel(di));
            end
        end
    end
end

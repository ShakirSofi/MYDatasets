%% Preapare weather dataset (lon x lat x year x month x day of month)
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

[reshaped, years, months, daysCount] = monthlyReshapeFunc(selData, startDate, endDate);

%clear data

%% Fill missing entries with moving avg.
mask = isnan(reshaped);
dim = 5; % dim. along which to compute moving AVEG. 
ma = movmean(reshaped, [5 0], dim, 'omitnan');
filledData = reshaped;
filledData(mask) = ma(mask);

% If there are still NaNs, issue a warning
if any(isnan(filledData(:)))
    warning('Some NaN values could not be filled.');
end

% add small noise
filledData  = noisy(filledData, 99);

%%{
%% save data
outFile = 'MonthlyData.nc';

nLon    = length(lons);
nLat    = length(lats);
nYears  = length(years);
nMonths = numel(months);
maxDays = max(daysCount(:));

% Create dimensions
nccreate(outFile, 'lon',   'Dimensions', {'lon', nLon});
nccreate(outFile, 'lat',   'Dimensions', {'lat', nLat});
nccreate(outFile, 'year',  'Dimensions', {'year', nYears});
nccreate(outFile, 'month', 'Dimensions', {'month', nMonths});
nccreate(outFile, 'day',   'Dimensions', {'day', maxDays});

% Create data variable
nccreate(outFile, 'TMAX', ...
    'Dimensions', {'lon', nLon, 'lat', nLat, 'year', nYears, 'month', nMonths, 'day', maxDays}, ...
    'Datatype', 'double');

% Write coordinate and data arrays
ncwrite(outFile, 'lon',   lons);
ncwrite(outFile, 'lat',   lats);
ncwrite(outFile, 'year',  years);
ncwrite(outFile, 'month', months);
ncwrite(outFile, 'day',   1:maxDays);
ncwrite(outFile, 'TMAX',  filledData);

% Add metadata attributes
ncwriteatt(outFile, 'lon',   'units',       'degrees_east');
ncwriteatt(outFile, 'lat',   'units',       'degrees_north');
ncwriteatt(outFile, 'year',  'description', 'Calendar year');
ncwriteatt(outFile, 'month', 'description', 'Calendar month (1=Jan,12=Dec)');
ncwriteatt(outFile, 'day',   'description', 'Day index within month');
ncwriteatt(outFile, 'TMAX',  'description', 'Monthly reshaped data array');

%}
%% Functions
function [reshaped, years, months, daysCount] = monthlyReshapeFunc(data, startDate, endDate)
%MONTHLYRESHAPEFUNC   Time-aware monthly reshape
%
%   [RESHAPED, YEARS, MONTHS, DAYSCOUNT] = monthlyReshapefunc(DATA, STARTDATE, ENDDATE)
%
%   Inputs:
%     DATA      - [nLon x nLat x nTime] spatiotemporal array
%     STARTDATE - datetime scalar (e.g., datetime(2004,11,16))
%     ENDDATE   - datetime scalar (e.g., datetime(2020,11,15))
%   Outputs:
%     RESHAPED   - [nLon x nLat x nYears x 12 x maxDays] monthly data
%     YEARS      - vector of full calendar years used
%     MONTHS     - vector 1:12
%     DAYSCOUNT  - [nYears x 12] number of days in each month per year

    %% 1. Build time vector and validate dimensions
    timeVec = startDate : endDate;
    nTime   = numel(timeVec);
    [nLon, nLat, tLen] = size(data);
    assert(tLen == nTime, 'Third dimension of data must match date range length');

    %% 2. Identify full calendar years
    allYears = unique(year(timeVec));
    years = allYears(allYears >= year(startDate) & allYears <= year(endDate));
    nYears = numel(years);

    %% 3. Define months
    months = 1:12;
    nMonths = numel(months);

    %% 4. Count days per month & determine max
    daysCount = zeros(nYears, nMonths);
    for yi = 1:nYears
        for mi = 1:nMonths
            mask = year(timeVec) == years(yi) & month(timeVec) == months(mi);
            daysCount(yi, mi) = sum(mask);
        end
    end
    maxDays = max(daysCount(:));

    %% 5. Preallocate and fill reshaped array
    reshaped = nan(nLon, nLat, nYears, nMonths, maxDays);
    for yi = 1:nYears
        for mi = 1:nMonths
            sel = find(year(timeVec) == years(yi) & month(timeVec) == months(mi));
            for di = 1:numel(sel)
                reshaped(:, :, yi, mi, di) = data(:, :, sel(di));
            end
        end
    end
end
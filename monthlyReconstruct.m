function reconstructed = monthlyReconstruct(reshaped, years, months, timeVec)
%MONTHLYRECONSTRUCT   Rebuild original time series from monthly reshape
%
%   reconstructed = monthlyReconstruct(RESHAPED, YEARS, MONTHS, TIMEVEC)
%
%   Inputs:
%     RESHAPED  - [nLon x nLat x nYears x 12 x maxDays] monthly array
%     YEARS     - vector of calendar years used
%     MONTHS    - vector 1:12 of month indices
%     TIMEVEC   - datetime vector used originally (startDate:endDate)
%
%   Output:
%     reconstructed - [nLon x nLat x numel(timeVec)] original series

    [nLon, nLat, nYears, nMonths, maxDays] = size(reshaped);
    nTime = numel(timeVec);
    reconstructed = nan(nLon, nLat, nTime);

    for yi = 1:nYears
        for mi = 1:nMonths
            sel = find(year(timeVec)==years(yi) & month(timeVec)==months(mi));
            for di = 1:numel(sel)
                reconstructed(:,:,sel(di)) = reshaped(:,:,yi,mi,di);
            end
        end
    end
end
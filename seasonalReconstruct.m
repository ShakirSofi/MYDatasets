% ----------------------
% Reconstruct original 3D time series from seasonal reshape
%
% USAGE:
%   reconstructed = seasonalReconstruct(reshaped, years, defSeasons, timeVec)
%
% INPUT:
%   reshaped    - [nLon x nLat x nYears x 4 x maxDays] seasonal array
%   years       - vector of full calendar years
%   defSeasons  - cell array of month vectors
%   timeVec     - datetime vector used originally (startDate:endDate)
%
% OUTPUT:
%   reconstructed - [nLon x nLat x numel(timeVec)] original series

function reconstructed = seasonalReconstruct(reshaped, years, defSeasons, timeVec)
    [nLon, nLat, nYears, nSeasons, maxDays] = size(reshaped);
    nTime = numel(timeVec);
    reconstructed = nan(nLon, nLat, nTime);

    for yi = 1:nYears
        for si = 1:nSeasons
            sel = find(year(timeVec)==years(yi) & ismember(month(timeVec), defSeasons{si}));
            for di = 1:numel(sel)
                reconstructed(:,:,sel(di)) = reshaped(:,:,yi,si,di);
            end
        end
    end
end

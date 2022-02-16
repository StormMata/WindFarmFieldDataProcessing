function [Indices] = IntFilter(LiCo, TuCo, InCo, TOI, Heights, data)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% LiDAR Conditions
    CondLiAv = zeros(length(Heights),length(data.DateTime));                                    % Preallocate array

    for i = 1:length(Heights)                                                                   % Retrieve Lidar availability
        CondLiAv(i,:) = data.Lidar.(strcat('H',num2str(Heights(i),'%02.f'),'m')).Available == LiCo.LiAv;
    end

% Turbine Conditions

    CondTuAv  = data.(strcat('BHR',num2str(TOI,'%02.f'))).Available         == TuCo.TAv;        % Retrieve TOI availability
    CondOp    = data.(strcat('BHR',num2str(TOI,'%02.f'))).OpSt              == TuCo.TOpSt;      % Retrieve TOI operational status
    CondStrat = data.(strcat('BHR',num2str(TOI,'%02.f'))).WFCStrategy       == TuCo.TWFCStrat;  % Retrieve TOI control strategy
    CondYaw   = data.(strcat('BHR',num2str(TOI,'%02.f'))).YawOffsetApplied  == TuCo.TOffsetAp;  % Retrieve TOI Yaw offset

% Inflow Conditions

    CondWiE   = (data.Lidar.H104m.WndDir >= InCo.W | data.Lidar.H104m.WndDir <= InCo.E);        % Retrieve only northern flow

% Create Filter

    FilterMatrix = CondLiAv .* CondTuAv .* CondOp .* CondStrat .* CondYaw .* CondWiE;           % Combine all individual filters into initial filter
    
    FullIndex = 1:1:size(FilterMatrix,2);                                                       % Generate inital full list of indices
    
    Indices = nonzeros(FullIndex .* FilterMatrix(1,:))';                                        % Indices after filter is applied

end
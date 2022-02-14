function [Indices] = IntFilter(C, T, data)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% LiDAR Conditions

    CondLiAv = zeros(length(T.Heights),length(data.DateTime));                                  % Preallocate array

    for i = 1:length(T.Heights)                                                                 % Retrieve Lidar availability
        CondLiAv(i,:) = data.Lidar.(strcat('H',num2str(T.Heights(i),'%02.f'),'m')).Available == C.LiAv;
    end

    CondLiAv = all(CondLiAv);                                                                   % Identifies only columns where all heights are available

% Turbine Conditions

    CondTuAv  = data.(strcat('BHR',num2str(T.TOI,'%02.f'))).Available         == C.TAv;         % Retrieve TOI availability
    CondOp    = data.(strcat('BHR',num2str(T.TOI,'%02.f'))).OpSt              == C.OpSt;        % Retrieve TOI operational status
    CondStrat = data.(strcat('BHR',num2str(T.TOI,'%02.f'))).WFCStrategy       == C.WFCStrat;    % Retrieve TOI control strategy
    CondYaw   = data.(strcat('BHR',num2str(T.TOI,'%02.f'))).YawOffsetApplied  == C.OffsetAp;    % Retrieve TOI Yaw offset

    if C.NacPos == 1
        CondNac = ~isnan(data.(strcat('BHR',num2str(T.TOI,'%02.f'))).NacellePosition_corrected);% Retrieve Nacelle Heading
    else 
        CondNac = ones(length(CondTuAV));                                                       % Ignore Nacelle Heading
    end

% Inflow Conditions

    CondWiE   = (data.Lidar.H104m.WndDir >= C.W | data.Lidar.H104m.WndDir <= C.E);              % Retrieve only northern flow

FilterVector = CondLiAv .* CondTuAv .* CondOp .* CondStrat .* CondYaw .* CondNac .* CondWiE;    % Combine all individual filters into initial filter

FullIndex = 1:1:size(FilterVector,2);                                                           % Generate inital full list of indices

Indices = nonzeros(FullIndex .* FilterVector)';                                                 % Indices after filter is applied

end
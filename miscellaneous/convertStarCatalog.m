function [sN, mag] = convertStarCatalog(RA_hms, DEC_dms, Vmag)

    % Convert RA: "hh mm ss" → rad
    H = RA_hms(:,1);
    M = RA_hms(:,2);
    S = RA_hms(:,3);
    RA_deg = 15 * (H + M/60 + S/3600);
    alpha = deg2rad(RA_deg);

    % Convert DEC: "dd mm ss" → rad
    D = DEC_dms(:,1);
    M = DEC_dms(:,2);
    S = DEC_dms(:,3);
    DEC_deg = D + M/60 + S/3600;
    delta = deg2rad(DEC_deg);

    % Convert to inertial unit vectors
    sN = [cos(delta).*cos(alpha), ...
          cos(delta).*sin(alpha), ...
          sin(delta)];

    % Magnitude
    mag = Vmag;

end

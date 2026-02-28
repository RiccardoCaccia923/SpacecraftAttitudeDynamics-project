function [starVecs_N , MAG] = starsFromCatalogue(filenametsv,maxMAG)

    T = readtable('starCatalogClean.tsv', 'FileType','text', 'Delimiter','\t');
    T.Properties.VariableNames(1) = "HR";
    T.Properties.VariableNames(2) = "VarID";
    T.Properties.VariableNames(3) = "RAJ2000";
    T.Properties.VariableNames(4) = "DEJ2000";
    T.Properties.VariableNames(5) = "Vmag";
    T.Properties.VariableNames(6) = "recno";
    
    idx_valid = ~(strcmp(T.RAJ2000,"") | strcmp(T.DEJ2000,"") | isnan(T.Vmag));
    T = T(idx_valid,:);

    idx_valid = ~(T.Vmag <= 6.00);
    T = T(idx_valid,:);
    
    RAstr  = string(T.RAJ2000);
    DEstr  = string(T.DEJ2000);
    MAG    = T.Vmag;
    
    N = length(RAstr);
    RA_deg = zeros(N,1);
    DE_deg = zeros(N,1);
    
    for i = 1:N
        RA_vals = sscanf(RAstr(i), '%d %d %f');
        h = RA_vals(1); 
        m = RA_vals(2); 
        s = RA_vals(3);
        RA_deg(i) = 15 * (h + m/60 + s/3600);
    
        DE_vals = sscanf(DEstr(i), '%d %d %f');
        deg = DE_vals(1); 
        fir = DE_vals(2); 
        sec = DE_vals(3);
        if deg < 0
            DE_deg(i) = deg - fir/60 - sec/3600;
        else
            DE_deg(i) = deg + fir/60 + sec/3600;
        end
    end
    
    RA  = deg2rad(RA_deg);
    DEC = deg2rad(DE_deg);
    
    % define inertial components
    sx = cos(DEC) .* cos(RA);
    sy = cos(DEC) .* sin(RA);
    sz = sin(DEC);
    
    starVecs_N = [sx sy sz];
    
end
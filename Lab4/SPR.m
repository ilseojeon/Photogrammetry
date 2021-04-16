% Single Photo Resection
% Impyeong Lee
% University of Seoul
% 2007. 6. 4

% Stage 4. Single Photo Resection
% Ilseo Jeon
% Dept. of Geoinformatics, The University of Seoul
% Adv. Photogrammetry Lab4 in Master Course
% 2020. 05. 8

function [EstiEO, SavedEO, SaveEt, SaveVct, SaveVEO, F0]= SPR(IP_ics, GCP, IO, ini_EO, delta)

EO=ini_EO;

for k = 1:100

    % Rotational Matrix    
    R = Rot3D(EO(4), EO(5), EO(6));

    % Derivatives of the Rotational Matrix
    [dR_om, dR_ph, dR_kp] = dRot3D(EO(4), EO(5), EO(6));

    % Observation Equations
    for n = 1:size(IP_ics,2)
        GC = GCP(:,n) - EO(1:3,1);
        ND = R * GC;
        F0(2*n-1:2*n,1) = IO(1:2,1) - IO(3) / ND(3) * ND(1:2,1);
    
        dND(:,1:3) = -R;
        dND(:,4) = dR_om * GC;
        dND(:,5) = dR_ph * GC;
        dND(:,6) = dR_kp * GC;
        
        A(2*n-1:2*n,:) = IO(3) / ND(3)^2 * [-ND(3) 0 ND(1); 0 -ND(3) ND(2)] * dND;
        Y(2*n-1:2*n,1) = [IP_ics(1,n) IP_ics(2,n)]';
    end
    EstiEO(:,k) = EO(:,1);
    
    % Estimate based on Least Squares Adjustment
    dEO = inv(A'*A) * A'*(Y-F0);

    % Update the Initial Appoximations of EO
    EO = EO + dEO;
    
    et = (Y-F0) - A * dEO;
    vct = ( et' * et ) / (size(A,1) - rank(A));
    sdt = sqrt(vct);
    VEO = vct * inv(A'*A);
    
    SavedEO(:,k)=dEO;
    SaveEt=et;
    SaveVct=vct;
    SaveVEO=VEO;
    if norm(dEO) < delta
        break
    end
k
end

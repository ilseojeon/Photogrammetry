function [EO, GP, stat] = RO_F(TP, GP, IO, EO)

% tolerance
delta = 1e-6; %언제쯤 수렴할지 설정, 차이 값
no_TP = size(TP,1); % # of tie points

Dy = eye(size(TP,1)*4); %4n개의 관측 값에 대한 분산, 공분산 행렬을 구하고자 함

for k = 1:100 %무한 루프가 되지 않도록 100까지만 반복
%    p16 관측방정식 세우기
%    keyboard
    for m = 1:2 %image pair
        % Rotational Matrix % omega, phi, kappa
        R{m} = Rot3D(EO{k}(m,4), EO{k}(m,5), EO{k}(m,6));
        % Derivatives of the Rotational Matrix
        [dR_om{m}, dR_ph{m}, dR_kp{m}] = dRot3D(EO{k}(m,4), EO{k}(m,5), EO{k}(m,6));
    end
    
    % Observation Equations %SPR colinearity
    % 변수 선언
    A{k} = zeros(size(TP,1)*4,12+size(TP,1)*3);
    F0{k} = zeros(size(TP,1)*4,1);
    Y{k} = zeros(size(TP,1)*4,1);
    dND_EO = zeros(3,6);
    dND_GP = zeros(3,3);
    
    for n = 1:no_TP
        for m = 1:2 %image pair
            % Row Index
            ri = (m-1)*no_TP*2+2*n-1:(m-1)*no_TP*2+2*n; %m=1일 때, 2n-1:2n, m=2일 때, #ofTP*2+2n-1:#ofTP*2+2n
            
            GC = GP{k}(n,:)' - EO{k}(m,1:3)'; %MP - perspective center
            ND = R{m} * GC; % model point coordinates in CCS
            F0{k}(ri,1) = IO(1:2,1) - IO(3) / ND(3) * ND(1:2,1); % Image point in CCS projected from model point (Computed)
            
            dF_ND = IO(3) / ND(3)^2 * [-ND(3) 0 ND(1); 0 -ND(3) ND(2)]; % Calculating A matrix (design matrix)

            dND_EO(:,1:3) = -R{m};
            dND_EO(:,4) = dR_om{m} * GC;
            dND_EO(:,5) = dR_ph{m} * GC;
            dND_EO(:,6) = dR_kp{m} * GC;

            ci_EO = 6*m-5:6*m;            
            A{k}(ri,ci_EO) =  dF_ND * dND_EO; %row index, eo index에 채워 넣음

            ci_GP = 12+3*n-2:12+3*n;
            dND_GP = R{m}; %-R이 아닌 R이 나옴
            A{k}(ri,ci_GP) = dF_ND * dND_GP;
            
            Y{k}(ri,1) = TP(n, 2*m:2*m+1)';
        end
    end
    y{k} = Y{k}-F0{k};
    % Estimate based on Least Squares Adjustment
%     ci_EO = [4 6 10:12]; % Independent RO img2opk
    ci_EO = [7 9:12]; % Dependent RO img2, x,z,o,p,k
    ci_GP = 13:12+3*no_TP;
    ci = [ci_EO ci_GP];
    [kt{k}, et{k}, vct(k), Dkt{k}] = LSE(A{k}(:,ci), y{k}, Dy); % A에서 실제로 갱신할 포인트에 대해서만 포함 고정한 점 7개

    % Update the Initial Appoximations of EO
    EO{k+1} = EO{k};
    for i = 1:length(ci_EO)
        r = floor((ci_EO(i)-1)/6)+1;
        c = mod(ci_EO(i)-1,6)+1;
        EO{k+1}(r,c) = EO{k}(r,c) + kt{k}(i);
    end
    GP{k+1} = GP{k} + reshape(kt{k}(6:5+3*no_TP),[3 no_TP])';
    
    if norm(kt{k}) < delta
        break
    end
end
stat = {kt, et, vct, Dkt};

end
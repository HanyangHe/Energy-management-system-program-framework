function [Aeq_grid,beq_grid]=LinearEqualityConstraintsForBFS_SOCP_powerbalance(Line,Kday,N_batt,Sinj_bus,BusOfBattery)
% Note: the root bus need to be defined as 1
LD=zeros(size(Line,1),5);
    LD(:,1)=1:size(Line,1);
    LD(:,2)=Line(:,1);
    LD(:,3)=Line(:,2);
    LD(:,4)=Line(:,3);
    LD(:,5)=Line(:,4);

Nbranch = size(LD,1);     
Nbus = max(max(LD(:,2:3)));

% index of the decision variable, which is the column index of Aeq
idx_P = 1:Nbranch;
idx_l = Nbranch+1:2*Nbranch;
idx_v = 2*Nbranch+1:3*Nbranch;

Aeq = zeros(2*Nbranch, 3*Nbranch); 
beq = zeros(2*Nbranch,1);

offset = 2*Kday*N_batt + 2*Kday;
totalVars=3*Nbranch;

Aeq_grid=[];
beq_grid=[];
for step=1:Kday
    for e = 1:Nbranch
	    INbus  = LD(e,2);
	    OUTbus = LD(e,3);
	    r_e    = LD(e,4);
    
	    % --- row number
	    row_pow = 2*e-1; % power balance equation
	    row_vlt = 2*e;   % voltage equation
    
	    % ----------（1）power balance equation----------
	    Aeq(row_pow, idx_P(e)) = 1;  % branch power at the current branch
    
	    % downstream branch start at OUTbus
	    child_branch = find(LD(:,2)==OUTbus);
	    for k = child_branch(:)'
		    Aeq(row_pow, idx_P(k)) = -1;
	    end
	    Aeq(row_pow, idx_l(e)) = -r_e;

    
	    % if there is node injection power at the right side (ejection side),
        % and them here with "-".
	    beq(row_pow) = -real(Sinj_bus(OUTbus,step)); % node injection power
    
	    % ----------（2）voltage drop equation-----------
        % save all downstream voltage 2~Nbus at the last Nbranch colcumn
	    Aeq(row_vlt, idx_v(OUTbus-1)) = 1;
        if INbus==1
            % root voltage is given at 1pu, and move into the right side
            % coefficient vector
		    beq(row_vlt) = 1;
        else
	        Aeq(row_vlt, idx_v(INbus-1))  = -1;
        end
	    Aeq(row_vlt, idx_l(e))      = -r_e^2;
	    Aeq(row_vlt, idx_P(e))      = 2*r_e;
    end

    Aeq_extd = zeros(size(Aeq,1), offset+totalVars*Kday); 
    Aeq_extd(:,offset+1+(step-1)*totalVars:offset+step*totalVars)=Aeq;

    Nbatt = length(BusOfBattery);
    for bat_idx = 1:Nbatt
        ii = BusOfBattery(bat_idx);  % battery node number
        base = 2*Kday*(bat_idx-1);   % start location of the battery variable sequence
        col_Pdis = base + step;
        col_Pch  = base + Kday + step;
    
        branchNum = find(Line(:,2)==ii,1); % battery branch (battery node is the downstream node, branch number is the row number of that branch in the Line matrix)
        Aeq_extd(2*branchNum-1, col_Pdis) = 1;
        Aeq_extd(2*branchNum-1, col_Pch)  = -1;
    end

    Aeq_grid=[Aeq_grid;Aeq_extd];

    beq_grid=[beq_grid;beq];
end

end
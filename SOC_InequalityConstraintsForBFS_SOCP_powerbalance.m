function [A_SOC,b_SOC]=SOC_InequalityConstraintsForBFS_SOCP_powerbalance(Kday,decisionLength,SOCmax,SOCmin,SOC0,EbattMax,inta_batt_ch,inta_batt_dis,N_batt,dt_sys)

idx_A=0;
idx_b=0;

b_SOC=zeros(Kday*2*N_batt,1);
for ii=1:N_batt
    b_SOC(idx_b+1:idx_b+Kday,1)=(SOCmax-SOC0(ii,1))*ones(Kday,1);idx_b=idx_b+Kday;
    b_SOC(idx_b+1:idx_b+Kday,1)=(SOC0(ii,1)-SOCmin)*ones(Kday,1);idx_b=idx_b+Kday;
    b_SOC(idx_b,1)=SOC0(ii,1)-SOC0(ii,1);%total charge >= total discharge in a day
end

A_SOC=zeros(length(b_SOC),decisionLength);
for ii=1:N_batt
    for i=1:Kday
        A_SOC(idx_A+i,idx_A+1:idx_A+i)=-1/(EbattMax*inta_batt_dis)*dt_sys;
        A_SOC(idx_A+i,idx_A+Kday+1:idx_A+Kday+i)=inta_batt_ch/(EbattMax)*dt_sys;
        A_SOC(idx_A+Kday+i,idx_A+1:idx_A+i)=1/(EbattMax*inta_batt_dis)*dt_sys;
        A_SOC(idx_A+Kday+i,idx_A+Kday+1:idx_A+Kday+i)=-inta_batt_ch/(EbattMax)*dt_sys;
    end
    idx_A=idx_A+2*Kday;
end
end
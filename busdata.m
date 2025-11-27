function [Bus_Prate,idx_P_PV,idx_Pbatt,idx_Pbusi,idx_Presi] = busdata(residential_node,business_node,BusOfPV,BusOfBattery,BusNum,PF,residentialPower,businessPower,PVInv_S_rate,BatteryInv_S_rate)% This file is useless for the EMS now

% S_base=x;%VA (1MVA)
% V_base=x;%V (low level distribution network)

% resiLoadNum=length(residential_node);
% busiLoadNum=length(business_node);
% SolarNum=length(BusOfPV);

% residentialPower=[x;x;x;...];%MW
% businessPower=[x;x;...];%MW

% PVInv_S_rate=x;%MW
% BatteryInv_S_rate=x;%MW

Bus_Prate=[];
idx_P_PV=5;
idx_Pbatt=7;
idx_Pbusi=9;
idx_Presi=11;
for ii=1:BusNum
    % [bus type Psys Qsys P_PV Q_PV P_batt Q_batt P_busi Q_busi P_resi Q_resi]
    Bus_Prate=[Bus_Prate;ii 0 0 0 0 0 0 0 0 0 0 0];
end
Bus_Prate(1,2)=1;
Bus_Prate(2:end,2)=3;

if length(residential_node)>0
Bus_Prate(residential_node,idx_Presi)=-residentialPower;
Bus_Prate(residential_node,idx_Presi+1)=-residentialPower*tan(acos(PF));
else
Bus_Prate(residential_node,idx_Presi)=-0;
Bus_Prate(residential_node,idx_Presi+1)=-0;
end

if length(business_node)>0
Bus_Prate(business_node,idx_Pbusi)=-businessPower;
Bus_Prate(business_node,idx_Pbusi+1)=-businessPower*tan(acos(PF));
else
Bus_Prate(business_node,idx_Pbusi)=-0;
Bus_Prate(business_node,idx_Pbusi+1)=-0;
end

if length(BusOfPV)>0
Bus_Prate(BusOfPV,idx_P_PV)=PVInv_S_rate;
else
Bus_Prate(BusOfPV,idx_P_PV)=0;
end

if length(BusOfBattery)>0
Bus_Prate(BusOfBattery,idx_Pbatt)=BatteryInv_S_rate;
else
Bus_Prate(BusOfBattery,idx_Pbatt)=0;
end
end

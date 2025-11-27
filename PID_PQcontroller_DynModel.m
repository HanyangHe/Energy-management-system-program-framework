function [Edq,dotX]=PID_PQcontroller_DynModel(Ref,FDBK,Xm_g,Zmv_g,controlParams,X)

% control parameters
KP_PLL_pq=controlParams(1,1);
KI_PLL_pq=controlParams(2,1);
KI_outer_p=controlParams(3,1);
KI_outer_q=controlParams(4,1);
KI_inner_p=controlParams(5,1);
KI_inner_q=controlParams(6,1);
KP_outer_p=controlParams(7,1);
KP_outer_q=controlParams(8,1);
KP_inner_p=controlParams(9,1);
KP_inner_q=controlParams(10,1);

Kcross=0.00;

% input
Pref = Ref(1,1);
Qref = Ref(2,1);

Vinv_g=FDBK(1,1);
Vpcc=FDBK(2,1);

% fdbk voltage and current calculation & xy->dq
Theta_g = mod(X(2,1),2*pi);
Vqs_g = -real(Vpcc)*sin(Theta_g) + imag(Vpcc)*cos(Theta_g);%Theta_g is I or V angle?
I_in_g = ((Vinv_g - Vpcc)/(1j*Xm_g+Zmv_g));
Id_g = real(I_in_g(1,1))*cos(Theta_g) + imag(I_in_g(1,1))*sin(Theta_g);
Iq_g = -real(I_in_g(1,1))*sin(Theta_g) + imag(I_in_g(1,1))*cos(Theta_g);
Pg=real(Vpcc*conj(I_in_g));
Qg=imag(Vpcc*conj(I_in_g));

% Define the xdot = f(x,y) functions
dotX(1,1) = KI_PLL_pq*Vqs_g;
dotX(2,1) = KP_PLL_pq*Vqs_g + X(1,1);

dotX(3,1) = (Pref - Pg)*KI_outer_p;
id_ref_g = KP_outer_p*(Pref - Pg) + X(3,1);

dotX(5,1) = (-Qref + Qg)*KI_outer_q;
iq_ref_g = KP_outer_q*(-Qref + Qg) + X(5,1);

% ----- reference current limitation -----
I_dq_ref_max = 1.2; % pu
I_ref_norm = sqrt(id_ref_g^2 + iq_ref_g^2);
if I_ref_norm > I_dq_ref_max
    scale = I_dq_ref_max / I_ref_norm;
    id_ref_g = id_ref_g * scale;
    iq_ref_g = iq_ref_g * scale;
end
% -----------------------

dotX(4,1) = (id_ref_g - Id_g)*KI_inner_p;

dotX(6,1) = (iq_ref_g - Iq_g)*KI_inner_q;

% virtual impedence voltage droop (improve the robustness)
% Rv_pu=real(Zmv_pu);
% Xv_pu=imag(Zmv_pu);
% Vd_v =  Rv_pu*(real(I_in_g)*cos(Theta_g) + imag(I_in_g)*sin(Theta_g)) ...
%        - Xv_pu*(-real(I_in_g)*sin(Theta_g)+imag(I_in_g)*cos(Theta_g));
% 
% Vq_v =  Rv_pu*(-real(I_in_g)*sin(Theta_g)+imag(I_in_g)*cos(Theta_g)) ...
%        + Xv_pu*( real(I_in_g)*cos(Theta_g)+imag(I_in_g)*sin(Theta_g) );

% output voltage modulation index
Edr_g = KP_inner_p*(id_ref_g - Id_g) + X(4,1)+Kcross*dotX(2,1)*Xm_g*Iq_g ;
Eqr_g = KP_inner_q*(iq_ref_g - Iq_g) + X(6,1)-Kcross*dotX(2,1)*Xm_g*Id_g ;
Edq=[Edr_g;Eqr_g];
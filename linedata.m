function [Line_Data, LineLimits ]= linedata(Topology,Z_base,I_base)

% cupper 3ph
z_AWG1000=0.010+1j*0.08;%ohm/km
z_AWG750=0.025+1j*0.08;%ohm/km
z_AWG600=0.03+1j*0.08;%ohm/km
z_AWG500=0.08+1j*0.08;%ohm/km
z_AWG350=0.11+1j*0.08;%ohm/km
z_AWG250=0.15+1j*0.08;%ohm/km
z_AWG4_0=0.17+1j*0.08;%ohm/km
z_AWG3_0=0.23+1j*0.08;%ohm/km
z_AWG2_0=0.29+1j*0.08;%ohm/km
z_AWG1_0=0.35+1j*0.08;%ohm/km
z_AWG2=0.55+1j*0.08;%ohm/km
z_AWG3=0.90+1j*0.08;%ohm/km
z_AWG4=1.13+1j*0.08;%ohm/km
z_AWG6=1.40+1j*0.08;%ohm/km
z_AWG8=2.18+1j*0.08;%ohm/km
z_AWG10=3.45+1j*0.08;%ohm/km

% branch current limit (75°C)
lmt_AWG1000=545;%A
lmt_AWG750=475;%A
lmt_AWG600=420;%A
lmt_AWG500=380;%A
lmt_AWG350=310;%A
lmt_AWG250=255;%A
lmt_AWG4_0=230;%A
lmt_AWG3_0=200;%A
lmt_AWG2_0=175;%A
lmt_AWG1_0=150;%A
lmt_AWG2=115;%A
lmt_AWG3=100;%A
lmt_AWG4=85;%A
lmt_AWG6=65;%A
lmt_AWG8=50;%A
lmt_AWG10=35;%A

BranchLength_long=120/1000;%km
BranchLength_main=80/1000;%km
BranchLength_sub=30/1000;%km

%            |From  |     |To    |     |R      |      |X     | 
%            |Bus   |     |Bus   |     |p.u.   |      |p.u.  | 
% Line_Data = [   1           3          0.2587500    0.0153750 ;
%                 3           6          0.4802811    0.0708750 ;
%                 3	        10	       0.0621250    0.0181562 ;
%                 10	        2	       0.6918750    0.0176250 ;
%                 3           4          0.124250     0.0363124 ;
%                 4           5          0.1633125    0.0151875 ;
%                 4           7          0.1863750    0.0544686 ;
%                 7           9          0.6918750    0.0176250 ;
%                 7           8          0.320875     0.0335312
%             ];


if strcmp(Topology, '10bus-case1')
Line_Data = [   1           3          BranchLength_main*real(z_AWG500)   BranchLength_main*imag(z_AWG500);
                3           6          BranchLength_main*real(z_AWG8)   BranchLength_main*imag(z_AWG8);
                3	        10	       BranchLength_main*real(z_AWG4_0)   BranchLength_main*imag(z_AWG4_0);
                10	        2	       BranchLength_main*real(z_AWG10)   BranchLength_main*imag(z_AWG10);
                3           4          BranchLength_main*real(z_AWG350)   BranchLength_main*imag(z_AWG350);
                4           5          BranchLength_main*real(z_AWG6)   BranchLength_main*imag(z_AWG6);
                4           7          BranchLength_main*real(z_AWG3_0)   BranchLength_main*imag(z_AWG3_0);
                7           9          BranchLength_main*real(z_AWG10)   BranchLength_main*imag(z_AWG10);
                7           8          BranchLength_main*real(z_AWG6)   BranchLength_main*imag(z_AWG6);
            ];
Line_Data(:,3:4)=Line_Data(:,3:4)/Z_base;

LineLimits = [lmt_AWG500;
            lmt_AWG8;
            lmt_AWG4_0;
            lmt_AWG10;
            lmt_AWG350;
            lmt_AWG6;
            lmt_AWG3_0;
            lmt_AWG10;
            lmt_AWG6
            ]/I_base;

elseif strcmp(Topology, '10bus-case2')

Line_Data = [   1           3          BranchLength_main*real(z_AWG2_0)/2   BranchLength_main*imag(z_AWG2_0)/2;
                3           6          BranchLength_main*real(z_AWG8)   BranchLength_main*imag(z_AWG8);
                3	        10	       BranchLength_main*real(z_AWG4_0)   BranchLength_main*imag(z_AWG4_0);
                10	        2	       BranchLength_main*real(z_AWG10)   BranchLength_main*imag(z_AWG10);
                3           4          BranchLength_main*real(z_AWG1_0)/2   BranchLength_main*imag(z_AWG1_0)/2;
                4           5          BranchLength_main*real(z_AWG6)   BranchLength_main*imag(z_AWG6);
                4           7          BranchLength_main*real(z_AWG2_0)   BranchLength_main*imag(z_AWG2_0);
                7           9          BranchLength_main*real(z_AWG10)   BranchLength_main*imag(z_AWG10);
                7           8          BranchLength_main*real(z_AWG6)   BranchLength_main*imag(z_AWG6);
            ];
Line_Data(:,3:4)=Line_Data(:,3:4)/Z_base;

LineLimits = [lmt_AWG2_0*2;
            lmt_AWG8;
            lmt_AWG4_0;
            lmt_AWG10;
            lmt_AWG1_0*2;
            lmt_AWG6;
            lmt_AWG2_0;
            lmt_AWG10;
            lmt_AWG6
            ]/I_base;

elseif strcmp(Topology, 'case name')
% you can design a new case follow the structure above
end

end


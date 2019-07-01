function Smem = calcSmem(Vint)

% Takes a given Vint to calculate the Smem term (volume of liposome to
% total membrane area).

% Based on measurements from liposomes transport assays of EmrE (Robinson
% et al. 2017, 10.1073/pnas.1708671114)

% (C) Hallie Hanson

Ldia = 220e-8;                                     % Dlip in units of dm
Lvol = (4/3)*3.14159*(((Ldia-(5e-8))/2)^3);        % Lvol in units of L
Lmem = 4*3.14159*((Ldia/2)^2+((Ldia/2)-(5e-8))^2); % Lmem in units of dm^2
Lnum = Vint/Lvol;
 
Smem = Lmem*Lnum;                                  % Smem in units of dm^2

end



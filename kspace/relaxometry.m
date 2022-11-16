function relaxometry

B0 = 0.55 *1e4
gamma_rad = 267.5222 * 1e6 * 1e-4  ; % rad/s/Gauss
gamma = gamma_rad/2/pi ; % Hz/Gauss

w0 = gamma*B0

% from 
% ﻿1. Rooney WD, Johnson G, Li X, Cohen ER, Kim SG, Ugurbil K, Springer CS. 
% Magnetic field and tissue dependencies of human brain longitudinal 1H2O relaxation in vivo.
% Magn. Reson. Med. 2007;57:308–318. 
% doi: 10.1002/mrm.21122.

putamenGM_T1 = 0.00116*(gamma*B0)^0.376
WM_T1 = 0.00071*(gamma*B0)^0.382
venous_T1 = 0.00335*(gamma*B0)^0.340 
function [rx2sat_range, ipp2rx_range, veff] = get_veff(time_utc, ...
    rx_traj_llh, rx_vel_ned, sat_traj_llh, sat_vel_ned, h_intercept, ...
    drift_vel_ned)
% References:
% [1] Jiao, Yu, Dongyang Xu, Charles L. Rino, Yu T. Morton, and Charles S.
%     Carrano. “A Multifrequency GPS Signal Strong Equatorial Ionospheric
%     Scintillation Simulator: Algorithm, Performance, and 
%     Characterization.” IEEE Transactions on Aerospace and Electronic 
%     Systems 54, no. 4 (August 2018): 1947–65. 
%     https://doi.org/10.1109/TAES.2018.2805232.
% [2] Rino, Charles. The Theory of Scintillation with Applications in 
%     Remote Sensing. John Wiley & Sons, 2011.
% [3] Vasylyev, Dmytro, Yannick Béniguel, Wilken Volker, Martin Kriegel, 
%     and Jens Berdermann. “Modeling of Ionospheric Scintillation.” 
%     Journal of Space Weather and Space Climate 12 (2022): 22.

%% Get satellite trajectory and velocity in ENU (origin in the receiver)
%{
    NOTE:
    The ENU (east-north-up) local frame is being used to represent the
    satellite trajectory. For out case, the origin is the receiver.

                 Z(upwards)
                 |    Y(northwards)
                 |   /
                 |  /
                 | /
                 |/
          ───────+───────X(eastwards)
                /|
               / |
              /  |
             /   |

    • X-axis: points eastwards
    • Y-axis: points nothwards
    • Z-axis: points upwards
    • Origin: is the receiver position, which is possibly moving
%}

% Setellite LLH to ENU with receiver position as origin of the TCS frame
sat_traj_enu_rx = llh2tcsT(sat_traj_llh, rx_traj_llh);

% Convert satellite velocity from NED to ENU
sat_vel_enu = zeros(size(sat_vel_ned));
sat_vel_enu(1,:) = sat_vel_ned(2,:);
sat_vel_enu(2,:) = sat_vel_ned(1,:);
sat_vel_enu(3,:) =-sat_vel_ned(3,:);

% Convert receiver velocity from NED to ENU
rx_vel_enu = zeros(size(rx_vel_ned));
rx_vel_enu(1,:) = rx_vel_ned(2,:);
rx_vel_enu(2,:) = rx_vel_ned(1,:);
rx_vel_enu(3,:) =-rx_vel_ned(3,:);

%% Get receiver-to-satellite parameters (azimuth, elevation, and range)

% elevation angle (θₛₐₜ) -> see Rino eq. (4.69)
rx2sat_elev = atan2(sat_traj_enu_rx(3,:),sqrt(sat_traj_enu_rx(1,:).^2+sat_traj_enu_rx(2,:).^2));
% TCS origin to sat range (rₛₐₜ) -> see Rino eq. (4.68)
rx2sat_range = sqrt(sat_traj_enu_rx(1,:).^2+sat_traj_enu_rx(2,:).^2+sat_traj_enu_rx(3,:).^2);
% azimuth (ϕₛₐₜ) -> see Rino eq. (4.69)
rx2sat_azim = atan2(sat_traj_enu_rx(1,:),sat_traj_enu_rx(2,:));

% receiver-to-satellite unit vector in TCS (origin in the receiver)
u_rx2sat_tcs = sat_traj_enu_rx./repmat(rx2sat_range,3,1);

%% Compute IPP in LLH and receiver trajectory in ENU (origin in the IPP)

% get IPP trajectory in LLH (geodetic)
ipp_traj_llh = findIntercept_LEO(h_intercept,u_rx2sat_tcs,rx2sat_range,rx_traj_llh);

% get satellite trajectory in ENU (origin is the IPP)
rx_traj_enu_ipp = llh2tcsT(rx_traj_llh, ipp_traj_llh);

%% Convert satellite, receiver and IPP trajectories and velocities from ENU to DES
%{
    NOTE:
    The DES (down-east-south) local frame is being used to represent the
    receiver trajectory, where the origin is the IPP.

                 |    
                 |   /
                 |  /
                 | /
                 |/
          ───────+───────Y(eastwards)
                /|
               / |
              /  |
             /   |
Z(southwards)    X(downwards)

    • X-axis: points eastwards
    • Y-axis: points nothwards
    • Z-axis: points upwards
    • Origin: is the IPP position, which is moving.
%}

% convert receiver trajectory from ENU to DES
rx_traj_des_ipp(1,:) = -rx_traj_enu_ipp(3,:);  %Component 1 is -z tcs (Downward)  
rx_traj_des_ipp(2,:) =  rx_traj_enu_ipp(1,:);  %Component 2 is  x tcs (Eastward)
rx_traj_des_ipp(3,:) = -rx_traj_enu_ipp(2,:);  %Component 3 is -y tcs (Southward)

% convert receiver trajectory from ENU to DES
rx_vel_des(1,:) =-rx_vel_enu(3,:);  %Component 1 is -z tcs (Downward)  
rx_vel_des(2,:) = rx_vel_enu(1,:);  %Component 2 is  x tcs (Eastward)
rx_vel_des(3,:) =-rx_vel_enu(2,:);  %Component 3 is -y tcs (Southward)

% convert satellite trajectory from ENU to DES
sat_vel_des(1,:) =-sat_vel_enu(3,:);  %Component 1 is -z tcs (Downward)  
sat_vel_des(2,:) = sat_vel_enu(1,:);  %Component 2 is  x tcs (Eastward)
sat_vel_des(3,:) =-sat_vel_enu(2,:);  %Component 3 is -y tcs (Southward)

% convert drift trajectory from NED to DES
drift_vel_des = zeros(size(drift_vel_ned));
drift_vel_des(1) = drift_vel_ned(3);
drift_vel_des(2) = drift_vel_ned(2);
drift_vel_des(3) =-drift_vel_ned(1);

%% Get receiver-to-IPP parameters (range and angles)

% range from the IPP to the receiver
ipp2rx_range  =sqrt(rx_traj_des_ipp(1,:).^2+rx_traj_des_ipp(2,:).^2+rx_traj_des_ipp(3,:).^2);
uk_xyzp  = rx_traj_des_ipp./repmat(ipp2rx_range,3,1);
% Propagation angles
% NOTE: See figure 2 in [1]
% θ TODO: check the angles
ipp2rx_east2south_angle = atan2(sqrt(rx_traj_des_ipp(2,:).^2+rx_traj_des_ipp(3,:).^2),rx_traj_des_ipp(1,:));
% ϕ
ipp2rx_to_down_angle   = atan2(rx_traj_des_ipp(3,:),rx_traj_des_ipp(2,:));

%% satellite signal scan velocity at the IPP
% Satellite velocity scale vector at penetration point
sat_vel_scale=ipp2rx_range./rx2sat_range;
% Receiver velocity scale vector at penetration point
rx_vel_scale=(rx2sat_range-ipp2rx_range)./rx2sat_range;
%{ 
    get IPP scan velocity
    NOTE: We didn't understand why it is like that
    SEE: [1] eq. (13) and (14)
    SEE: Eq. (4.4), (4.41), (4.42) in Rino's book
%}

ipp_scan_vel_des(1,:)=(sat_vel_des(1,:).*sat_vel_scale)+(rx_vel_des(1,:).*rx_vel_scale);
ipp_scan_vel_des(2,:)=(sat_vel_des(2,:).*sat_vel_scale)+(rx_vel_des(2,:).*rx_vel_scale);
ipp_scan_vel_des(3,:)=(sat_vel_des(3,:).*sat_vel_scale)+(rx_vel_des(3,:).*rx_vel_scale);

%% Get the megnetic field's unit vector in NED
nsamps=length(time_utc);
mag_field_ned=zeros(3,nsamps);

% NOTE: igrf() outputs the magnetic field in NED
dtr=pi/180;
[mag_field_ned(1,:), mag_field_ned(2,:), mag_field_ned(3,:)] =igrf(datenum(time_utc(1)), ...
               ipp_traj_llh(1,:)/dtr, ipp_traj_llh(2,:)/dtr, ipp_traj_llh(3,:)/1000);  

% Get absolute velue of the magnitic field vector
mag_field_abs=sqrt(mag_field_ned(1,:).^2+mag_field_ned(2,:).^2+mag_field_ned(3,:).^2);
% Get magnetic field's unit vector from NED to DES
u_mag_field_ned(1,:)=  mag_field_ned(3,:)./mag_field_abs; % northward to downward
u_mag_field_ned(2,:)=  mag_field_ned(2,:)./mag_field_abs; % Eastward to Eastward
u_mag_field_ned(3,:)= -mag_field_ned(1,:)./mag_field_abs; % Downward to southward

%% Get magnetic field's parmeters

thetaB = atan2(sqrt(u_mag_field_ned(2,:).^2+u_mag_field_ned(3,:).^2),u_mag_field_ned(1,:));
phiB   = atan2(u_mag_field_ned(3,:),u_mag_field_ned(2,:));
cosBP  = abs(dot(uk_xyzp,u_mag_field_ned));

%% Get the effective scan velocity
%{
    SEE: Refer to [2, Appendix A.3]
    SEE: Refer to [3, Section 4.2]
%}

% SEE: Refer to [1, Equation 13]
vky = drift_vel_des(2)-ipp_scan_vel_des(2,:)+tan(ipp2rx_east2south_angle).*cos(ipp2rx_to_down_angle).*(-drift_vel_des(1)+ipp_scan_vel_des(1,:));
% SEE: Refer to [1, Equation 14]
vkz = drift_vel_des(3)-ipp_scan_vel_des(3,:)+tan(ipp2rx_east2south_angle).*sin(ipp2rx_to_down_angle).*(-drift_vel_des(1)+ipp_scan_vel_des(1,:));

% NOTE: Unidentified paramater name. Review this later.
gam_b=0; 
% Principal axis enlongation
a=50; 
% Transverse axis enlongation
b=1;
% Change magnetic angle ref to horizontal.  CLR 2/13/2016
thetaB=thetaB+pi/2;
% Get the anisotropy factors (A, B, C)
[A,B,C] = ABC(ipp2rx_east2south_angle,ipp2rx_to_down_angle,thetaB,phiB,gam_b,a,b);

% SEE: (4.48) in Rino's book
veff=sqrt((C.*vky.^2-B.*vky.*vkz+A.*vkz.^2)./(A.*C-B.^2/4));
end
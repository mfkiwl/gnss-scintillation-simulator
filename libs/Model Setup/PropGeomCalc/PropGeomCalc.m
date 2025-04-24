function        satGEOM_struct=PropGeomCalc(GPSTime,UT_time,eph,rx_traj_llh,h_intercept,rx_vel_enu,varargin)
%USAGE:   satGEOM_struct=PropGeomCalc(GPSTime,UT_time,eph,rx_traj_llh,h_intercept,rx_v,varargin)
%
%INPUTS:
%   GPSTIme = [GPS week; seconds of week]
%   UT_time    =  6 element date time for IGRF 
%   eph             = RINEX navagation format ephemeris file
%   rx_traj_llh    = station  [latitude (rad); longitude( rad); height (m)]
%   Drift            = [downward,eastward,southward] drift mps  
%
%OUTPUTS:
%   *************ecef Coordinates**********************************************************
%   xsat_ecef, vsat_ecef = satellite state vector in ecef coordinates from spg4 (3XN)
%   ************GPS Coordinates***********************************************************
%   sat_llh            = satellite geodetic coordinates (3XN)
%   ************Station TCS coordinates***************************************************
%   sat_tcs,  vsat_tcs = satellite state vector in receiver tcs system at rx_traj_llh (3XN)
%   sat_rng,  sat_rdot = satellite rang & range rate (=> sat_tcs) (3XN)
%   sat_elev, sat_phi  = satellite elevation & true bearing (NX1)
%   *************Propagation Reference Coordinates at penetration point*******************
%   xyzp               = propagation coordinates with origin at h_intercept (3xN)
%   thetap, phip       = polar angles wrt x  (phip  cw from y-axis (NX1)
%                                          theta cw from x-axis
%   rngp               = range from receiver to intercept point (Nx1)
%   uk_xyzp         = unit vector pointing along propagation direction (3XN)
%   s                     = unit magnetic field vector xp,yp,zp system (3xN)
%   thetaB,psiB   = polar angles wrt xp (Nx1)
%   vp                 = penetration point velocity <= satellite motion 
%   vk                 = apparent velocity in measurement plane
%   sat_utsec          = time (sec)
%
%Libraries:   GPS_CoordinateXforms IGRF_Compston 
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
%Author:
%Charles Rino
%Rino Consulting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%VERSION: February 13, 2016%%%%%%%%%%%%%%%%%%%%

dtr=pi/180;

%% Initalize drift velocity
if isempty(varargin)
    drift_vel = [0;0;0];
else
    drift_vel = varargin{1};
end
nsamps=length(GPSTime(1,:));

%% Get satellite trajectory and velocity in ENU (origin in the receiver)
%{
    NOTE:
    The ENU (east-north-up) local frame is being used to represent the
    satellite trajectory, where the origin is the receiver.

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

    • X-axis: points to eastwards
    • Y-axis: points to nothwards
    • Z-axis: points to upwards
    • Origin: is the receiver position, which is moving
%}

sat_traj_ecef=zeros(3,nsamps);
sat_vel_ecef=zeros(3,nsamps);

% Calculate satellite ECEF position & velocity
for nsamp=1:nsamps
    [sat_traj_ecef(:,nsamp), sat_vel_ecef(:,nsamp)] = satposvel(GPSTime(2,nsamp),eph);
end

% Satellite ECEF to LLH (geodetic)
sat_traj_llh=ecef2llhT(sat_traj_ecef);
% Setellite LLH to ENU with receiver position as origin of the TCS frame
sat_traj_enu_rx = llh2tcsT(sat_traj_llh, rx_traj_llh);

% Convert satellite ECEF velocity to ENU
sat_vel_enu = zeros(3, size(rx_traj_llh,2));
for i = 1:size(rx_traj_llh, 2)
    D = Rotate_ecef2tcs(rx_traj_llh(:,i));
    sat_vel_enu(:,i) = D*sat_vel_ecef(:,i) - rx_vel_enu;
end

%% Get receiver-to-satellite parameters (azimuth, elevation, and range)

% elevation angle (θₛₐₜ) -> see Rino eq. (4.69)
rx2sat_elev = atan2(sat_traj_enu_rx(3,:),sqrt(sat_traj_enu_rx(1,:).^2+sat_traj_enu_rx(2,:).^2));
% TCS origin to sat range (rₛₐₜ) -> see Rino eq. (4.68)
rx2sat_range = sqrt(sat_traj_enu_rx(1,:).^2+sat_traj_enu_rx(2,:).^2+sat_traj_enu_rx(3,:).^2);
% azimuth (ϕₛₐₜ) -> see Rino eq. (4.69)
rx2sat_azim =atan2(sat_traj_enu_rx(1,:),sat_traj_enu_rx(2,:));

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

    • X-axis: points to eastwards
    • Y-axis: points to nothwards
    • Z-axis: points to upwards
    • Origin: is the IPP position, which is moving.
%}

% convert receiver trajectory from ENU to DES
rx_traj_des_ipp(1,:)=-rx_traj_enu_ipp(3,:);  %Component 1 is -z tcs (Downward)  
rx_traj_des_ipp(2,:)= rx_traj_enu_ipp(1,:);  %Component 2 is  x tcs (Eastward)
rx_traj_des_ipp(3,:)=-rx_traj_enu_ipp(2,:);  %Component 3 is -y tcs (Southward)

% convert receiver trajectory from ENU to DES
rx_vel_des(1,:)=-rx_vel_enu(3,:);  %Component 1 is -z tcs (Downward)  
rx_vel_des(2,:)= rx_vel_enu(1,:);  %Component 2 is  x tcs (Eastward)
rx_vel_des(3,:)=-rx_vel_enu(2,:);  %Component 3 is -y tcs (Southward)

% convert satellite trajectory from ENU to DES
sat_vel_des(1,:)=-sat_vel_enu(3,:);  %Component 1 is -z tcs (Downward)  
sat_vel_des(2,:)= sat_vel_enu(1,:);  %Component 2 is  x tcs (Eastward)
sat_vel_des(3,:)=-sat_vel_enu(2,:);  %Component 3 is -y tcs (Southward)

%% Get receiver-to-IPP parameters (range and angles)

% range from the IPP to the receiver
rp  =sqrt(rx_traj_des_ipp(1,:).^2+rx_traj_des_ipp(2,:).^2+rx_traj_des_ipp(3,:).^2);
uk_xyzp  = rx_traj_des_ipp./repmat(rp,3,1);
% Propagation angles
% NOTE: See figure 2 in [1]
% ϕ
ipp2rx_east2south_angle = atan2(sqrt(rx_traj_des_ipp(2,:).^2+rx_traj_des_ipp(3,:).^2),rx_traj_des_ipp(1,:));
% θ
ipp2rx_to_down_angle   = atan2(rx_traj_des_ipp(3,:),rx_traj_des_ipp(2,:));

%% satellite signal scan velocity at the IPP
% Satellite velocity scale vector at penetration point
sat_vel_scale=rp./rx2sat_range;
% Receiver velocity scale vector at penetration point
rx_vel_scale=(rx2sat_range-rp)./rx2sat_range;
%{ 
    get IPP scan velocity
    NOTE: We didn't understand why it is like that
    SEE: [1] eq. (13) and (14)
    SEE: Eq. (4.4), (4.41), (4.42) in Rino's book
%}

ipp_scan_vel(1,:)=(sat_vel_des(1,:).*sat_vel_scale)+(rx_vel_des(1,:).*rx_vel_scale);
ipp_scan_vel(2,:)=(sat_vel_des(2,:).*sat_vel_scale)+(rx_vel_des(2,:).*rx_vel_scale);
ipp_scan_vel(3,:)=(sat_vel_des(3,:).*sat_vel_scale)+(rx_vel_enu(3,:).*rx_vel_scale);

%% Get the megnetic field's unit vector in NED
%Magnetic field at penetration points 
time=datenum(UT_time(1),UT_time(2),UT_time(3));
mag_field=zeros(3,nsamps);

% NOTE: igrf() outputs the magnetic field in NED
[mag_field(1,:), mag_field(2,:), mag_field(3,:)] =igrf(time,...
               ipp_traj_llh(1,:)/dtr, ipp_traj_llh(2,:)/dtr, ipp_traj_llh(3,:)/1000);  

% Get absolute velue of the magnitic field vector
mag_field_abs=sqrt(mag_field(1,:).^2+mag_field(2,:).^2+mag_field(3,:).^2);
% Get magnetic field's unit vector from NED to DES
u_mag_field(1,:)=  mag_field(3,:)./mag_field_abs; % northward to downward
u_mag_field(2,:)=  mag_field(2,:)./mag_field_abs; % Eastward to Eastward
u_mag_field(3,:)= -mag_field(1,:)./mag_field_abs; % Downward to southward

%% Get magnetic field's parmeters

thetaB = atan2(sqrt(u_mag_field(2,:).^2+u_mag_field(3,:).^2),u_mag_field(1,:));
phiB   = atan2(u_mag_field(3,:),u_mag_field(2,:));
cosBP  = abs(dot(uk_xyzp,u_mag_field));

%% Get the effective scan velocity
%{
    SEE: Refer to [2, Appendix A.3]
    SEE: Refer to [3, Section 4.2]
%}

% SEE: Refer to [1, Equation 13]
vky = drift_vel(2)-ipp_scan_vel(2,:)+tan(ipp2rx_east2south_angle).*cos(ipp2rx_to_down_angle).*(-drift_vel(1)+ipp_scan_vel(1,:));
% SEE: Refer to [1, Equation 14]
vkz = drift_vel(3)-ipp_scan_vel(3,:)+tan(ipp2rx_east2south_angle).*sin(ipp2rx_to_down_angle).*(-drift_vel(1)+ipp_scan_vel(1,:));


% NOTE: (Rodrigo): Unidentified paramater name. Review this later.
gam_b=0; 
% Principal axis enlongation
a=50; 
% Transverse axis enlongation
b=1;
thetaB=thetaB+pi/2;    %Change magnetic angle ref to horizontal.  CLR 2/13/2016
    [A,B,C]=ABC(ipp2rx_east2south_angle,ipp2rx_to_down_angle,thetaB,phiB,gam_b,a,b);

% SEE: (4.48) in Rino's book
veff=sqrt((C.*vky.^2-B.*vky.*vkz+A.*vkz.^2)./(A.*C-B.^2/4));

satGEOM_struct=struct('eph',eph,...
    'UT_time',UT_time,'rx_traj_llh',rx_traj_llh',...
    'h_intercept',h_intercept','a',a,'b',b,'gam_b',gam_b,'Drift',drift_vel,...
    'xsat_ecef',sat_traj_ecef,'vsat_tcs',sat_vel_enu,'sat_llh',sat_traj_llh, 'sat_tcs',sat_traj_enu_rx,...
    'rs',rx2sat_range, 'sat_elev',rx2sat_elev, 'sat_phi',rx2sat_azim,...
    'satp_llh',ipp_traj_llh, 'xyzp',rx_traj_des_ipp,'uk_xyzp',uk_xyzp,...
    'vkyz',vkyz,'rngp',rp, 'ipp2rx_east2south_angle',ipp2rx_east2south_angle, 'phip',ipp2rx_to_down_angle,'veff',veff,'A',A,'B',B,'C',C,...
    's',u_mag_field,'thetaB',thetaB,'phiB',phiB,'cosBP',cosBP);
return
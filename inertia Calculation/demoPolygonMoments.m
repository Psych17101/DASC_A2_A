close all

%% Demo 

disp('-------------------------------------')
disp('Demo: Moments of a circle with r=10')
disp('-------------------------------------')
r=10; 
th=(0:1:359)'; 
xy=r*[cos(th*pi/180),sin(th*pi/180)];
Polygon=PolygonMoments(xy,[],-4)
disp 'In this example, because of the 360 degree symmetry,'
disp 'any x-y axes are principal axes'
disp 'hit enter to continue'

pause

% Polygon = 
%                     Area: 314.14
%                      MAx: -5.6029e-14
%                      MAy: -1.7645e-13
%                      Ixx: 7853.2
%                      Iyy: 7853.2
%                      Ixy: -3.2611e-13
%                      Izz: 15706
%          XGirationRadius: 4.9999
%          YGirationRadius: 4.9999
%          ZGirationRadius: 7.0709
%                      ACx: -5.6169e-16
%                      ACy: -1.7836e-16
%                    IxxAC: 7853.2
%                    IyyAC: 7853.2
%                    IxyAC: -3.2611e-13
%                    IzzAC: 15706
%        XACGirationRadius: 4.9999
%        YACGirationRadius: 4.9999
%        ZACGirationRadius: 7.0709
%     PrincAxesRotationDeg: 4.081
%      MaxInertiaPrincAxes: 7853.2
%      MinInertiaPrincAxes: 7853.2
%            xyInPrincAxes: [360x2 double]

%% Demo 
disp('--------------------------------------------------')
disp('Demo: Moments of an elipse parallel to x,y axes')
disp('--------------------------------------------------')
r=10; 
th=(0:1:359)'; 
xy=r*[cos(th*pi/180),3*sin(th*pi/180)];
Polygon=PolygonMoments(xy,[],-4)
disp 'hit enter to continue'
pause
% Polygon = 
% 
%                     Area: 942.43
%                      MAx: 1.9304e-12
%                      MAy: -1.7124e-12
%                      Ixx: 2.1204e+05
%                      Iyy: 23560
%                      Ixy: -1.8829e-12
%                      Izz: 2.356e+05
%          XGirationRadius: 15
%          YGirationRadius: 4.9999
%          ZGirationRadius: 15.811
%                      ACx: -1.817e-15
%                      ACy: 2.0483e-15
%                    IxxAC: 2.1204e+05
%                    IyyAC: 23560
%                    IxyAC: -1.8829e-12
%                    IzzAC: 2.356e+05
%        XACGirationRadius: 15
%        YACGirationRadius: 4.9999
%        ZACGirationRadius: 15.811
%     PrincAxesRotationDeg: 5.724e-16
%      MaxInertiaPrincAxes: 2.1204e+05
%      MinInertiaPrincAxes: 23560
%            xyInPrincAxes: [360x2 double]
%% Demo 

disp('--------------------------------------------------')
disp('Demo: Moments of an elipse rotated 30 degress')
disp('--------------------------------------------------')
r=10; 
rot=30*pi/180;
th=(0:1:359)'; 
xy=r*[cos(th*pi/180),3*sin(th*pi/180)];
xy=([cos(rot) -sin(rot);sin(rot) cos(rot)]*xy')';
Polygon=PolygonMoments(xy,[],-4)
disp 'hit enter to continue'
pause
% Polygon = 
% 
%                     Area: 942.43
%                      MAx: 3.8215e-12
%                      MAy: -2.4964e-12
%                      Ixx: 1.6492e+05
%                      Iyy: 70679
%                      Ixy: -81613
%                      Izz: 2.356e+05
%          XGirationRadius: 13.228
%          YGirationRadius: 8.66
%          ZGirationRadius: 15.811
%                      ACx: -2.6489e-15
%                      ACy: 4.055e-15
%                    IxxAC: 1.6492e+05
%                    IyyAC: 70679
%                    IxyAC: -81613
%                    IzzAC: 2.356e+05
%        XACGirationRadius: 13.228
%        YACGirationRadius: 8.66
%        ZACGirationRadius: 15.811
%     PrincAxesRotationDeg: 30
%      MaxInertiaPrincAxes: 2.1204e+05
%      MinInertiaPrincAxes: 23560
%            xyInPrincAxes: [360x2 double]

%% Demo 

disp('---------------------------------------------------------------')
disp('Demo: Moments of an elipse rotated 30 degress, translated [3,5]')
disp('---------------------------------------------------------------')
r=10; 
rot=30*pi/180;
th=(0:1:359)'; 
xy=r*[cos(th*pi/180),3*sin(th*pi/180)];
xy=([cos(rot) -sin(rot);sin(rot) cos(rot)]*xy')';
xy=xy+repmat([3,5],size(xy,1),1);
Polygon=PolygonMoments(xy,[],-4)
disp 'hit enter to continue'
pause

% Polygon = 
% 
%                     Area: 942.43
%                      MAx: 942.43
%                      MAy: 2827.3
%                      Ixx: 1.6586e+05
%                      Iyy: 79161
%                      Ixy: -78785
%                      Izz: 2.4502e+05
%          XGirationRadius: 13.266
%          YGirationRadius: 9.1649
%          ZGirationRadius: 16.124
%                      ACx: 3
%                      ACy: 1
%                    IxxAC: 1.6492e+05
%                    IyyAC: 70679
%                    IxyAC: -81613
%                    IzzAC: 2.356e+05
%        XACGirationRadius: 13.228
%        YACGirationRadius: 8.66
%        ZACGirationRadius: 15.811
%     PrincAxesRotationDeg: 30
%      MaxInertiaPrincAxes: 2.1204e+05
%      MinInertiaPrincAxes: 23560
%            xyInPrincAxes: [360x2 double]
%% Demo 
% JM
disp('---------------------------------')
disp('Demo: Moments of JM areas')
disp('---------------------------------')
xy=[...
    53.3000   183.4000
    53.3000    36.7000
    68.1000    36.7000
    68.1000   176.4000
    74.4000   181.3000
    80.1000   181.3000
    80.1000    36.7000
    94.9000    36.7000
    94.9000   176.4000
   101.2000   181.3000
   106.9000   181.3000
   106.9000    36.7000
   121.6000    36.7000
   121.6000   181.3000
   113.2000   194.7000
    98.7000   194.7000
    94.9000   191.2000
    85.0000   191.2000
    80.1000   194.7000
    68.1000   194.7000
    53.3000   183.4000
    37.4000   194.7000
    22.2000   194.7000
    22.2000    43.0000
    22.2000    39.5000
    19.1000    36.7000
    16.9000    39.1000
    15.2000    43.0000
     3.5000    43.0000
     3.5000    33.1000
    17.3000    21.8000
    28.6000    21.8000
    37.4000    29.6000
    37.4000   194.7000
];
JMPolygon=PolygonMoments(xy,[0,0],2)

% JMPolygon = 
% 
%                     Area: 10045
%                      MAx: 1.1328e+06
%                      MAy: 7.0763e+05
%                      Ixx: 1.5133e+08
%                      Iyy: 6.0631e+07
%                      Ixy: 8.2061e+07
%                      Izz: 2.1196e+08
%          XGirationRadius: 122.74
%          YGirationRadius: 77.693
%          ZGirationRadius: 145.27
%                      ACx: 70.449
%                      ACy: 112.78
%                    IxxAC: 2.3572e+07
%                    IyyAC: 1.078e+07
%                    IxyAC: 2.2549e+06
%                    IzzAC: 3.4352e+07
%        XACGirationRadius: 48.443
%        YACGirationRadius: 32.76
%        ZACGirationRadius: 58.48
%     PrincAxesRotationDeg: -9.7097
%      MaxInertiaPrincAxes: 2.3958e+07
%      MinInertiaPrincAxes: 1.0394e+07
%                    M_0_0: 10045
%            xyInPrincAxes: [34x2 double]
 disp 'Only gray areas are used for computation of moments'
 disp 'Notice that areas don''t need to be connected'
disp 'hit enter to continue'
pause

%% Demo 
disp('---------------------------------')
disp('Demo: Moments of disjoint areas')
disp('---------------------------------')
origin=[0 0];
Lx=200;
Ly=400;
rect=[  origin
   Lx   origin(2)
   Lx   Ly
   origin(1) Ly
   origin];
   
elipse=[    19.305   111.447
    20.669   128.476
    23.851   136.383
    28.700   151.587
    34.004   157.669
    40.823   165.576
    49.460   172.266
    57.340   178.956
    73.705   188.687
    96.586   191.119
   113.709   191.119
   130.378   186.254
   142.803   181.388
   159.472   168.008
   170.534   152.804
   176.747   142.465
   180.989   127.868
   183.111   112.055
   180.535    96.851
   177.201    83.471
   173.261    74.348
   165.836    64.009
   153.259    52.453
   138.561    42.722
   124.014    34.816
   102.648    33.600
    81.585    35.424
    54.157    46.980
    41.126    58.535
    30.518    72.523
    26.276    80.430
    21.730    94.418
    20.063   100.500
    19.305   111.447];
tongue=[   110.073    72.523
   114.467    80.430
   120.680    90.161
   125.983   100.500
   131.742   104.149
   139.167    95.026
   146.743    83.471
   151.441    79.822
   155.987    79.822
   159.775    87.728
   160.078    97.459
   156.744   107.798
   152.653   111.447
   147.349   114.488
   140.530   115.096
   134.924   116.921
   128.559   122.395
   123.407   124.827
   120.680   124.827
   117.346   118.745
   113.406   106.582
   109.770    92.594
   107.497    80.430
   107.345    74.956
   110.073    72.523];
   
   star=[  119.01   260.45
   134.32   308.50
   121.13   359.59
   145.83   341.34
   162.81   386.35
   162.81   324.31
   186.14   303.03
   162.05   282.96
   159.02   225.18
   144.02   274.44
   119.01   260.45];
   J=[    50.521   310.931
    47.945   296.943
    46.581   281.739
    43.247   242.815
    40.217   230.651
    37.337   228.219
    33.549   223.961
    29.609   228.219
    25.215   229.435
    22.790   232.476
    22.336   239.166
    20.366   254.370
    19.457   266.534
    19.154   275.657
    20.972   276.265
    22.942   275.657
    26.276   273.224
    28.700   271.400
    27.488   260.452
    26.427   253.762
    25.821   248.289
    25.973   238.558
    28.700   240.382
    32.488   244.031
    34.458   264.101
    35.367   285.996
    34.761   302.417
    32.337   312.756
    37.640   310.931
    40.974   307.282
    46.732   307.282
    50.521   310.931];
    J(:,2)=J(:,2)+50;
    elipse(:,2)=elipse(:,2)+30;
    tongue(:,1)=tongue(:,1)-20;
    xy=[rect;
    elipse;origin
    tongue;origin
    star;origin
    J+repmat([20,0],size(J,1),1);origin
    flipud(J)+repmat([200,0],size(J,1),1);origin
    flipud(star)+repmat([100,-200],size(star,1),1)
    ];
    xy(:,1)=5*xy(:,1);

 Polygon=PolygonMoments(xy,[0,2;2,0;1,1],2)
 disp 'Only gray areas are used for computation of moments'
 disp 'Notice that areas don''t need to be connected'
disp 'hit enter to continue'
pause
% Polygon = 
% 
%                     Area: 3.0552e+05
%                      MAx: 6.2423e+07
%                      MAy: 1.6806e+08
%                      Ixx: 1.7603e+10
%                      Iyy: 1.3309e+11
%                      Ixy: 3.3074e+10
%                      Izz: 1.5069e+11
%          XGirationRadius: 240.04
%          YGirationRadius: 660.01
%          ZGirationRadius: 702.31
%                      ACx: 550.07
%                      ACy: 204.32
%                    IxxAC: 4.8489e+09
%                    IyyAC: 4.0645e+10
%                    IxyAC: -1.2634e+09
%                    IzzAC: 4.5494e+10
%        XACGirationRadius: 125.98
%        YACGirationRadius: 364.74
%        ZACGirationRadius: 385.89
%     PrincAxesRotationDeg: -2.0189
%      MaxInertiaPrincAxes: 4.0689e+10
%      MinInertiaPrincAxes: 4.8044e+09
%                    M_0_2: 1.7603e+10
%                    M_2_0: 1.3309e+11
%                    M_1_1: 3.3074e+10
%            xyInPrincAxes: [155x2 double]
 
%% Demo 
disp('---------------------------------')
disp('Demo: Figure #2 in J. Marin paper')
disp('---------------------------------')
 point1=[262.619   406.203];
rectwhole=[ 262.619   406.203
   209.043   449.409
   108.806   449.409
   108.806   556.559
   209.043   556.559
   209.043   449.409
   262.619   406.203
   262.619   606.678
    57.823   606.678
    57.823   406.203
   262.619   406.203];
rect=[631.60   449.409
   631.60   607.79
   474.33   607.79
   474.33   449.409
   631.60   449.409];
pie=[ 406.59   720.73
   406.59   886.64
   362.52   879.72
   325.37   860.71
   296.85   836.52
   271.79   801.09
   256.24   762.20
   248.46   720.73
   406.59   720.73];

   xy=[rectwhole;
   rect;
   point1;
   pie];   
Polygon=PolygonMoments(xy,[],2)
% Polygon = 
% 
%                     Area: 75151
%                      MAx: 4.4303e+07
%                      MAy: 2.544e+07
%                      Ixx: 2.7433e+10
%                      Iyy: 1.0936e+10
%                      Ixy: 1.5135e+10
%                      Izz: 3.8369e+10
%          XGirationRadius: 604.19
%          YGirationRadius: 381.47
%          ZGirationRadius: 714.54
%                      ACx: 338.52
%                      ACy: 589.53
%                    IxxAC: 1.3153e+09
%                    IyyAC: 2.3237e+09
%                    IxyAC: 1.3693e+08
%                    IzzAC: 3.639e+09
%        XACGirationRadius: 132.29
%        YACGirationRadius: 175.84
%        ZACGirationRadius: 220.05
%     PrincAxesRotationDeg: 7.5966
%      MaxInertiaPrincAxes: 2.342e+09
%      MinInertiaPrincAxes: 1.297e+09
%            xyInPrincAxes: [26x2 double]
disp 'Figure #2 in J. Marin''s paper'
disp 'hit enter to continue'
pause

%% Demo 
disp('-----------------------------------------------------------------------')
disp('Demo: Volume and Volume-moments of a polynomial over a rectangular area')
disp('-----------------------------------------------------------------------')
 Amn =[5    10    0    0
      20   -30    0    0
       0     0    0   40 ];   
 xy=[0,0; 10, 0; 10,20;0,20]; 
 [Volume, VolMom_x, VolMom_y, PressureCenter_X, PressureCenter_Y]=VolumeAndMomentsUnderPolynomial(Amn,xy)
% Volume =
% 
%    5.3307e+08
% 
% 
% VolMom_x =
% 
%    8.5298e+09
% 
% 
% VolMom_y =
% 
%    3.9982e+09
% 
% 
% PressureCenter_X =
% 
%        7.5003
% 
% 
% PressureCenter_Y =
% 
%        16.001
% 
% 
% Polygon = 
% 
%                     Area: 200
%                      MAx: 2000
%                      MAy: 1000
%                      Ixx: 26667
%                      Iyy: 6666.7
%                      Ixy: 10000
%                      Izz: 33333
%          XGirationRadius: 11.547
%          YGirationRadius: 5.7735
%          ZGirationRadius: 12.91
%                      ACx: 5
%                      ACy: 10
%                    IxxAC: 6666.7
%                    IyyAC: 1666.7
%                    IxyAC: 0
%                    IzzAC: 8333.3
%        XACGirationRadius: 5.7735
%        YACGirationRadius: 2.8868
%        ZACGirationRadius: 6.455
%     PrincAxesRotationDeg: 0
%      MaxInertiaPrincAxes: 6666.7
%      MinInertiaPrincAxes: 1666.7
%                    M_0_0: 200
%                    M_0_1: 2000
%                    M_0_2: 26667
%                    M_1_0: 1000
%                    M_1_1: 10000
%                    M_1_2: 1.3333e+05
%                    M_2_0: 6666.7
%                    M_2_1: 66667
%                    M_2_3: 1.3333e+07
%                    M_2_4: 2.1333e+08
%                    M_3_3: 100000000
%            xyInPrincAxes: [4x2 double]
disp 'hit enter to continue'
pause
%% Demo
disp('-----------------------------------------------------------------------')
disp('Demo: Moments of an elipse rotated 30 degress translated to point (3,5)')
disp('-----------------------------------------------------------------------')
r=10; 
rot=30*pi/180;
th=(0:1:359)'; 
% xy vertices without the translation:
xy=r*[cos(th*pi/180),3*sin(th*pi/180)];
xy=([cos(rot) -sin(rot);sin(rot) cos(rot)]*xy')';

% requesting the Ixx wrt axis passing through [3,5]
p=0;
q=2;
X0=3;
Y0=5;
% compute necessary moments before calling TranslateMoment
mn=[];
for m=0:p
    for n=0:q
        mn=[mn;m,n];
    end
end
Polygon=PolygonMoments(xy,mn);

% perform the translation
TMpq=TranslateMoment(p,q,X0,Y0,xy,Polygon)

% the same result if Polygon is not precomputed:
TMpq=TranslateMoment(p,q,X0,Y0,xy)


%% 
disp '--------------------------------------------------------------'
disp 'For more application of these formulas read the source paper:'
disp 'Marin, J., 1984,'
disp 'Computing Columns, Footings and Gates through Moments of Area,'
disp 'Computers & Structures, Vol. 18., Nr. 2, pp. 343-349, 1984.'
disp '--------------------------------------------------------------'

% Group difference e.g. CTL and ASD based on Surfstat

clear ;
clc ;

load('C:\Users\Budha\Desktop\ABIDE\Analysis\Codes_Used\Used_in_paper\constraintdata_used.mat') ;
load avsurf_latest ;
load mask_latest ;

age3=age2-mean(age2);      
Age  = term(age3) ;
% Site = term(site) ;
Div  = term(div) ;

M = 1 + Age + Div + Age*Div ;
%M = 1 + Age + Age^2 + Div + Age*Div + (Age^2)*Div ;
%M = 1 + Age + Age^2 + Age^3 + Div + Age*Div + (Age^2)*Div + (Age^3)*Div ;

slm = SurfStatLinMod( Y, M, avsurf ) ;



slm = SurfStatT( slm, age3.*Div.ASD - age3.*Div.CTL ) ;

SurfStatView(SurfStatQ(slm,mask), avsurf)



[ pval, peak, clus ] = SurfStatP( slm, mask );

figure
[ a, cb ] = SurfStatView( pval, avsurf, 'RFT P Vertex' );
    
    
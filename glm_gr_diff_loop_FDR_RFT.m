% Group difference e.g. CTL and ASD based on Surfstat

clear ;
clc ;

load('C:\Users\Budha\Desktop\ABIDE\Analysis\Codes_Used\Used_in_paper\constraintdata_used.mat') ;
load avsurf_latest ;
load mask_latest ;
load('site.mat')

age_center = 6:35;
agep = strread(num2str(age_center), '%s'); 

t_matrix = zeros(length(age_center),81924) ;
q_matrix = zeros(length(age_center),81924) ;

for i = 1:length(age_center);
    
    d1 = [ agep{i}, 'ASD_CTL_RFT'];
    
    age2 = age - age_center(i) ;
    
    Age  = term(age2) ;
    Site = term(site) ;
    Div  = term(div) ;
    
%    M = 1 + Age + Div + Age*Div ;
%    M = 1 + Age + Age^2 + Div + Age*Div + (Age^2)*Div ;
   M = 1 + Age + Site + Age^2 + Age^3 + Div + Age*Div + (Age^2)*Div + (Age^3)*Div ;
    
    slm = SurfStatLinMod( Y, M, avsurf ) ;
    
    
    %% Effect of Diagnostic Division on Cortical Thickness
    contrast = Div.ASD - Div.CTL ;
    
    slm = SurfStatT( slm, contrast ) ;
    
    t_temp = slm.t ;
    
    qval = SurfStatQ(slm, mask) ;
    q_temp = qval.Q ;
    
    
    h=0.39;
    w=0.4;
    agen = age_center(i);
    

    h=figure(1); set(h,'OuterPosition',[600 100 800 800]);
%     [ a, cb ] = SurfStatViewJDL( slm.t.*mask, avsurf, 'T (571 df) for ASD-CTL' ) ;
%    [ a, cb ] = SurfStatViewJDL( SurfStatP( slm, mask ), avsurf, 'ASD-CTL, RFT-corrected' );
%     [ a, cb ] = SurfStatViewJDL(SurfStatQ(slm, mask), avsurf, 'ASD-CTL, FDR-corrected') ;

%     [ pval, peak, clus, clusid ] = SurfStatP( slm, mask ) ;
%     struct.P = pval.P;
%     struct.mask = mask ;
%     [ a, cb ] = SurfStatViewJDL( struct, avsurf, 'ASD-CTL, RFT-corrected' );
    
     [ pval, peak, clus ] = SurfStatP( slm, mask, 0.001 );
%      if agen>27
%         pval.C = pval.P ;
%      end
%      pvalP = rmfield(pval,'C');
    % [ a, cb ] = SurfStatViewJDL2( pvalP, avsurf, 'RFT P Vertex' );
     [ a, cb ] = SurfStatView( pval.P, avsurf, 'RFT P Vertex' );
     
     au=autumn(256); au(1,:)=0.75;
     SurfStatColormap(flipud(au)); SurfStatColLim([0 0.05/30]);  
     
%     % SurfStatColLimJDL( [-4 4] );
%     fprintf('length(a)=%d\n', length(a));
%     a(9)=axes('position',[(0.5-w) 0.07 (w.*2) 0.0001]);
%     fprintf('length(a)=%d\n', length(a));
%     axis([6 35 0 1])
%     set(a(length(a)),'Xgrid','on');
%     set(a(length(a)),'XAxisLocation','top');
%     set(a(length(a)),'TickLength',[0.0 0.0]);
%     set(a(length(a)),'Ygrid','off');
%     set(a(length(a)),'YTick',[]);
%     % set(a(length(a)),'XTick',[-0.0000001,(age-6)/(58-6),1.0000001]);
%     set(a(length(a)),'XTick',[5.8,agen,58.2]);
%     set(a(length(a)),'TickDir','out');
%     set(a(length(a)),'XTickLabel',{'',sprintf('%d',agen),''});
%     set(a(length(a)),'TickLength',[0.02 0.0]);
%     set(a(length(a)),'FontSize',12);
%     set(a(length(a)),'FontWeight','bold');
%     set(a(length(a)),'linewidth',4);
%     
%     agetitle=get(a(length(a)),'Title');
%     set(agetitle,'String','Age','FontSize',12,'FontWeight','bold');
%     agetitpos=get(agetitle,'Position')
%     fprintf('agetitpos=[%d %d %d %d]\n',agetitpos(1), agetitpos(2), agetitpos(3));
%     set(agetitle,'Position',[agetitpos(1), -50, agetitpos(3)])
%     
% %     SurfStatColormap('jet');
% %     cb = SurfStatColLimJDL([-4 4]);
% %     set(cb, 'XTick', [-4 4]);
% %     set(cb, 'XTickLabel', {'-4' '4'});
%     
%     set(cb,'FontSize',8, 'FontWeight', 'Bold');
%     cbtit=get(cb,'Title');
%     set(cbtit,'FontSize',10, 'FontWeight', 'Bold');
%     xl=get(cb,'XLabel');
%     set(xl,'FontWeight','bold');
%     
%     cbpos = get(cb,'OuterPosition');
%     fprintf('cbpos=[%d %d %d %d]\n',cbpos(1), cbpos(2), cbpos(3), cbpos(4));
%     set(cb,'OuterPosition', [cbpos(1) (cbpos(2)+0.075) cbpos(3) cbpos(4)] );
    
    print('-dtiff',d1);
    close

    t_matrix(i,:) = t_temp ;
    q_matrix(i,:) = q_temp ;    
end


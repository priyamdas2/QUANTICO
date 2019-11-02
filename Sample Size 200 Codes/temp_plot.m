QR_Y_sim_TRUEs = csvread('QR_Y_sim_TRUEs.csv');
ts = csvread('QR_ts.csv');
csvwrite('QR_Y_sim_ESTs.csv',QR_Y_sim_ESTs);
QR_Y_sim_ESTs = csvread('QR_Y_sim_ESTs.csv');


for jjj = 1:10
    
    figure
    hold on
    col=hsv(2);
    h = zeros(2,1);
    h(1) = plot(ts, QR_Y_sim_TRUEs(jjj,:),'color',col(1,:));
    h(2) = plot(ts, QR_Y_sim_ESTs(jjj,:),'color',col(2,:));
    axis([0 1 -20 20])
    legend([h(1) h(2)],{'True', 'BEHAVIOR_{QR}'},'Location','northwest');
    title(['True and Estimated Quantile values for Patient number = ' num2str(jjj)])
    xlabel('\tau')
    ylabel('Y')
    hold off
    
end
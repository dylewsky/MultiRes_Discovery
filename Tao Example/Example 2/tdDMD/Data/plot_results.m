close all;
ranks = 4:4:20;
for j = 1:length(ranks)
    r = ranks(j);
    inFile = ['varied_r_' num2str(r) '.fig'];
    openfig(inFile)
    title(['Rank ' num2str(r)])
    if j ~= 1
        ylim([0 6000])
    end
end

%%
close all;
delays = floor(5 * 10.^(0:0.5:2.5));
ymaxes = [3*10^5 3*10^5 1.2*10^5 1.2*10^4 5*10^3 4*10^3];
for j = 1:length(delays)
    inFile = ['varied_delay_' num2str(j) '.fig'];
    openfig(inFile)
    title(['Delay ' num2str(delays(j)) ' steps'])
    ylim([0 ymaxes(j)])
end
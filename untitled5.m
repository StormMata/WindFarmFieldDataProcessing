G = 8.3:0.0001:8.4;

K = 0.003:-0.0001:0;

z = flip(T.Heights)';

R = zeros(length(K),length(G));

index = 6345;

Ref = D.Shear(:,index);

for i = 1:length(K)
    for j = 1:length(G)
        Line = sqrt((G(j) * (1 - exp(-z * sqrt(5.6985e-05/(2 * K(i)))).*cos(z * sqrt(5.6985e-05/(2 * K(i)))))).^2 + ((G(j) * (exp(-z * sqrt(5.6985e-05/(2 * K(i)))) .* sin(z * sqrt(5.6985e-05/(2 * K(i)))))).^2));

        R(i,j) = 1 - sum((Ref - Line).^2)/sum((Ref - mean(Ref)).^2);
    end
end

[M,i] = max(R,[],'all');

[r,c] = ind2sub(size(R,1),i);

figure;
    imagesc(G,K,R);axis xy
    colorbar
    caxis([0 M])

figure;
    plot(Ref,z)
    hold on

z = 43:1:200;

uE = sqrt((Ekman.G(index) * (1 - exp(-z * sqrt(5.6985e-05/(2 * Ekman.K(index)))).*cos(z * sqrt(5.6985e-05/(2 * Ekman.K(index)))))).^2 + ((Ekman.G(index) * (exp(-z * sqrt(5.6985e-05/(2 * Ekman.K(index)))) .* sin(z * sqrt(5.6985e-05/(2 * Ekman.K(index)))))).^2));
uM = sqrt((G(c) * (1 - exp(-z * sqrt(5.6985e-05/(2 * K(r)))).*cos(z * sqrt(5.6985e-05/(2 * K(r)))))).^2 + ((G(c) * (exp(-z * sqrt(5.6985e-05/(2 * K(r)))) .* sin(z * sqrt(5.6985e-05/(2 * K(r)))))).^2));

    plot(uE,z)
    plot(uM,z)

    legend('Profile','Toolbox','Manual')


figure;
    plot(datetime(D.Time,'ConvertFrom','datenum'),Ekman.R)
    ylabel('R^2')
    title('New Ekman Fit')
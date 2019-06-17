function [xk_filt, bordi] = wiener_filter2019(noisy, cycle)

y=noisy+0.0001;

% x1=1; 
% x2=1000;
% y1=1;
% y2=1000;
%y=y(x1:x2, y1:y2);

% Image Size
[M, N]= size(y);

%Assume the Lambda/mean of the noise = 1
ECon = 0.5772156649;
varw = (2 - ECon.^2);

%Homomorphic Filtering
%---------------------------------------------------------------------------------------------

%Apply log to convert multiplicative noise to Additive Noise
ylog= log(y);

%Normalization
[m1,n1]=size(ylog);
meanylog = mean(ylog(:));
 
% %Padarray circular
ylog2=padarray(ylog-meanylog,[M/2 N/2],'circular');

%Power spectrum of the degraded Image
Yf =fftshift( fft2(ylog2));
absYf=abs(Yf);

%Initializing
xmedf=medfilt2(y,[5,5])+1e-3;
xmedlog= zeros(2*m1,2*n1);
xmedlog(m1/2+1:m1/2+m1,n1/2+1:n1/2+n1)=log(xmedf)-mean(log(xmedf(:)));
Xmed=fftshift(fft2(xmedlog));

X0=Xmed;

mu = meanylog + ECon;


beta= 1;

min_alpha=1;
max_alpha=cycle;
cicli=100;  %80
ALPHA=linspace(min_alpha, max_alpha, cicli).^beta; %beta=1.7 %2.7
v=zeros([M,N,8]);
for krt = 1 : 100

    Wk=(abs(X0).^2/numel(X0))./(abs(X0).^2/numel(X0)+varw/2);
    X0=Wk .* Yf;

end

Xlogk =ifft2(ifftshift(X0))+ mu;
xk=abs( exp(Xlogk(m1/2+1:m1/2+m1,n1/2+1:n1/2+n1)));

for krt = 1 : cicli

    Wk=(abs(X0).^2/numel(X0))./(abs(X0).^2/numel(X0)+(ALPHA(krt))*varw/2);
    X1=Wk .* Yf;
    Xlogk =ifft2(ifftshift(X1))+ mu;
    xk=abs( exp(Xlogk(m1/2+1:m1/2+m1,n1/2+1:n1/2+n1)));
  
   
    b=padarray(xk, [1 1],'replicate');
    v(:,:,1)=xk-b(3:end, 2:end-1);     % sotto
    v(:,:,2)=xk-b(1:end-2, 2:end-1);   % sopra
    v(:,:,3)=xk-b(2:end-1, 3:end);     % destra
    v(:,:,4)=xk-b(2:end-1, 1:end-2);   % sinistra
    v(:,:,5)=xk-b(3:end, 1:end-2);     % sotto sinistra
    v(:,:,6)=xk-b(3:end, 3:end);       % sotto destra
    v(:,:,7)=xk-b(1:end-2, 1:end-2);   % sopra sinistra
    v(:,:,8)=xk-b(1:end-2, 3:end);     % sopra destra
    totale=(sum(v.^2,3)./xk);          %iperparametri
    totale=totale./max(totale(:));
        
    xk_stack(:,:, krt)=xk;
    totale_stack(:,:, krt)=totale;
end

% Enhanced wiener Filter
% totale=sum(totale_stack,3);
% totale=totale./max(totale(:));
% bordi=1-totale.^0.1;  %0.1 per squares e homogeneous
% bordi=round(medfilt2((cicli)*bordi));
% bordi(bordi<1)=1;

%%%%%
totale=sum(totale_stack,3);
totale=medfilt2(totale, [5 5]);
totale=totale./max(totale(:));
perc=prctile(totale(:),95);
bordi=(1-totale/perc);  %0.1 per squares e homogeneous
bordi(bordi<=0.01)=0.05;
bordi=round(medfilt2(cicli*bordi)+0.5);

xk_filt=zeros(M,N);

for m=1:M
    for n=1:N
        xk_filt(m,n)=xk_stack(m,n,bordi(m,n));
    end
end

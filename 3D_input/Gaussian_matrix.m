function G=Gaussian_matrix(D,C,kernel)
% by Dr. Xinhuan Zhou from Imperial College London,zhouxinhuan0205@126.com

G=zeros(3*size(D,1),3*size(C,1));
alpha=kernel.alpha;
for i=1:3:3*size(D,1)
    for j=1:3:3*size(C,1)
        r_input=[D(floor(i/3)+1,1) D(floor(i/3)+1,2) D(floor(i/3)+1,3)]-[C(floor(j/3)+1,1) C(floor(j/3)+1,2) C(floor(j/3)+1,3)];
        r2=norm(r_input,2)^2;
        G(i:i+2,j:j+2)=((4*alpha-4*alpha^2*r2)*eye(3,3)+4*alpha^2*r_input'*r_input)*exp(-alpha*r2);
    end
end


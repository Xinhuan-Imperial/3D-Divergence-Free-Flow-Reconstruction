function fun=v_poiseulle(P,cylinder,Q)
% by Xinhuan Zhou, x.zhou16@ic.ac.uk

if any(sqrt((P(:,1)-cylinder.center(1)).^2+(P(:,2)-cylinder.center(2)).^2)>cylinder.radius)||any(P(:,3)>cylinder.length_limit(2))||any(P(:,3)<cylinder.length_limit(1))
    error('Point position out of the pipe');
end

fun=zeros(size(P));
idx=sqrt((P(:,1)-cylinder.center(1)).^2+(P(:,2)-cylinder.center(2)).^2)>cylinder.radius;
fun(idx,:)=0;
r_square=(P(~idx,1)-cylinder.center(1)).^2+(P(~idx,2)-cylinder.center(2)).^2;
w=2*Q/pi/(cylinder.radius)^2*(1-r_square/(cylinder.radius)^2);
fun(~idx,:)=[zeros(sum(~idx),1) zeros(sum(~idx),1) w];
end
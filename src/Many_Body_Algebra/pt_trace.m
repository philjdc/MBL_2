function BT = pt_trace(n,B)
% takes the partial trace over site i and then tensors back in the identity
% to preserve the dimensionality.

I=(diag(n)==1);
J=(diag(n)==0);

BT=zeros(size(B));

BTSmall=(1/2)*(B(I,I)+B(J,J));

BT(I,I)=BTSmall;
BT(J,J)=BTSmall;

end


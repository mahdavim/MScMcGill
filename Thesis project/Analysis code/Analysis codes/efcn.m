function e_ = efcn ( q )

global designmx

PF = (1 ./ ( 1 + exp( -(designmx(:,1:end-1)*q) ) ));

e_ = -sum(log(PF(designmx(:,end)==1))) - sum(log(1-PF(designmx(:,end)==0)));
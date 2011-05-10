function  [Q_pre, obj_pre, A_pre, b_pre, LB_pre, UB_pre, idx_neq_pre, idx_eq_pre, offset, A_back, b_back]= qp_simplify(Q, obj, A, b, LB, UB, idx_neq, idx_eq) ;
% [Q_pre, obj_pre, A_pre, b_pre, LB_pre, UB_pre, idx_neq_pre, idx_eq_pre]= qp_simplify(Q, obj, A, b, LB, UB, idx_neq, idx_eq) ;

A_back=0 ;
b_back=0 ;

num_used=full(sum(A~=0)) ;
idx1s=find(num_used==1) ;

Q_orig=Q ;
Q_diag=full(diag(Q)) ;

obj_orig=obj ;
offset=0 ;

con_keep_map=ones(1, size(A,1)) ;
var_keep_map=ones(1, size(A,2)) ;

Q_new=Q ;
Q_new(idx1s,:)=[] ;
Q_new(:,idx1s)=[] ;
Q_new=full(Q_new) ;

for idx1=idx1s,
    idx1_line=find(A(:,idx1)) ;
    assert(length(idx1_line)==1) ;

    % (a'*x) + c*z = b1    => z = (b1-(a'*x))/c
    % q*z^2/2 + p*z =  q*(b1-(a'*x))^2/c^2/2 + p*(b1-(a'*x))/c = q*b1^2/c^2/2 - q*b1*(a'*x)/c^2/2 + q*(a'*x)^2/c^2/2 + p*b1/c - p*(a'*x)/c 
    %               = (a'*x)^2*(q/c^2)/2 - (a'*x)*(q*b1/c^2/2 + p/c) + (q*b1^2/c^2/2 + p*b1/c)
    
    idx1_var=setdiff(find(A(idx1_line,:)~=0),idx1) ;
    assert(isempty(intersect(idx1_var, idx1s))) ;

    con_keep_map(idx1_line)=0 ;
    var_keep_map(idx1) = 0 ;

    q=Q_diag(idx1) ;
    p=obj(idx1) ;
    c=full(A(idx1_line, idx1)) ;
    a=full(A(idx1_line, idx1_var)) ;
    b1=b(idx1_line) ;

    if ~isempty(idx1_var),
        Q_new(idx1_var,idx1_var) = Q_new(idx1_var,idx1_var) + (a'*a)*q/c^2 ;
        obj(idx1_var) = obj(idx1_var) - a'*(q*b1/c^2/2 + p/c) ;
    end ;
    offset = offset + q*b1^2/c^2/2 + p*b1/c ;
end ;

Q_pre=Q_new(find(var_keep_map), find(var_keep_map)) ;
obj_pre=obj(find(var_keep_map)) ;
A_pre=A(find(con_keep_map), find(var_keep_map)) ;
b_pre=b(find(con_keep_map)) ;

LB_pre=LB(find(var_keep_map)) ;
UB_pre=UB(find(var_keep_map)) ;


idx_eq_pre=[] ;
cum_con_keep_map=cumsum(con_keep_map) ;
for i=idx_eq,
    if con_keep_map(i),
        idx_eq_pre(end+1)=cum_con_keep_map(i) ;
    end ;
end ;

idx_neq_pre=[] ;
cum_con_keep_map=cumsum(con_keep_map) ;
for i=idx_neq,
    if con_keep_map(i),
        idx_neq_pre(end+1)=cum_con_keep_map(i) ;
    end ;
end ;



% fix upper and lower bounds
LB_pre(LB_pre>1e6)=1e6 ;
LB_pre(LB_pre<-1e6)=-1e6 ;
UB_pre(UB_pre>1e6)=1e6 ; 
UB_pre(UB_pre<-1e6)=-1e6 ;

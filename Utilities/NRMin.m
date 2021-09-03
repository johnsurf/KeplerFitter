function [pinc] = NRMin(m, nfree, opzero, hitzero)
%//NRmin Matrix utility used in the minimization procedure

%//Determine Pivots to arrange diagonal elements of Hessian in decreasing order.

%//s vector is local to this function
%vector<double> s;

%//pivot declaration can be moved here for security does not need to be global etc...

pivot = zeros(nfree+1,1);

for i=1:nfree+1
    pivot( i) =  i;
end

hitzero = 0;

% for i  = 1:nfree
%     ii = pivot(i);
%     for j  = i+1:nfree
%         jj = pivot(j);
%         if abs(m(jj,jj)) > abs(m(ii,ii))
%             pivot(j) = ii;
%             pivot(i) = jj;
%             ii = jj;
%         end
%     end
% end

%//Triagonlize Hessian
%//cout << " Check out Symmetry of Hessian " << endl;
%//cout << " m(0)(1) - m(1)(0) = "  << m(0)(1) - m(1)(0) << endl;
%//cout << " m(0)(2) - m(2)(0) = "  << m(0)(2) - m(2)(0) << endl;
%//cout << " m(0)(3) - m(3)(0) = "  << m(0)(3) - m(3)(0) << endl;
%//cout << " m(1)(2) - m(2)(1) = "  << m(1)(2) - m(2)(1) << endl;
%//cout << " m(1)(3) - m(3)(1) = "  << m(1)(3) - m(3)(1) << endl;
%//cout << " m(2)(3) - m(3)(2) = "  << m(2)(3) - m(3)(2) << endl;

%//cout << "Before Triagonalization" << endl;

for i  = 1:nfree
    ii = pivot(i);
    for j  = i+1:nfree
        jj = pivot(j);
        if abs(m(ii,ii)) > 0
            fact = m(jj,ii)/m(ii,ii);
            for k  = i:nfree+1
                kk = pivot(k);
                m(jj,kk) = m(jj,kk) - fact*m(ii,kk);
            end
        end
    end
end

%//cout << "after Triagonalization" << endl;
%//cout << " " << m(0)(0) << ", " << m(0)(1) << ", " << m(0)(2) << ", " << m(0)(3) << ", " << m(0)(4) << endl;
%//cout << " " << m(1)(0) << ", " << m(1)(1) << ", " << m(1)(2) << ", " << m(1)(3) << ", " << m(1)(4) << endl;
%//cout << " " << m(2)(0) << ", " << m(2)(1) << ", " << m(2)(2) << ", " << m(2)(3) << ", " << m(2)(4) << endl;
%//cout << " " << m(3)(0) << ", " << m(3)(1) << ", " << m(3)(2) << ", " << m(3)(3) << ", " << m(3)(4) << endl;

%//Backsubstitution

pinc = zeros(nfree,1);

%//Backsubstitution step of algorithm
%//cout << "Back-Substitute for Solution to Linear System " << endl;
for l  = nfree:-1:1
    ll = pivot(l);
    r = 0.0;
    if ll < nfree
        for j  = l+1:nfree
            jj = pivot(j);
            r = r + m(ll,jj)*pinc(jj);
        end
    end
    %//cout<<"The M matrix"<<endl;
    %//for(int i = 0; i < 4; i++){
    %//      for(int j = 0; j < 5; j++){
    %//      cout<<m(i)(j)<<" ";
    %// }
    %//       cout<<endl;
    %//}
    if abs(m(ll,ll)) < opzero
        disp(["yes 1 and m(ll)(ll) is ",num2str(m(ll,ll))])
        pinc(ll) = 0.0;
        hitzero = 1;
    end
    if abs(m(ll,ll)) >= opzero
         %disp("yes 2")
         pinc(ll) = (m(ll,nfree+1) - r)/m(ll,ll);
    end
end

%//cout<<"pinc 1 OG is : "<<pinc(1)<<endl;

%for( int i = 0; i<nfree; i++){
%    s.push_back( pinc(i) );
%    }
%    return s;
%end
       
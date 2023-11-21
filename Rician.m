% Create the Rician channel matrices
function Hout = Rician(Hlos,FSPL,no_mat,K)                     
    Hlos = repmat(Hlos,no_mat,1);
    Hnlos = sqrt(1/2)*(randn(size(Hlos))+1i*randn(size(Hlos)));
    Htot = FSPL/sqrt(K+1)*(Hlos*sqrt(K)+Hnlos);
    dim = size(Hlos,1)/no_mat;
    for ind = 1:no_mat
       Hout{ind} = Htot((ind-1)*dim+1:ind*dim,:); 
    end
end
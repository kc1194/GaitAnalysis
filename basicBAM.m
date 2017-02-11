function D = basicBAM()
% Given a set of inputs, encode them in a BAM matrix 
	A = zeros(4,15);
	A(1,:) = [1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1];
	A(2,:) = [1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1];
	A(3,:) = [1,1,1,-1,-1,-1,1,1,1,-1,-1,-1,1,1,1];
	A(4,:) = [1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1];

	B = zeros(4,10);
	B(1,:) = [1,1,1,1,-1,-1,-1,-1,1,1];
	B(2,:) = [1,1,1,-1,-1,-1,1,1,1,-1];
	B(3,:) = [1,1,-1,-1,1,1,-1,-1,1,1];
	B(4,:) = [1,-1,1,-1,1,-1,1,-1,1,-1];

	M = encoding(A,B);
	 
	C = [1,-1,1,-1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1];
	D = decoding(C,M);
end

function M = encoding(A,B)
	%calculates encoding of heteroassociative A,B samples
	rowsize = size(A,2);
	columnsize = size(B,2);
	itrns = size(A,1);
	M = zeros(rowsize,columnsize);
	for i = 1:itrns
		M = M + (A(i,:).')*B(i,:);
	end
end

function D = decoding(C,M)
    %returns association of C wrt matrix M
    D = C;
    F = D*M;
    F = normalise(F);
    G = F*(M.');
    G = normalise(G);
    itrns = 1;
    while not(isequal(D,G))
        D = G;
        F = G*M;
        F = normalise(F);
        G = F*(M.');
        G = normalise(G);
        itrns = itrns + 1;
    end
end

function E = normalise(C)
	%normalises C to a binary matrix
	E = C;
	E(E(:)>0) = 1;
	E(E(:)<0) = -1;
end

%* ------------------------------------------------------------------------
%*   Bashar Tahir, bashar.tahir@tuwien.ac.at
%*   (c) 2019 Institute of Telecommunications, TU Wien.
%*   www.nt.tuwien.ac.at 
%*
%*   The algorithm is published under
%*   B. Tahir, S. Schwarz, and M. Rupp, "Constructing Grassmannian Frames by
%*   an Iterative Collision-Based Packing," in IEEE Signal Processing
%*   Letters, 2019. DOI: 10.1109/LSP.2019.2919391 .
%* -----------------------------------------------------------------------
%*   The following function calculates a lower bound on the frame coherence, 
%*   which is described in the publication above.
%*
%%
function bound = getBound(N,M)
    if M <= N^2
            bound = sqrt((M-N)/(N*(M-1)));
    elseif N^2 < M && M <= 2*(N^2-1)
            bound = max(max(sqrt(1/N), sqrt((2*M-N^2-N)/((N+1)*(M-N)))), 1 - 2*M^(-1/(N-1)));
    elseif M > 2*(N^2-1)
            bound = max(sqrt((2*M-N^2-N)/((N+1)*(M-N))), 1 - 2*M^(-1/(N-1)));
    end
end

%</
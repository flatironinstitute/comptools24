function chnkr = squarechunker(n, nref, rect)

    if ( nargin < 2 )
        nref = 0;
    end

    if ( nargin < 3 )
        rect = [-1 1 -1 1];
    end

    m = 2^nref+1;
    [xa, xb, ya, yb] = dealm(rect);
    xab = linspace(xa, xb, m).';
    yab = linspace(ya, yb, m).';

    south = [ xab repmat(ya, m, 1) ];
    east  = [ repmat(xb, m, 1) yab ];
    north = [ flip(xab) repmat(yb, m, 1) ];
    west  = [ repmat(xa, m, 1) flip(yab) ];

    verts = [ south(1:m-1,:) ; east(1:m-1,:) ; north(1:m-1,:) ; west(1:m-1,:) ];
    verts = [ verts ; verts(1,:) ];

    nch = size(verts,1) - 1;

    pref = [];
    pref.k = n;
    chnkr = chunker(pref);
    chnkr = chnkr.addchunk(nch);

    t = legpts(n, [0 1]);

    for k = 1:nch
        chnkr.r(1,:,k) = (1-t)*verts(k,1) + t*verts(k+1,1);
        chnkr.r(2,:,k) = (1-t)*verts(k,2) + t*verts(k+1,2);
        len = sqrt((verts(k,1)-verts(k+1,1)).^2 + (verts(k,2)-verts(k+1,2)).^2);
        vx = (verts(k+1,1)-verts(k,1));
        vy = (verts(k+1,2)-verts(k,2));
        chnkr.d(1,:,k) = vx/2*ones(1, n);
        chnkr.d(2,:,k) = vy/2*ones(1, n);
        chnkr.d2(1,:,k) = zeros(1, n);
        chnkr.d2(2,:,k) = zeros(1, n);
        chnkr.adj(1,k) = k-1;
        chnkr.adj(2,k) = k+1;
    end

    chnkr.adj(1,1) = nch;
    chnkr.adj(2,nch) = 1;

    chnkr.n = normals(chnkr);
    chnkr.wts = weights(chnkr);

    chnkr.nstor(:,:,1:nch) = normals(chnkr);
    chnkr.wtsstor(:,1:nch) = weights(chnkr);

end
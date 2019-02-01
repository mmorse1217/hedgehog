# Pseudocode for adaptive quadrature refinement

The constraints we impose on the collection of patches defining the surface are
as follows:
1. The error incurred from the interpolation of the function defining the parent
patch, $X(u,v)$, it's partial derivatives, and the boundary data $f(u,v)$ at
$p^2$ points to the $p^2$ points of each of it's four children is less than some
absolute tolerance $\eps_a$.

2. the QBKIX points generated from the patch's $p^2$ sample points are positioned in $\Omega$ such that smooth quadrature is accurate to some user-specified tolerance.
  2a. the QBKIX point farthest from the target must be contained in the interior of $\Omega$
  2b. for a $p$th order QBKIX expansion generated from $x \in \partial
    \Omega$, the closest point on $\partial \Omega$ to the $\left \lfloor{p/2}\right \rfloor $-th QBKIX point must be $x$.
  2c. each QBKIX point must be in the far-field of the refined set of patches.

3. the size of each of it's neighboring patches is no more than twice it's own size
  3a. Enforced as a post-processing via p4est 


step 1. is addressed in FaceMap::resolve_function(). It accepts a black-box
evaluator for $f(u,v)$, a set of patches and a threshold tolerance. Currently it does NOT
compute step 1., but rather computes a bidegree p polynomial approximation
$f_p(u,v)$ of $f$ on the patch, and checks that 
$\|f(u,v) - f_p(u,v)\|_\inf/\|f(u,v)\|_\inf \leq \eps_a$. Changing this to the
proper criteria is straightforward. step 3. is a p4est API call. The work is in the
implementation of step 2.

Two entry points into refinement algorithm defining step 2.: refine() and
refine_patches_for_fixed_qbkix_points(). refine() implements steps 2a. and 2b.,
producing a set of patches that serve as the coarse surface discretization.
refine_patches_for_fixed_qbkix_points() implements 2c., which computes the
upsampled patch set to be used to evaluate the solution at the qbkix points.

For the remainder of this document:
..* A patch is a top level patch defined in face-map as a bidegree polynomial of
degree p on [0,1]^2. 
..* A subpatch is a domain restriction of original top level polynomial.
..* A quad is a leaf box of p4est, defining a hierarchical domain restriction of
the top-level patch domain
..* \Omega is the PDE domain

# bdry3d/patch_surf_face_map.cpp
## PatchSurfFaceMap::refine():

# bdry3d/p4est_interface.cpp

## Compute the set of patches defining fixed qbkix points to evaluate on the
## coarse grid

### refine_patches_for_fixed_qbkix_points(p4est, patches):
    refine_patches_point_location(p4est, patches) 
    refine_patches_midpoint_near_medial_axis(p4est, patches)

## Given an initial set of patches and p4est defining connectivity, refine the
## p4est until the corresponding set of new patches are such that their qbkix 
## points generated are all inside \Omega (we only check the last point in the
## expansion i.e. point furthest from the target point).

### refine_patches_point_location(p4est, patches) # needs renaming :)
    // mark all leafs quads in the p4est for refinement
    mark_all_patches_for_refinement(p4est)
    final_patches = None 
    
    while there are patches that need refinement:
        // Use current p4est tree and initial patch set to construct the subpatches
        // on the leaf level; store in final_patches. Generate the last qbkix
        // point in the expansion  for 
        // all subpatches that need refinement, store in qbkix_points
        
        qbkix_points = None
        update_and_resample_face_map(p4est, patches, final_patches, qbkix_points)
        closest_points = mark_target_points(qbkix_points,
                                final_patches)
        for each patch p that needs refinement:
            C_p = closest points to all qbkix points generated from p
            if any of C_p is outside of \Omega:
                break
            else: 
                mark p as valid

        Refine marked p4est quads one additional level

## Given an initial set of patches and p4est defining connectivity, refine the
## p4est until the corresponding set of new patches are such that their qbkix 
## points generated are closest to their 

### refine_patches_midpoint_near_medial_axis(p4est, patches)
    // mark all leafs quads in the p4est for refinement
    mark_all_patches_for_refinement(p4est)
    final_patches = None 
    
    while there are patches that need refinement:
        // Use current p4est tree and initial patch set to construct the subpatches
        // on the leaf level; store in final_patches. Generate the middle 
        // qbkix point for all subpatches that need refinement, 
        // store in qbkix_points
        
        qbkix_points = None
        update_and_resample_face_map(p4est, patches, final_patches, qbkix_points)
        
        // Same functionality as above, but returns all candidate closest points
        // instead of strictly the closest one
        closest_points_map = mark_target_points(qbkix_points,
                                final_patches)

        // Since one qbkix point can be near multiple patches at equal distances
        // (qbkix points on patches edges for example), we need to check all
        // computed closest points on nearby patches
        for each patch p that needs refinement:
            mark p as valid
            
            for each qbkix point q generated from p:
                C_q = all closest points to q
                C_q,min = min( d(c, q)) for c in C_q
                C_q,p = closest point to q on patch p
                
                if C_q,p is None:
                    // if there is no closest point on patch p, we're too close
                    // to the other side of the surface.

                    mark p for refinement
                    break;
                else if C_q,p != C_q,min:
                    // if the absolute closest point is not equal to the closest
                    // point on p, we're too close to some other patch (likely
                    // in a region of high-ish curvature; split
                    
                    mark p for refinement
                    break;
        Refine marked p4est quads one additional level


## Given an initial set of patches and p4est defining connectivity, refine the
## p4est until the corresponding set of new patches are such that the qbkix 
## points generated from the original input patch set are in the far-field of
## the final patch set.

### refine_patches_for_fixed_qbkix_points(p4est, patches)
    // mark all leafs quads in the p4est for refinement
    mark_all_patches_for_refinement(p4est)
    final_patches = None 
    // Using the original p4est tree, which is assumed to define the coarse
    // surface discretization, generate the qbkix points needed to evaluate at
    // the collocation points on the coarse discretization. Throw out the
    // generated patches
    target_qbkix_points = None
    update_and_resample_face_map(p4est, patches,  ____, target_qbkix_points)

    
    while there are patches that need refinement:
        // Use current p4est tree and initial patch set to construct the subpatches
        // on the leaf level; store in final_patches. Destroy generated qbkix
        // points 
        
        update_and_resample_face_map(p4est, patches, final_patches, ____)
        mark_all_patches_as_valid(p4est)     
        
        // in practice this loop is accelerated by a spatial grid queries
        // and near-zone bounding boxes
        for each qbkix point q in target_qbkix_points:
            for each patch p in final_patches:
                if q is in the near-zone of p:
                    mark p for refinement
        Refine marked p4est quads one additional level



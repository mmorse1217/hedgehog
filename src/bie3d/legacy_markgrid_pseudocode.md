marknode(B):
    # compute the nearest surface sample point for each box point of B
    
    #one for each point on an FMM box ( point per face,
    #edge, corner and one in center) we call these "box points"
    #NOTE: domainant := sample is within radius .5 of the center of patch 
    
    for each box point b_i (i=1,...,27):
        Find the nearest dominant sample point s_i on the surface boundary to b_i
        b_i.nearest = s_i
    

     //p = vector pointing from the ith box point of the FMM box to 
     //    closest sample point                                   
     /*                                                             
      *         |                                                   
      *  x ---> o <--- x                                            
      *         |                                                   
      *    <=== |                                                   
      *         |                                                   
      *                                                             
      *  x = possible box point locations                           
      *  o = nearest boundary sample point to x                     
      *  --> = p                                                    
      *  <== = surface normal at o                                  
      *                                                             
      * Below: n = normal at o.                                     
      * If: p \dot n > 0, x is outside the boundary hence outside   
      *     p \dot n = 0, p is perp to n, hence x is closer to the  
      *                   surface at another point assuming geometry 
      *                   is nice. NOTE this is important to        
      *                   consider for complex geometry as an edge case;
      *                   but this will never be exactly zero in practice
      *     p \dot n < 0, x is in the interior                      
      *    ^ depends on convention, logic might need to flipped in the code 
      *      since these seem reversed
      */  

    for each box point b_i:
        s_i = b_i.nearest
        n_i = surface normal at s_i

        p = b_i - s_i
        n = outward normal at s_i
        
        if p \dot n >= 0:
            #b_i is in the exterior
        else
            #b_i is in the interior
    
    if all b_i are in exterior:
        B.color = out
    else if all b_i are in interior:
        B.color = in
    else:
        # mark for further processing; box intersects boundary
        B.color = boundary 


markleafgrid(B):
    # perform a variant of marknode on leaf boxes that intersect the boundary
    k = # samples per dimension
    Gr = k^3 uniform samples in a grid in box B
    #bis3dov.cpp:490
    for g_i in Gr:
        
        Find the nearest dominant sample point s_i on the surface boundary to g_i
        g_i.nearest = s_i
        
        p = g_i - s_i
        n = outward normal at s_i
        while .....:
            TODO
                .......
        if p \dot n >= 0:
            #g_i is in the exterior
            g_i.color = out
        else
            #g_i is in the interior
            g_i.color = in
    

    #### Flood-fill remaining grid points 
    #### This might be buggy, doesnt quite seem right to do it once from each direction
    #### bis3dov.cpp:645
    for g_i in Gr from left->right, bottom->top, back->front:
        #propagate in/out labels forwards: 
        if g_{i+1} is unmarked:
            g_{i+1}.color = g_i.color
    
    for g_i in Gr in reverse order:
        #propagate in/out labels backwards: 
        if g_{i-1} is unmarked:
            g_{i-1}.color = g_i.color

    If all g_i.color are equal:
        B.color = g_i.color
    else:
        B.color = boundary
    
    #Store Gr with B for later use.
    return Gr

# mark each FMM box as in/out/intersecting the boundary
marknodes():
    for each box B in reverse level-order:
        if B is a leaf:
            marknode(B)
    for each leaf box B such that B.color = boundary:
        markleafgrid(B)

markgrid(targets, bs):
    marknodes()
    
    # Mark points that are within bs of the boundary for more careful processing
    for each target point t_i \in targets:
        Find the closest surface sample point s_i to t_i
        if ||t_i - s_i|| < bs
            t_i.nearest_sample = s_i
            t_i.near = true
    
    # Dense computation to mark near surface targets
    # Copied from markleafgrid
    for each target point t_i:
        Find the nearest dominant sample point s_i on the surface boundary to g_i
        g_i.nearest = s_i
        
        p = t_i - s_i
        n = outward normal at s_i
        while .....:
            TODO
                .......
        if p \dot n >= 0:
            #t_i is in the exterior
            t_i.color = out
        else
            #t_i is in the interior
            t_i.color = in
    
    
    for each leaf box B with B.color = boundary:
        for each target t_i \in B:
            Find the nearest grid point g_i to t_i  (see markleafgrid for context)
            t_i.nearest_grid = g_i
            t_i.color = g_i.color
    
    Flood-fill remaining unmarked targets based on marking of containing leaf
    box (computed in markleafgrid). 
    
    return targets


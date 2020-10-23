#ifndef __ON_SURFACE_POINT_HPP__
#define __ON_SURFACE_POINT_HPP__

#include <unistd.h>
#include <common/ebi.hpp>
#include <vec2t.hpp>
BEGIN_EBI_NAMESPACE
enum Region {
    UNMARKED = 0,
    NEAR = 1,
    FAR = 2,
    ON = 3, 
};

enum DomainMembership {
    NOWHERE= 0, // MJM TODO update to scoped enums when moving to c++11 
    OUTSIDE = 1,
    INSIDE = 2,
    OUTSIDE_SPATIAL_GRID = 3, // hacky... sorry future me
    ON_SURFACE = 4 
};

class OnSurfacePoint {
    public:
        // index of patch P containing closest point
        int parent_patch;

        // distance to P(u,v) from ith sample
        double distance_from_target;

        // (u,v) coordinates on patch for closest point
        Point2 parametric_coordinates;
        
        // 3d position of closest point (coord = P(parametric_coordinates) )
        double coord[3];

        // interior normal vector at closest point
        double direction[3];

        // characteristic length of patch with id parent_patch.
        double patch_char_length;

        // Region identifier
        Region region;

        // If point is inside domain
        DomainMembership inside_domain;

        // Target index
        int target_index;

        OnSurfacePoint():
            parent_patch(-1),
            distance_from_target(DBL_MAX), // TODO check this doesn't break things
            parametric_coordinates(Point2(-1,-1)),
            region(UNMARKED),
            inside_domain(NOWHERE), 
            target_index(-1) {;}

        OnSurfacePoint(int parent_patch_,
                double distance_from_target_,
                Point2 parametric_coordinates_,
                Region region_, 
                int target_index_):
            parent_patch(parent_patch_),
            distance_from_target(distance_from_target_),
            parametric_coordinates(parametric_coordinates_),
            region(region_),
            target_index(target_index_) {;}
        
        OnSurfacePoint(const OnSurfacePoint &o){
            parent_patch = o.parent_patch;
            distance_from_target = o.distance_from_target;
            parametric_coordinates = o.parametric_coordinates;
            region = o.region;
            inside_domain = o.inside_domain;
            target_index = o.target_index;
            for(int i=0; i<3; i++)
            {
                coord[i] = o.coord[i];
                direction[i] = o.direction[i];
            }
            patch_char_length = o.patch_char_length;
        }
        OnSurfacePoint& operator =(const OnSurfacePoint &o){
            parent_patch = o.parent_patch;
            distance_from_target = o.distance_from_target;
            parametric_coordinates = o.parametric_coordinates;
            region = o.region;
            inside_domain = o.inside_domain;
            target_index = o.target_index;
            for(int i=0; i<3; i++)
            {
                coord[i] = o.coord[i];
                direction[i] = o.direction[i];
            }
            patch_char_length = o.patch_char_length;
            return *this;
        }


};

inline ostream& operator<<( ostream& os, const OnSurfacePoint& p) {
    cout << "OnSurfacePoint [" << endl;
    cout << "target_index: " << p.target_index<< endl;
    cout << "; parent_patch: " << p.parent_patch << endl;
    cout << "; distance_from_target: " << p.distance_from_target<< endl;
    cout << "; parametric_coordinates: " << p.parametric_coordinates<< endl;
    cout << "; region: " << p.region<< endl;
    cout << "; inside_domain: " << p.inside_domain << endl << "]" << endl;

    return os;
}

END_EBI_NAMESPACE
#endif

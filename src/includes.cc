#include <includes.hh>



OM::Vec2d toVecOM(const Vec2d &vec){
    return {vec[0], vec[1]};
}

Vec2d toVecEigen(const OM::Vec2d vec){
    return {vec[0], vec[1]};
}

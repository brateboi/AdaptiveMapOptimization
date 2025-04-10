#ifndef INCLUDES_HH
#define INCLUDES_HH


#include <iostream>


#include <Eigen/Dense>
//#include <eigen3/Eigen/Dense>


#include <OpenMesh/Core/Mesh/TriMeshT.hh>
#include <OpenMesh/Core/Mesh/DefaultTriMesh.hh>

#include <OpenMesh/Core/IO/MeshIO.hh>

#include <CoMISo/Utils/StopWatch.hh>
#include <CoMISo/Config/config.hh>


#include <CoMISo/NSolver/NPDerivativeChecker.hh>
#include <CoMISo/NSolver/NPTiming.hh>
#include <CoMISo/NSolver/NewtonSolver.hh>
#include <CoMISo/NSolver/TruncatedNewtonPCG.hh>
#include <CoMISo/NSolver/ConstraintTools.hh>
#include <CoMISo/NSolver/LinearConstraint.hh>
#include <CoMISo/NSolver/FiniteElementProblem.hh>
#include <OptimizationElements2D.hh>

#endif // INCLUDES_HH

namespace OM = OpenMesh;
using TM = OM::TriMesh;

using OptimizationTarget = std::vector<std::pair<OM::VertexHandle, OM::Vec3d>>;

using Vec2d = Eigen::Vector2d;
using Mat2d = Eigen::Matrix<double, 2, 2>;


// OM::Vec2d toVec(const Vec2d &vec){
//     return {vec[0], vec[1]};
// }

// Vec2d toVec(const OM::Vec2d vec){
//     return {vec[0], vec[1]};
// }


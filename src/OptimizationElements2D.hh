#ifndef DEFORMATIONELEMENTS2D_HH
#define DEFORMATIONELEMENTS2D_HH


//== INCLUDES =================================================================

#include <fstream>

#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<Eigen/StdVector>
#include<cmath>

#include <CoMISo/NSolver/FiniteElementTinyAD.hh>
#include <CoMISo/NSolver/FiniteElementHessianProjection.hh>
#include <CoMISo/Utils/Polynomials.hh>



//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

class SquaredDistanceElement {
public :


    // define dimensions
    const static int NV = 2; // [vertex_x, vertex_y]
    const static int NC = 2; // [target_x, target_y]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec2T   = Eigen::Matrix<ScalarT,2,1,Eigen::ColMajor>;

        //std::cout<<" - x: "<<_x.transpose()<<std::endl;
        //std::cout<<" - c: "<<_c.transpose()<<std::endl;

        Vec2T d = (_x - _c);

        //std::cout<<" - e = "<<(d.dot(d))<<std::endl;

        return d.dot(d);

    }

    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        return std::numeric_limits<double>::max();
    }

};

using QD_2D_AD = COMISO::FiniteElementTinyAD<SquaredDistanceElement>;
using QD_2D_PH_AD = COMISO::FiniteElementHessianProjection< QD_2D_AD >;


class PullingElement : public SquaredDistanceElement {
public :

    // define dimensions
    const static int NV = 2; // [vertex_x, vertex_y]
    const static int NC = 2; // [direction_x, direction_y]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        VecC xc;
        xc << TinyAD::to_passive(_x[0]), TinyAD::to_passive(_x[1]);

        return SquaredDistanceElement::eval_f(_x, xc + _c);
    }
};


using PULL_2D_AD = COMISO::FiniteElementTinyAD<PullingElement>;
using PULL_2D_PH_AD = COMISO::FiniteElementHessianProjection< PULL_2D_AD >;




// ================================================================================
// ============================= BASE CLASSES
// ================================================================================


class RegularRefElement{

public:

    RegularRefElement(){
        S_inv_ << 2, 1,
            0, std::sqrt(3.0);
        S_inv_ *= 0.5;
        S_inv_ = S_inv_.inverse().eval();
    }


protected:

    Eigen::Matrix2d S_inv_;

};


class JacobianElement2D : public RegularRefElement {
public :

    JacobianElement2D() : RegularRefElement() {}

    template<class ScalarT>
    using JacobianT = Eigen::Matrix<ScalarT, 2, 2, Eigen::ColMajor>;


    template<class ScalarT>
    JacobianT<ScalarT> jacobian(const Eigen::Matrix<ScalarT, 6, 1> &_x) const {
        using VecT = Eigen::Matrix<ScalarT, 2, 1, Eigen::ColMajor>;
        using MatT = Eigen::Matrix<ScalarT, 2, 2, Eigen::ColMajor>;

        Eigen::Map <VecT> A((ScalarT *) &(_x[0]));
        Eigen::Map <VecT> B((ScalarT *) &(_x[2]));
        Eigen::Map <VecT> C((ScalarT *) &(_x[4]));

        //std::cout<<" - S-1: "<<S_inv_<<std::endl;

        MatT J;
        J.col(0) = B - A;
        J.col(1) = C - A;

        return J * S_inv_;
    }

    template<class ScalarT>
    static JacobianT<ScalarT> compute_jacobian(const Eigen::Matrix<ScalarT, 6, 1> &_x) {
        JacobianElement2D elem;
        return elem.jacobian(_x);
    }
};





// ================================================================================
// ============================= SYMMETRIC DIRICHLET WITH BARRIER
// ================================================================================


class SymmetricDirichletWithBarrierElement2D : public JacobianElement2D {
public:

    // define dimensions
    const static int NV = 6; // [ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz] = 2D points A, B, C of a triangle
    //const static int NC = 12; // [weight, barrier_w, barrier_min, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]
    const static int NC = 7; // [weight, barrier_w, barrier_min, E_ref00, E_ref01, E_ref10, E_ref11]

    using VecV = Eigen::Matrix<double, NV, 1>;
    using VecC = Eigen::Matrix<double, NC, 1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const {
        using Vec2T = Eigen::Matrix<ScalarT, 2, 1, Eigen::ColMajor>;
        using Mat2x2T = Eigen::Matrix<ScalarT, 2, 2, Eigen::ColMajor>;
        using Mat2x3T = Eigen::Matrix<ScalarT, 2, 3, Eigen::ColMajor>;
        using Mat3x2T = Eigen::Matrix<ScalarT, 3, 2, Eigen::RowMajor>;
        using Vec2dT = Eigen::Matrix<double, 2, 1, Eigen::ColMajor>;
        using Mat2x2dT = Eigen::Matrix<double, 2, 2, Eigen::ColMajor>;
        using Mat2x3dT = Eigen::Matrix<double, 2, 3, Eigen::RowMajor>;


        // get constants
        const double &w = _c[0];
        const double &barrier_w = _c[1];
        const double &barrier_min = _c[2];


        // get ref edge  matrix
        Mat2x2dT E_ref;

        E_ref.col(0) = _c.segment<2>(3);
        E_ref.col(1) = _c.segment<2>(5);

        // get points of triangle
        Mat3x2T P;
        P <<_x[0], _x[1],_x[2],_x[3],_x[4],_x[5];

        //std::cout<<" P: "<<std::endl<<P<<std::endl;

        // get edge vectors
        Mat2x2T E;
        E.col(0) = P.row(1) - P.row(0);
        E.col(1) = P.row(2) - P.row(0);

        //std::cout<<" E: "<<std::endl<<E<<std::endl;

        // calc  Jacobi Matrix of map
        Mat2x2T J = E * E_ref.inverse();

        //regular ref element

        // calculate determinant
        ScalarT d = J.determinant();
        //std::cout<<" d = "<<d<<std::endl;



        /*std::cout<<" - x: "<<_x.transpose()<<std::endl;
        std::cout<<" - E: "<<std::endl<<E<<std::endl;
        std::cout<<" - E_ref: "<<std::endl<<E_ref<<std::endl;
        std::cout<<" - det = "<<d<<std::endl;
        std::cout<<" - w = "<<w<<std::endl;
        std::cout<<" - barrier w = "<<barrier_w<<std::endl;
        std::cout<<" - barrier min = "<<barrier_min<<std::endl;*/

        if (d <= 0.0) {
            //std::cout << "NEGATIVE DETERMINANT FOUND< INFINITE ENERGY,,     D " << d << std::endl;
            return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
        } else {
            //auto E_sd = (1.0+d)*(J.squaredNorm() + J.inverse().squaredNorm() - 4.0);
            ScalarT E_sd = (J.squaredNorm() + J.inverse().squaredNorm() - 4.0);
            //std::cout<<" E_sd = "<<E_sd<<std::endl;

            //return w * E_sd;

            if (barrier_w && barrier_min - E_sd <= 0.0) {
                std::cout<<" inf barrier term"<<std::endl;
                return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
            }
            //return w * E_sd;

            return w * E_sd + barrier_w / (barrier_min - E_sd);
        }
    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c) {

        // HACK (add exact computation or better line search)
        double t = 2.0;
        for (int i = 0; i < 50; ++i) {
            VecV xn = _x + t * _v;
            //std::cout<<" ------------------------- t"<<i<<" = "<<t<<std::endl;
            if (std::isfinite(eval_f(xn, _c)))
                //if(eval_f(xn, _c) < _c[2]) //alternative, simply not letting the energy reach the maximum
                return t;
            else
                t *= 0.7;
        }

        return 0;
    }


    static int compute_constants(const double _w, const double _barrier_w, const double _barrier_c,
                                 //const Eigen::Vector3d& _p0, const Eigen::Vector3d& _p1, const Eigen::Vector3d& _p2, const Eigen::Vector3d& _p3,
                                 const std::vector <Eigen::Vector2d> &_p,
                                 VecC &_c){

        _c[0] = _w;
        _c[1] = _barrier_w;
        _c[2] = _barrier_c;
        Eigen::Matrix<double ,2,2,Eigen::ColMajor> E;
        E.col(0) = _p[1]-_p[0];
        E.col(1) = _p[2]-_p[0];

        double d  = E.determinant();
        double di = 1.0/d;
        if(!std::isfinite(di) || d < 0.0)
        {
            std::cerr << "ERROR: determinant of edge vectors numerically not valid ---> " << d << std::endl;
            di = 0.0;
            return -1;
        }

        _c.segment<2>(3) = E.col(0);
        _c.segment<2>(5) = E.col(1);

        return 0;

    }
};

using SDEB2D_AD = COMISO::FiniteElementTinyAD<SymmetricDirichletWithBarrierElement2D>;
using SDEB2D_PH_AD = COMISO::FiniteElementHessianProjection< SDEB2D_AD >;


// ================================================================================
// ============================= QUADRATIC DISTANCE PENALTY ELEMENT
// ================================================================================

class QuadraticPenaltyElement2D {
public:
    // define dimensions
    const static int NV = 2; // single vertex, two free variables [x,y]
    const static int NC = 3; // constants, [penalty_weight, target_x, target_y]

    using VecV = Eigen::Matrix<double, NV, 1>;
    using VecC = Eigen::Matrix<double, NC, 1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        const double w = _c[0];
        const double target_x = _c[1];
        const double target_y = _c[2];

        // Compute Squared Distance
        ScalarT dx = _x[0] - target_x;
        ScalarT dy = _x[1] - target_y;

        return w * (dx * dx + dy * dy);
    }

    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        return std::numeric_limits<double>::max();
    }
};

using QP2D_AD = COMISO::FiniteElementTinyAD<QuadraticPenaltyElement2D>;
using QP2D_AD_PH = COMISO::FiniteElementHessianProjection<QP2D_AD>;



// ================================================================================
// ============================= EDGE GLUEING ELEMENT
// ================================================================================


//probably make this generic
class EdgeGlueingElement2D
{
public:

    // define dimensions
    const static int NV = 8; // [] = two 2D points A_left, B_left, A_right, B_right, of two edges to glue
    const static int NC = 1; // [weight]

    using VecV = Eigen::Matrix<double,NV,1>;
    using VecC = Eigen::Matrix<double,NC,1>;

    template<class ScalarT>
    inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
    {
        using Vec2T    = Eigen::Matrix<ScalarT,2,1,Eigen::ColMajor>;

        // get constants
        const double& w         = _c[0];

        Eigen::Map<Vec2T> Al( (ScalarT*) &(_x[0]) );
        Eigen::Map<Vec2T> Bl( (ScalarT*) &(_x[2]) );
        Eigen::Map<Vec2T> Ar( (ScalarT*) &(_x[4]) );
        Eigen::Map<Vec2T> Br( (ScalarT*) &(_x[6]) );

        /*std::cout<<" Al: "<<Al.transpose()<<std::endl;
        std::cout<<" Bl: "<<Bl.transpose()<<std::endl;
        std::cout<<" Ar: "<<Ar.transpose()<<std::endl;
        std::cout<<" Br: "<<Br.transpose()<<std::endl;*/

        Vec2T left = Bl - Al;
        Vec2T right = Br - Ar;

        //std::cout<<" - left edge: "<<left.transpose()<<std::endl;
        //std::cout<<" - right edge: "<<right.transpose()<<std::endl;

        Vec2T diff = left - right;

        //std::cout<<" - diff = "<<diff.transpose()<<std::endl;
        //std::cout<<" - e = "<<(left - right).squaredNorm()<<std::endl;

        ScalarT norm_w = 1.0;//0.5 * (left.norm() + right.norm());

        //return w * (left - right).template lpNorm<1>();
        return (w * (left - right).squaredNorm()) / norm_w;

    }


    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c){
        return std::numeric_limits<double>::max();
    }

};

using EGLUE2D_AD = COMISO::FiniteElementTinyAD<EdgeGlueingElement2D>;
using EGLUE2D_PH_AD = COMISO::FiniteElementHessianProjection< EGLUE2D_AD >;


#endif

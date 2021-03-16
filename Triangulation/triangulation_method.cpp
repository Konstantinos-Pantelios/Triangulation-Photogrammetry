/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>


using namespace easy3d;


/// convert a 3 by 3 matrix of type 'Matrix<double>' to mat3
mat3 to_mat3(Matrix<double> &M) {
    mat3 result;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
            result(i, j) = M(i, j);
    }
    return result;
}


/// convert M of type 'matN' (N can be any positive integer) to type 'Matrix<double>'
template<typename mat>
Matrix<double> to_Matrix(const mat &M) {
    const int num_rows = M.num_rows();
    const int num_cols = M.num_columns();
    Matrix<double> result(num_rows, num_cols);
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j)
            result(i, j) = M(i, j);
    }
    return result;
}



mat34 Rt_selection( mat3 R, vec3 t) {

    mat34 Rt; // 3 by 4
    Rt.set_col(0, R.col(0)); // R is the first column of Rt
    Rt.set_col(1, R.col(1));
    Rt.set_col(2, R.col(2));
    Rt.set_col(3, t);

    return Rt;
}
// ****************************** Triangulation function : returns the 3D Point **************************************//
vec3 triangulation1(vec3 p0, vec3 p1, mat3 K, mat34 Rt, mat34 Rt0){
    mat34   M1 = K*Rt0; // 3 by 4 M // for the second camera
    mat34   M2 = K*Rt;

    Matrix<double> A(4,4,0.0);
    vec4 fr, sr,tr,lr;

    fr = p0.x*M1.row(2) - M1.row(0);
    sr = p0.y*M1.row(2) - M1.row(1);
    tr = p1.x*M2.row(2) - M2.row(0);
    lr = p1.y*M2.row(2) - M2.row(1);
    A.set_row({fr[0],fr[1],fr[2],fr[3]},0);
    A.set_row({sr[0],sr[1],sr[2],sr[3]},1);
    A.set_row({tr[0],tr[1],tr[2],tr[3]},2);
    A.set_row({lr[0],lr[1],lr[2],lr[3]},3);

    Matrix<double> Ua(4, 4, 0.0);   // initialized with 0s
    Matrix<double> Sa(4, 4, 0.0);   // initialized with 0s
    Matrix<double> Va(4, 4, 0.0);   // initialized with 0s

    svd_decompose(A, Ua, Sa, Va);

    vec4 P;
    P = vec4(Va.get_column(Va.cols()-1).data()); //read data from std::vector

    vec3 P_homogeneuous = vec3(P[0]/P[3], P[1]/P[3], P[2]/P[3]);

    return P_homogeneuous;
}

bool Triangulation::triangulation(
        float fx, float fy,     /// input: the focal lengths (same for both cameras)
        float cx, float cy,     /// input: the principal point (same for both cameras)
        const std::vector<vec3> &points_0,    /// input: image points (in homogenous coordinates) in the 1st image.
        const std::vector<vec3> &points_1,    /// input: image points (in homogenous coordinates) in the 2nd image.
        std::vector<vec3> &points_3d,         /// output: reconstructed 3D points
        mat3 &R,   /// output: recovered rotation of 2nd camera (used for updating the viewer and visual inspection)
        vec3 &t    /// output: recovered translation of 2nd camera (used for updating the viewer and visual inspection)

        ) const


{    for (int i = 0; i < points_0.size(); i++){
        if (points_0[i].size() !=points_1.size() && points_0[i].z ==1 && points_1[i].z==1 && points_0[i].size() ==3  && points_1[i].size()==3){
            continue;
        }
        else{
            std::cout << "Input files are not correct. \n";
            return false;
        }
    }

    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for each sub-task. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       triangulation(), or feel free to put them in one or multiple separate files.

    std::cout << "\nTODO: I am going to implement the triangulation() function in the following file:" << std::endl
              << "\t    - triangulation_method.cpp\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tFeel free to use any data structure and function offered by Easy3D, in particular the following two\n"
                 "\tfiles for vectors and matrices:\n"
                 "\t    - easy3d/core/mat.h  Fixed-size matrices and related functions.\n"
                 "\t    - easy3d/core/vec.h  Fixed-size vectors and related functions.\n"
                 "\tFor matrices with unknown sizes (e.g., when handling an unknown number of corresponding points\n"
                 "\tstored in a file, where their sizes can only be known at run time), a dynamic-sized matrix data\n"
                 "\tstructure is necessary. In this case, you can use the templated 'Matrix' class defined in\n"
                 "\t    - Triangulation/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\tPlease refer to the corresponding header files for more details of these data structures.\n\n"
                 "\tIf you choose to implement the non-linear method for triangulation (optional task). Please refer to\n"
                 "\t'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations. \n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (original code framework + your implementation).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without any modification.\n\n" << std::flush;

    /// Easy3D provides fixed-size matrix types, e.g., mat2 (2x2), mat3 (3x3), mat4 (4x4), mat34 (3x4).
    /// To use these matrices, their sizes should be known to you at the compile-time (i.e., when compiling your code).
    /// Once defined, their sizes can NOT be changed.
    /// In 'Triangulation/matrix.h', another templated 'Matrix' type is also provided. This type can have arbitrary
    /// dimensions and their sizes can be specified at run-time (i.e., when executing your program).
    /// Below are a few examples showing some of these data structures and related APIs.

    /// ----------- fixed-size matrices

    /// define a 3 by 4 matrix M (you can also define 3 by 4 matrix similarly)
    mat34 M(1.0f);  /// entries on the diagonal are initialized to be 1 and others to be 0.

    /// set the first row of M
    M.set_row(0, vec4(1,1,1,1));    /// vec4 is a 4D vector.

    /// set the second column of M
    M.set_col(1, vec4(2,2,2,2));

    /// get the 3 rows of M
    vec4 M1 = M.row(0);
    vec4 M2 = M.row(1);
    vec4 M3 = M.row(2);

    /// ----------- fixed-size vectors

    /// how to quickly initialize a std::vector
    std::vector<double> rows = {0, 1, 2, 3,
                                4, 5, 6, 7,
                                8, 9, 10, 11};
    /// get the '2'-th row of M
    const vec4 b = M.row(2);    // it assigns the requested row to a new vector b

    /// get the '1'-th column of M
    const vec3 c = M.col(1);    // it assigns the requested column to a new vector c

    /// modify the element value at row 2 and column 1 (Note the 0-based indices)
    M(2, 1) = b.x;

    /// apply transformation M on a 3D point p (p is a 3D vector)
    vec3 p(222, 444, 333);
    vec3 proj = M * vec4(p, 1.0f);  // use the homogenous coordinates. result is a 3D vector

    /// the length of a vector
    float len = p.length();
    /// the squared length of a vector
    float sqr_len = p.length2();

    /// the dot product of two vectors
    float dot_prod = dot(p, proj);

    /// the cross product of two vectors
    vec3 cross_prod = cross(p, proj);

    /// normalize this vector
    cross_prod.normalize();

    /// a 3 by 3 matrix (all entries are intentionally NOT initialized for efficiency reasons)
    mat3 F;
    /// ... here you compute or initialize F.
    /// compute the inverse of K
    mat3 invF = inverse(F);

    /// ----------- dynamic-size matrices

    /// define a non-fixed size matrix
    Matrix<double> W2(2, 3, 0.0); // all entries initialized to 0.0.

    /// set its first row by a 3D vector (1.1, 2.2, 3.3)
    W2.set_row({ 1.1, 2.2, 3.3 }, 0);   // here "{ 1.1, 2.2, 3.3 }" is of type 'std::vector<double>'

    /// get the last column of a matrix
    std::vector<double> last_column = W2.get_column(W2.cols() - 1);

    // TODO: delete all above demo code in the final submission

//******************************** Normalization of points: image1 and image2 ****************************************//

    float x_1= 0; // image 1
    float y_1 = 0;
    float x_2= 0; // image 2
    float y_2 = 0;
    std::vector<float> x_m1;
    std::vector<float> y_m1;
    std::vector<float> x_m2;
    std::vector<float> y_m2;

    for (vec3 p1: points_0){ //sum of all x and y:  from points_0 dataset
        x_1 = x_1 + p1[0];
        y_1 = y_1 + p1[1];
    }

    for (vec3 p2: points_1){ //sum of all x and y:  from points_1 dataset
        x_2 = x_2 + p2[0];
        y_2 = y_2 + p2[1];
    }

    vec3 centroid1 = {x_1 / points_0.size(),y_1 / points_0.size(),1 }; // centroid computation for the first camera
    vec3 centroid2 = {x_2 / points_1.size(),y_2 / points_1.size(),1 }; // centroid computation for the second camera

    float s1=0;
    float s2=0;

    for (int i=0; i<points_0.size(); i++) { // computation of scale factor for each camera

        float xmean1 = points_0[i].x - centroid1[0];
        float ymean1 = points_0[i].y - centroid1[1];
        x_m1.push_back(xmean1);
        y_m1.push_back(ymean1);
        s1 = s1 + sqrt(x_m1[i]*x_m1[i]+y_m1[i]*y_m1[i]);

        float xmean2 = points_1[i].x - centroid2[0];
        float ymean2 = points_1[i].y - centroid2[1];
        x_m2.push_back(xmean2);
        y_m2.push_back(ymean2);
        s2 = s2 + sqrt(x_m2[i]*x_m2[i]+y_m2[i]*y_m2[i]);
    }
    s1 = s1/points_0.size();
    s2 = s2/ points_1.size();

    s1 = sqrt(2)/(s1); // scale factor for image 1 : 0.00846696
    s2 = sqrt(2)/(s2); // scale factor for image 2 : 0.00983743


    std::vector<vec3> points1; // normalized points image 1 // type: float
    std::vector<vec3> points2; // normalized points image 2 // type: float

    mat3 T1; //Transformation matrix image 1 //type float
    mat3 T2; // Transformation matrix image 2 //type float

    T1.set_row(0,vec3(s1, 0 , -(s1*centroid1[0])));
    T1.set_row(1,vec3(0, s1 , -(s1*centroid1[1])));
    T1.set_row(2,vec3(0, 0 , 1));

    T2.set_row(0,vec3(s2, 0 , -(s2*centroid2[0])));
    T2.set_row(1,vec3(0, s2 , -(s2*centroid2[1])));
    T2.set_row(2,vec3(0, 0 , 1));

//
//    for (int i =0; i<points_0.size(); i++){
//        points1.push_back((points_0[i] - centroid1)*s1);
//        points2.push_back((points_1[i] - centroid2)*s2);
//    }

    for (int i =0; i < points_0.size(); i++){

        points1.push_back(T1.operator*({points_0[i].x,points_0[i].y,points_0[i].z}));
        points2.push_back(T2.operator*({points_1[i].x,points_1[i].y,points_1[i].z}));
    }
 std::cout<< points1<<std::endl;
//********************************************   Foundamental matrix construction ************************************//

    Matrix<double> W1(points_0.size(),9);
    for (int i = 0; i < points_0.size(); i++){
        auto a = points1[i].x  * points2[i].x;
        auto b1 = points1[i].y * points2[i].x;
        auto c1 = points2[i].x;
        auto d = points1[i].x * points2[i].y;
        auto e = points1[i].y * points2[i].y;
        auto f = points2[i].y;
        auto g = points1[i].x;
        auto h = points1[i].y;
        W1.set(i,0,a); W1.set(i,1,b1); W1.set(i,2,c1); W1.set(i,3,d);
        W1.set(i,4,e); W1.set(i,5,f); W1.set(i,6,g); W1.set(i,7,h);
        W1.set(i,8,1);
    }


    // Compute the SVD decomposition of W1

    Matrix<double> U(W1.rows() , W1.rows() , 0.0);   // initialized with 0s
    Matrix<double> S(W1.rows() , W1.cols() , 0.0);   // initialized with 0s
    Matrix<double> V(W1.cols() , W1.cols() , 0.0);   // initialized with 0s

    svd_decompose(W1, U, S, V);

    std::vector<double> f = V.get_column(V.cols()-1);

    Matrix<double> Fq (3,3,0.0);
    Fq.set_row({f[0], f[1], f[2]}, 0);
    Fq.set_row({f[3], f[4], f[5]}, 1);
    Fq.set_row({f[6], f[7], f[8]}, 2);

    Matrix<double> Uf(Fq.rows(), Fq.rows(), 0.0);   // initialized with 0s
    Matrix<double> Sf(Fq.rows(), Fq.cols(), 0.0);   // initialized with 0s
    Matrix<double> Vf(Fq.cols(), Fq.cols(), 0.0);   // initialized with 0s
    svd_decompose(Fq, Uf, Sf, Vf);

    Sf.set(Sf.rows()-1,Sf.cols()-1,0.0); //enforce rank-2 constraint

    auto Ftemp =Uf*Sf*transpose(Vf); // SVD derived from equations USVtranspose

    F = to_mat3(Ftemp);
    F = transpose(T2)*F*T1;   //de-normalization
    F = F/(F(2,2));  //set the last column's element =1 : second constraint

    std::cout << "Foundamental matrix F = " << std::endl << F << std::endl;

    //matrix K same for both images
    mat3 K;
    K.set_row (0, vec3(fx, 0, cx));
    K.set_row(1, vec3(0, fy,cy));
    K.set_row(2, vec3(0, 0, 1));

    mat3 Etemp;
    Etemp = transpose(K)*F*K; // transpose K * E * K

    Matrix<double> Ue(3, 3, 0.0);   // initialized with 0s
    Matrix<double> Se(3, 3, 0.0);   // initialized with 0s
    Matrix<double> Ve(3, 3, 0.0);   // initialized with 0s

    auto E = to_Matrix(Etemp);
    svd_decompose(E, Ue, Se, Ve );

    std::cout << "Essential matrix E = " << std::endl <<  E << std::endl;

// ****************************************** PART 2 retrieve R and t *********************************************** //
    mat3 W;
    W.set_row(0,vec3(0,-1,0));
    W.set_row(1,vec3(1,0,0));
    W.set_row(2,vec3(0,0,1));
//
//    mat3 Z;
//    Z.set_row(0,vec3(0,1,0));
//    Z.set_row(1,vec3(-1,0,0));
//    Z.set_row(2,vec3(0,0,0));

    std::vector<double> t_temp;
    auto t1 = t;
    auto t2 = t;
    auto R1 = R;
    auto R2 = R;

    t_temp = Ue.get_column(Ue.cols()-1); //vec3 float
    t1 = vec3(t_temp.data()); // read it from the std::vector // +u3: positive value of vector t
    t2 = - vec3(t_temp.data()); // read it from the std::vector // -u3: negative value of vector t

    R1 = to_mat3(Ue)*W*transpose(to_mat3(Ve))*determinant(to_Matrix(to_mat3(Ue)*W*transpose(to_mat3(Ve)))); //mat3 float : first case where R = det(UWVT)UWVT
    R2 =  to_mat3(Ue)*transpose(W)*transpose(to_mat3(Ve))*determinant(to_Matrix(to_mat3(Ue)*transpose(W)*transpose(to_mat3(Ve)))); //mat3 float: second case where R = det(UWVTT)UWTVT

    std::vector<vec3> casest = {t1,t2,t1,t2};
    std::vector<mat3> casesR = {R1,R1,R2,R2};
    auto thers1 = determinant(R1);
    auto thers2 = determinant(R2);

    mat34 Rt0;     //set M for the origin R=I t=0
    Rt0.load_identity(1);
    vec3 t3(0,0,0);
    Rt0.set_col(3,t3);

    auto points_3d_1 = points_3d;
    auto points_3d_2 = points_3d;
    auto points_3d_3 = points_3d;
    auto points_3d_4 = points_3d;

//  *************************************************** Triangulation of points **************************************//
    for (int i = 0; i < points_0.size(); i++) {
        auto Rt1 = Rt_selection(R1,t1);
        points_3d_1.push_back(triangulation1(points_0[i], points_1[i], K, Rt1, Rt0));
    }
    for (int i = 0; i < points_0.size(); i++) {
        auto Rt2 = Rt_selection(R1,t2);
        points_3d_2.push_back(triangulation1(points_0[i], points_1[i], K, Rt2, Rt0));
    }
    for (int i = 0; i < points_0.size(); i++) {
        auto Rt3 = Rt_selection(R2,t1);
        points_3d_3.push_back(triangulation1(points_0[i], points_1[i], K, Rt3, Rt0));
    }
    for (int i = 0; i < points_0.size(); i++) {
        auto Rt4 = Rt_selection(R2,t2);
        points_3d_4.push_back(triangulation1(points_0[i], points_1[i], K, Rt4, Rt0));
    }

//  *********************************** Estimation of relative pose of the views *************************************//
    int positivez1=0;
    int positivez2=0;
    int positivez3=0;
    int positivez4=0;
    std::vector<std::vector<vec3>> points;
    points.push_back(points_3d_1);
    points.push_back(points_3d_2);
    points.push_back(points_3d_3);
    points.push_back(points_3d_4);

    for(int i=0; i<points_3d_1.size(); i++){ // check for points with positive z values
        vec3 Rtp=R1*points_3d_1[i] + t1;
        if (points_3d_1[i].z > 0 and Rtp.z > 0){
            positivez1++ ;
        }
    }
    for(int i=0; i<points_3d_2.size(); i++){ // check for points with positive z values
        vec3 Rtp=R1*points_3d_2[i] + t2;
        if (points_3d_2[i].z > 0 and Rtp.z > 0){
            positivez2++ ;
        }
    }
    for(int i=0; i<points_3d_3.size(); i++){ // check for points with positive z values
        vec3 Rtp=R2*points_3d_3[i] + t1;
        if (points_3d_3[i].z > 0 and Rtp.z > 0){
            positivez3++ ;
        }
    }
    for(int i=0; i<points_3d_4.size(); i++){ // check for points with positive z values
        vec3 Rtp=R2*points_3d_4[i] + t2;
        if (points_3d_4[i].z > 0 and Rtp.z > 0){
            positivez4++ ;
        }
    }

    std::vector<int> z_coords = {positivez1,positivez2,positivez3,positivez4};
    int max = z_coords[0];
    int counter0=0;
    int counter1=0;
    for (auto const& i : z_coords){
        if (i >= max){
            max = i;
            counter1 = counter0;
        }
        counter0++;
    }

    points_3d = points[counter1];
    R = casesR[counter1];
    t = casest[counter1];
    std::vector<std::string> cases = {"R = det((U W V^T)(U W V^T)), t = +u3","R = det((U W V^T)(U W V^T)), t = -u3","R = det((U W^T V^T)(U W^T V^T)), t = +u3","R = det((U W^T V^T)(U W^T V^T)), t = -u3"};
    std::cout<< "Rotation matrix R "<< R << " "<< " Determinant of matrix R" << std::endl;
    std::cout<< " Vector t " << t << std::endl;
    std::cout<<" Points with positive z values" << " "<< z_coords[counter1] << std::endl;
    std::cout<< " Case number " << counter1+1 << " -> "<< cases[counter1]<< std::endl;

    // TODO: check if the input is valid (always good because you never known how others will call your function).


    // TODO: Don't forget to
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       However, there are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.
    return points_3d.size() > 0;
}
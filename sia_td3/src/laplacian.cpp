
#include "laplacian.h"
#include "mesh.h"
#include <Eigen/SparseCholesky>

using namespace Eigen;
using namespace pmp;

typedef SparseMatrix<float> SpMat;
typedef PermutationMatrix<Dynamic> Permutation;
typedef Eigen::Triplet<double> T;

double cotan_weight(const SurfaceMesh &mesh, pmp::Halfedge he) {
  auto points = mesh.get_vertex_property<Point>("v:point");

  // TODO

  return 1;
}

/// Computes the Laplacian matrix in matrix \a L using cotangent weights or the
/// graph Laplacian if useCotWeights==false.
void create_laplacian_matrix(const SurfaceMesh &mesh, SpMat &L,
                             bool useCotWeights) {
  // number of vertices in mesh
  int n = (int)mesh.n_vertices();

  // TODO
  std::vector<T> tripletList;
  int diag[n];
  float weight;
  for(int i = 0; i<n ; i++){
    diag[i] = 0;
  }
  auto points = mesh.get_vertex_property<Point>("v:point");

  for(auto vi : mesh.vertices()){
    //std::cout << vi.idx() << "\n";
    int index = 0;
    mesh.vertices(vi).begin();
    for(SurfaceMesh::VertexAroundVertexCirculator vj = mesh.vertices(vi).begin(); vj !=mesh.vertices(vi).end(); ++vj){
      Vertex vertex1 = *(--vj);
      ++vj;
      Vertex vertex2 = *(++vj);
      --vj;

      Eigen::Vector3f a = static_cast<Eigen::Vector3f>(points[vertex1]);
      Eigen::Vector3f b = static_cast<Eigen::Vector3f>(points[vertex2]);
      Eigen::Vector3f j = static_cast<Eigen::Vector3f>(points[*vj]);
      Eigen::Vector3f i = static_cast<Eigen::Vector3f>(points[vi]);

      Eigen::Vector3f v1 = i-a;
      Eigen::Vector3f v2 = j-a;
      Eigen::Vector3f v3 = i-b;
      Eigen::Vector3f v4 = j-b;

      v1 = v1.normalized();
      v2 = v2.normalized();
      v3 = v3.normalized();
      v4 = v4.normalized();

      /*std::cout << "V1:" << v1 << "\n\n";
      std::cout << "V2:" << v2 << "\n\n";
      std::cout << "V3:" << v3 << "\n\n";
      std::cout << "V4:" << v4 << "\n\n";*/

      std::cout << vi.idx() << " " << (*vj).idx() <<" "<<vertex1.idx()<<" " << vertex2.idx() <<"\n";

      double sin = v1.cross(v2).norm();

      //std::cout << "Sin:" << v1.cross(v2).norm() << "\nNorm :" << v1.norm() * v2.norm() << "\n\n";

      double cos = v1.dot(v2);
      double cota = cos / sin;

      sin = v3.cross(v4).norm() / (v3.norm() * v4.norm() );
      cos = v3.dot(v4);
      double cotb = cos / sin;

      
      
      //std::cout << vi.idx() << " " << (*vj).idx() << " "  << vertex1.idx()<< " "  << vertex2.idx() <<"\n";

      if(!useCotWeights){
        weight = 1;
      }else{
        weight = (cota-cotb)/2;
      }

      tripletList.push_back(T(vi.idx(),(*vj).idx(),weight));
      //tripletList.push_back(T((*vj).idx(),vi.idx(),weight));
      diag[vi.idx()]-= weight;
      //diag[(*vj).idx()]-= weight;
      
    }
  }

  for(int i = 0; i<n ; i++){
    tripletList.push_back(T(i,i,diag[i]));
  }
L.setFromTriplets(tripletList.begin(), tripletList.end());
}

/// Computes the permutation putting selected vertices (mask==1) first, and the
/// others at the end. It returns the number of selected vertices.
int create_permutation(const SurfaceMesh &mesh, Permutation &perm) {
  auto masks = mesh.get_vertex_property<int>("v:mask");

  // number of vertices in mesh
  int n = (int)mesh.n_vertices();
  // TODO
  perm.resize(n);

  int i = 0;
  int u = 1;
  for(auto vi : mesh.vertices()){
    if(masks[vi] == 1){
      perm.indices()[vi.idx()] = i;
      i++;
    }else{
      perm.indices()[vi.idx()] = n-u;
      u++;
    }
  }

  return i;
}

/// Performs the poly-harmonic interpolation (order k) over the selected
/// vertices (mask==1) of the vertex attributes u. For each vertex V of index i,
///     if  mask[V]!=1 then u.col(i) is used as input constraints,
///     otherwise, mask[V]==1, and u.col(i) is replaced by the poly-harmonic
///     interpolation of the fixed values.
void poly_harmonic_interpolation(const SurfaceMesh &mesh,
                                 Ref<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > u,
                                 int k) {
  // Number of vertices in the mesh
  int n = (int)mesh.n_vertices();

  // 1 - Create the sparse Laplacian matrix
  SpMat L(n,n);
  create_laplacian_matrix(mesh, L, true);

  // 2 - Create the permutation matrix putting the fixed values at the end,
  //     and the true unknown at the beginning
  Permutation perm;
  int nb_unknowns = create_permutation(mesh, perm);

  // 3 - Apply the permutation to both rows (equations) and columns (unknowns),
  //     i.e., L = P * L * P^-1

  // TODO
  
  if(k == 2){
    L=L*L;
  }
  if(k == 3){
    L=L*L*L;
  }

  L = perm * L * perm.inverse();

  // 4 - solve L * [x^T u^T]^T = 0, i.e., L00 x = - L01 * u

  // TODO
  u = perm * u;
  MatrixXf Ubar = u.bottomRows(n-nb_unknowns) ;
  SpMat L00 = L.topLeftCorner(nb_unknowns,nb_unknowns);
  SpMat L01 = L.topRightCorner(nb_unknowns,n-nb_unknowns);

  MatrixXf b = - L01 * Ubar;

  SimplicialLDLT<SpMat> solver;
  solver.compute(L00);
  
  if(solver.info()!=Success) {
    // decomposition failed
    return; 
  }
  MatrixXf x = solver.solve(b);
  if(solver.info()!=Success) {
    // solving failed
    return;
  }
  
  // 5 - Copy back the results to u

  // TODO
  u <<x,Ubar;
  u = perm.inverse() * u;
}

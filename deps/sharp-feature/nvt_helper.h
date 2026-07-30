[[nodiscard]] Eigen::Matrix3d computeTensor(const Vector& n)
{
    // Convert to Eigen
    Eigen::Vector3d ne(static_cast<double>(n.x()),
                       static_cast<double>(n.y()),
                       static_cast<double>(n.z()));

    // Outer product n * n^T
    Eigen::Matrix3d outer = ne * ne.transpose();

    return outer;
}
double solveTensor(const Eigen::Matrix3d& A)
{
    // Self-adjoint eigen solver
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(A);
    std::array<double,3> l;
    if (solver.info() != Eigen::Success) {
        l = {0.0, 0.0, 0.0};
    }
    // Eigenvalues and eigenvectors
    Eigen::Vector3d eigvals = solver.eigenvalues();       // ascending order
    Eigen::Matrix3d eigvecs = solver.eigenvectors();     // columns = eigenvectors

    // Sort descending (lambda1 >= lambda2 >= lambda3)
    std::array<int,3> idx = {0,1,2};
    std::sort(idx.begin(), idx.end(),
              [&](int a, int b){ return eigvals[a] > eigvals[b]; });

    // Assign sorted eigenvalues
    l[0] = eigvals[idx[0]];
    l[1] = eigvals[idx[1]];
    l[2] = eigvals[idx[2]];
    for (int k = 0; k < 3; ++k)
        l[k] = std::max(0.0, l[k]);

    // Normalize eigenvalues
    double sum = l[0] + l[1] + l[2] + 1e-12;
    l[0] /= sum;
    l[1] /= sum;
    l[2] /= sum;

    auto metric = 2*l[2] + l[1];
    if(metric < 0)
    {
        std::cout<<l[0]<<" "<<l[1]<<" "<<l[2]<<"\n";
    }
    return metric;
}

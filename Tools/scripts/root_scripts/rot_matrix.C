



void rot_matrix(){

    // Unit vector for target to detector direction in det coordinates
    double v1 = 0.46237185;
    double v2 = 0.048854068;
    double v3 = 0.88533923;

    // Origin of det coordinates
    double z1 = 0.0;
    double z2 = 0.0;
    double z3 = 1.0;

    double ux = v2;
    double uy = -1.0*v1;
    double uz = 0.0;

    // cos theta
    double ctheta = v3;
    double t = 1 - ctheta;

    // Sine theta
    double stheta = 0.46494564;


    TRotation Rot;             // Rotations
    std::vector<double> rotmatrix;     // Inputs

    // From beam to detector rotation matrix
    rotmatrix = {
        t*ux*ux+ctheta,    t*ux*uy-stheta*uz, t*ux*uz+stheta*uy,
        t*ux*uy+stheta*uz, t*uy*uy+ctheta,    t*uy*uz-stheta*ux,
        t*ux*uz-stheta*uy, t*uy*uz+stheta*ux, t*uz*uz+ctheta
     };

    for (unsigned int k =0; k < rotmatrix.size(); k++){
        std::cout << rotmatrix.at(k) << std::endl;
    }

    // Return the TRotation
    TVector3 newX, newY, newZ;
    newX = TVector3(rotmatrix[0], rotmatrix[1], rotmatrix[2]);
    newY = TVector3(rotmatrix[3], rotmatrix[4], rotmatrix[5]);
    newZ = TVector3(rotmatrix[6], rotmatrix[7], rotmatrix[8]);

    Rot.RotateAxes(newX, newY, newZ); // Return the TRotation now det to beam
    Rot.Invert(); // Invert back to the beam to det

    double p1 = 0.4036 ;
    double p2 = 0.0440307;
    double p3 = 0.772186;

    double px = rotmatrix[0]*p1 + rotmatrix[1]*p2 + rotmatrix[2]*p3;
    double py = rotmatrix[3]*p1 + rotmatrix[4]*p2 + rotmatrix[5]*p3;
    double pz = rotmatrix[6]*p1 + rotmatrix[7]*p2 + rotmatrix[8]*p3;

    double p = std::sqrt(px*px + py*py + pz*pz);
    double theta = std::acos(pz/p) * 180 / 3.1415;

    std::cout << theta << std::endl; 





}
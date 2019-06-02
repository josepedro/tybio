/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file
 * Implementation of a stationary, pressure-driven 2D channel flow, and
 * comparison with the analytical Poiseuille profile. The velocity is initialized
 * to zero, and converges only slowly to the expected parabola. This application
 * illustrates a full production cycle in a CFD application, ranging from
 * the creation of a geometry and definition of boundary conditions over the
 * program execution to the evaluation of results and production of instantaneous
 * graphical snapshots. From a technical standpoint, this showcase is not
 * trivial: it implements for example hypbrid velocity/pressure boundaries, 
 * and uses an analytical profile to set up the boundary and initial conditions,
 * and to compute the error. As a first Palabos example, you might prefer to 
 * look at a more straightforward code, such as cavity2d.
 **/
 
#include "palabos2D.h"
#include "palabos2D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9Descriptor

/// Velocity on the parabolic Poiseuille profile
T poiseuilleVelocity(plint iY, IncomprFlowParam<T> const& parameters) {
    T y = (T)iY / parameters.getResolution();
    return 4.*parameters.getLatticeU() * (y-y*y);
}

/// Linearly decreasing pressure profile
T poiseuillePressure(plint iX, IncomprFlowParam<T> const& parameters) {
    T Lx = parameters.getNx()-1;
    T Ly = parameters.getNy()-1;
    return 8.*parameters.getLatticeNu()*parameters.getLatticeU() / (Ly*Ly) * (Lx/(T)2-(T)iX);
}

/// Convert pressure to density according to ideal gas law
T poiseuilleDensity(plint iX, IncomprFlowParam<T> const& parameters) {
    return poiseuillePressure(iX,parameters)*DESCRIPTOR<T>::invCs2 + (T)1;
}

/// A functional, used to initialize the velocity for the boundary conditions
template<typename T>
class PoiseuilleVelocity {
public:
    PoiseuilleVelocity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, Array<T,2>& u) const {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }
private:
    IncomprFlowParam<T> parameters;
};

/// A functional, used to initialize the density for the boundary conditions
template<typename T>
class PoiseuilleDensity {
public:
    PoiseuilleDensity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    T operator()(plint iX, plint iY) const {
        return poiseuilleDensity(iX,parameters);
    }
private:
    IncomprFlowParam<T> parameters;
};

/// A functional, used to create an initial condition for with zero velocity,
///   and linearly decreasing pressure.
template<typename T>
class PoiseuilleDensityAndZeroVelocity {
public:
    PoiseuilleDensityAndZeroVelocity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, T& rho, Array<T,2>& u) const {
        rho = poiseuilleDensity(iX,parameters);
        u[0] = T();
        u[1] = T();
    }
private:
    IncomprFlowParam<T> parameters;
};

enum InletOutletT {pressure, velocity};

//Creating a function to draw the star_electrode
void createStar(plint sizeStar, Array<plint, 2> positionStar, MultiBlockLattice2D<T,DESCRIPTOR>& lattice) {
    //Defining Central point of the star
    plint initX = positionStar[0] - sizeStar/2;
    plint initY = positionStar[1] - sizeStar/2;

    //Defining Block 1 positions
    plint x0 = 0.4*sizeStar + initX;
    plint x1 = 0.8*sizeStar + initX;
    plint y0 = 0.4*sizeStar + initY;
    plint y1 = 0.8*sizeStar + initY;
    Box2D block1(x0, x1, y0, y1);
    defineDynamics(lattice, block1, new BounceBack<T,DESCRIPTOR>((T)-99999991));

   //Defining Block 2 positions
    plint halfBlock1 = x1 - x0;
    x0 = x1;
    x1 = x0 + halfBlock1/2;
    y0 = y1 - halfBlock1/2;
    Box2D block2(x0, x1, y0, y1);
    //std::cout << x0 << " " << x1 << " " << y0 << " " << y1 << " " << std::endl;
    defineDynamics(lattice, block2, new BounceBack<T,DESCRIPTOR>((T)-99999991));    
    //Defining Block 3 positions
    x1 = x0;
    x0 = x0 - halfBlock1/2;
    y0 = y0 - halfBlock1;
    y1 = y0 + halfBlock1/2;
    Box2D block3(x0 - halfBlock1/6, x1, y0, y1);
    //std::cout << x0 << " " << x1 << " " << y0 << " " << y1 << " " << std::endl;
    defineDynamics(lattice, block3, new BounceBack<T,DESCRIPTOR>((T)-99999991));
    //Defining Block 4 positions
    x1 = 0.4*sizeStar + initX;
    x0 = x1 - halfBlock1/2;
    y0 = 0.4*sizeStar + initY;
    y1 = y0 + halfBlock1/2;
    Box2D block4(x0, x1, y0, y1);
    //std::cout << x0 << " " << x1 << " " << y0 << " " << y1 << " " << std::endl;
    defineDynamics(lattice, block4, new BounceBack<T,DESCRIPTOR>((T)-99999991));

    //Defining Block 5 positions
    x0 = 0.4*sizeStar + initX;
    x1 = x0 + halfBlock1/2;
    y0 = 0.8*sizeStar + initY;
    y1 = y0 + halfBlock1/2;
    Box2D block5(x0, x1 + halfBlock1/6, y0, y1);
    defineDynamics(lattice, block5, new BounceBack<T,DESCRIPTOR>((T)-99999991)); 
}

void channelSetup( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                   IncomprFlowParam<T> const& parameters,
                   OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition,
                   InletOutletT inletOutlet )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();

    pcout << "Nx" << endl;
    pcout << nx << endl;

    pcout << "Ny" << endl;
    pcout << ny << endl;
    // Note: The following approach illustrated here works only with boun-
    //   daries which are located on the outmost cells of the lattice. For
    //   boundaries inside the lattice, you need to use the version of
    //   "setVelocityConditionOnBlockBoundaries" which takes two Box2D
    //   arguments.

    // Velocity boundary condition on bottom wall. 
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, nx-1, 0, 0) );
    // Velocity boundary condition on top wall. 
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, nx-1, ny-1, ny-1) );

    // Pressure resp. velocity boundary condition on the inlet and outlet.
    if (inletOutlet == pressure) {
        // Note: pressure boundary conditions are currently implemented
        //   only for edges of the boundary, and not for corner nodes.
        boundaryCondition.setPressureConditionOnBlockBoundaries (
                lattice, Box2D(0,0, 1,ny-2) );
        boundaryCondition.setPressureConditionOnBlockBoundaries (
                lattice, Box2D(nx-1,nx-1, 1,ny-2) );
    }
    else {
        boundaryCondition.setVelocityConditionOnBlockBoundaries (
                lattice, Box2D(0,0, 1,ny-2) );
        boundaryCondition.setVelocityConditionOnBlockBoundaries (
                lattice, Box2D(nx-1,nx-1, 1,ny-2) );
    }

    // Define the value of the imposed density on all nodes which have previously been
    //   defined to be pressure boundary nodes.
    setBoundaryDensity (
            lattice, lattice.getBoundingBox(),
            PoiseuilleDensity<T>(parameters) );
    // Define the value of the imposed velocity on all nodes which have previously been
    //   defined to be velocity boundary nodes.
    setBoundaryVelocity (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocity<T>(parameters) );
    // Initialize all cells at an equilibrium distribution, with a velocity and density
    //   value of the analytical Poiseuille solution.
    initializeAtEquilibrium (
           lattice, lattice.getBoundingBox(),
           PoiseuilleDensityAndZeroVelocity<T>(parameters) );

    
    //Inserting Periodic Boundary Condition
    //lattice.periodicity().toggleAll(true); // Use periodic boundaries.
  
    //Inserting BounceBack Boundary Condition
    //Inserting Star BounceBack 
    plint sizeStar = 50;
    //First Electrode line
    createStar(sizeStar, Array<plint, 2>(175, 105), lattice);
    createStar(sizeStar, Array<plint, 2>(235, 87), lattice);
    createStar(sizeStar, Array<plint, 2>(296, 105), lattice);
    createStar(sizeStar, Array<plint, 2>(355, 87), lattice);
    createStar(sizeStar, Array<plint, 2>(415, 97), lattice);
    createStar(sizeStar, Array<plint, 2>(475, 96), lattice);
    //Second Electrode line
    createStar(sizeStar, Array<plint, 2>(175, 247), lattice);
    createStar(sizeStar, Array<plint, 2>(246, 256.5), lattice);
    createStar(sizeStar, Array<plint, 2>(296, 256.5), lattice);
    createStar(sizeStar, Array<plint, 2>(355, 256), lattice);
    createStar(sizeStar, Array<plint, 2>(416, 255), lattice);
    createStar(sizeStar, Array<plint, 2>(475.2, 236.7), lattice);
    //Third Electrode line
    createStar(sizeStar, Array<plint, 2>(185.2, 386), lattice);
    createStar(sizeStar, Array<plint, 2>(245.7, 386), lattice);
    createStar(sizeStar, Array<plint, 2>(305.5, 376.5), lattice);
    createStar(sizeStar, Array<plint, 2>(355.6, 377.5), lattice);
    createStar(sizeStar, Array<plint, 2>(425, 386.5), lattice);
    createStar(sizeStar, Array<plint, 2>(484.2, 376.8), lattice);
    //Fourth Electrode line
    createStar(sizeStar, Array<plint, 2>(185.2, 526), lattice);
    createStar(sizeStar, Array<plint, 2>(245, 536), lattice);
    createStar(sizeStar, Array<plint, 2>(293.2, 536), lattice);
    createStar(sizeStar, Array<plint, 2>(354.3, 536), lattice);
    createStar(sizeStar, Array<plint, 2>(415, 526.5), lattice);
    createStar(sizeStar, Array<plint, 2>(485, 524), lattice);    
        //Inserting Upper Wall BounceBack 
    Box2D above_bounce_back(0, nx - 1, ny - 20, ny - 1);
    defineDynamics(lattice, above_bounce_back, new BounceBack<T,DESCRIPTOR>((T)-999999));
    //Inserting Down Wall BounceBack 
    Box2D down_bounce_back(0, nx - 1, 0 , 20);
    defineDynamics(lattice, down_bounce_back, new BounceBack<T,DESCRIPTOR>((T)-999999));
     


    // Call initialize to get the lattice ready for the simulation.
    lattice.initialize();
}

/// Produce a GIF snapshot of the velocity-norm.
void writeGif(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, plint iter)
{
    const plint imSize = 600;

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("u", iter, 6),
                               *computeVelocityNorm(lattice),
                               imSize, imSize );
}

/// Write the full velocity and the velocity-norm into a VTK file.
void writeVTK(MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX()*60;
    T dt = parameters.getDeltaT()*60;
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<2,float>(*computeVelocity(lattice), "velocity", dx/dt);
}

T computeRMSerror ( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                    IncomprFlowParam<T> const& parameters )
{
    MultiTensorField2D<T,2> analyticalVelocity(lattice);
    setToFunction( analyticalVelocity, analyticalVelocity.getBoundingBox(),
                   PoiseuilleVelocity<T>(parameters) );
    MultiTensorField2D<T,2> numericalVelocity(lattice);
    computeVelocity(lattice, numericalVelocity, lattice.getBoundingBox());

           // Divide by lattice velocity to normalize the error
    return 1./parameters.getLatticeU() *
           // Compute RMS difference between analytical and numerical solution
               std::sqrt( computeAverage( *computeNormSqr(
                              *subtract(analyticalVelocity, numericalVelocity)
                         ) ) );
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");

    IncomprFlowParam<T> parameters(
            (T) 1.758e-2,  // uMax
            (T) 113.8,    // Re
            600,        // N
            1.,        // lx
            1.0833333333333333        // ly 
    );
    const T logT     = (T)0.1;
    const T imSave   = (T)0.2;
    const T vtkSave  = (T)10.;
    const T maxT     = (T)20.1;
    // Change this variable to "pressure" if you prefer a pressure boundary
    //   condition with Poiseuille profile for the inlet and the outlet.
    const InletOutletT inletOutlet = velocity;

    writeLogFile(parameters, "Poiseuille flow");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
              parameters.getNx(), parameters.getNy(),
              new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    channelSetup(lattice, parameters, *boundaryCondition, inletOutlet);

    // Main loop over time iterations.
    for (plint iT=0; iT*parameters.getDeltaT()<maxT; ++iT) {
       if (iT%parameters.nStep(imSave)==0) {
            pcout << "Saving Gif ..." << endl;
            writeGif(lattice, iT);
        }

        //if (iT%parameters.nStep(vtkSave)==0 && iT>0) {
            pcout << "Saving VTK file ..." << endl;
            writeVTK(lattice, parameters, iT);
        //}

        if (iT%parameters.nStep(logT)==0) {
            pcout << "step " << iT
                  << "; t=" << iT*parameters.getDeltaT()
                  << "; RMS error=" << computeRMSerror(lattice, parameters);
            Array<T,2> uCenter;
            lattice.get(parameters.getNx()/2,parameters.getNy()/2).computeVelocity(uCenter);
            pcout << "; center velocity=" << uCenter[0]/parameters.getLatticeU() << endl;
        }

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();
    }

    delete boundaryCondition;
}
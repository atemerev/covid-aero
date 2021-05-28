
#include<cmath>
#include<random>
#include "olb3D.h"

#ifndef OLB_PRECOMPILED // Unless precompiled version is used

#include "olb3D.hh"     // Include full template code

#endif

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;

// Choose your turbulent model of choice
//#define RLB
#define Smagorinsky //default
//#define ConsitentStrainSmagorinsky
//#define ShearSmagorinsky
//#define Krause

#ifdef ShearSmagorinsky
typedef D3Q19<AV_SHEAR> DESCRIPTOR;
#else
typedef D3Q19<> DESCRIPTOR;
#endif

// Parameters for the simulation setup
const int N = 4;                 // resolution of the model, for RLB N>=5, others N>=2, but N>=5 recommended
const int M = 25;                 // time discretization refinement
const int inflowProfileMode = 0; // block profile (mode=0), power profile (mode=1)
const T maxPhysT = 90.;         // max. simulation time in s, SI unit
const T OUT_STEPS = 360;

const T BLOCK_SIZE = 0.5; // m
const T HEIGHT = 4; // m
const T LENGTH = 20; // m
const T INLET_RADIUS = 0.5; // m
const T INLET_LENGTH = 0.5; // m
const T MAX_VELOCITY = 1.3; // m/s

constexpr int NUM_PARTICLES = 1000;
constexpr T DENSITY = 1000.; // kg/m^3

constexpr T MEAN_RADIUS_B = 1.2e-6; // m
constexpr T DEVIATION_B = 8.4;
constexpr T MEAN_RADIUS_L = 4e-5; // m
constexpr T DEVIATION_L = 8.4;
constexpr T MEAN_RADIUS_O = 1e-5; // m
constexpr T DEVIATION_O = 8.4;

constexpr double prob_B = 0.2;
constexpr double prob_L = 0.4;

//std::lognormal_distribution<T> dist(std::log(MEAN_RADIUS), std::log(8.4));

// distributions for the BLO model

std::uniform_real_distribution<double> uni(0, 1);
std::lognormal_distribution<T> dist_B(std::log(MEAN_RADIUS_B), std::log(DEVIATION_B));
std::lognormal_distribution<T> dist_L(std::log(MEAN_RADIUS_L), std::log(DEVIATION_L));
std::lognormal_distribution<T> dist_O(std::log(MEAN_RADIUS_O), std::log(DEVIATION_O));

template<typename T, typename _DESCRIPTOR>
class TurbulentVelocity3D : public AnalyticalF3D<T, T> {

protected:
    // block profile (mode=0), power profile (mode=1)
    int _mode;
    T rho;
    T nu;
    T u0;
    T p0;
    T charL;
    T dx;


public:
    TurbulentVelocity3D(UnitConverter<T, _DESCRIPTOR> const &converter, int mode = 0) : AnalyticalF3D<T, T>(3) {
        _mode = mode;
        u0 = converter.getCharLatticeVelocity();
        rho = converter.getPhysDensity();
        nu = converter.getPhysViscosity();
        charL = converter.getCharPhysLength();
        p0 = converter.getCharPhysPressure();
        dx = converter.getConversionFactorLength();

        this->getName() = "turbulentVelocity3d";
    };

    // TODO check the turbulence operator correctness (what method is used)

    bool operator()(T output[], const BaseType<T> input[]) override {
        T y = input[1];
        T z = input[2];
        // block profile inititalization
        T u_calc = u0;
        // power profile inititalization
        if (_mode == 1) {
            T obst_y = 5.5 + dx;
            T obst_z = 5.5 + dx;
            T obst_r = 0.5;

            T B = 5.5;
            T kappa = 0.4;
            T ReTau = 183.6;

            u_calc = u0 / 7. * (2. * nu * ReTau / (charL * kappa) *
                                log(fabs(2. * ReTau / charL * (obst_r - sqrt(pow(y - obst_y, 2.)
                                                                             + pow(z - obst_z, 2.))) * 1.5 *
                                         (1 + sqrt(pow(y - obst_y, 2.)
                                                   + pow(z - obst_z, 2.)) / obst_r) /
                                         (1 + 2. * pow(sqrt(pow(y - obst_y, 2.)
                                                            + pow(z - obst_z, 2.)) / obst_r, 2.))) + B));
        }
        T a = -1., b = 1.;
        T nRandom = rand() / (T) RAND_MAX * (b - a) + a;

        output[0] = u_calc + 0.15 * u0 * nRandom;
        output[1] = 0.15 * u0 * nRandom;
        output[2] = 0.15 * u0 * nRandom;
        return true;
    };
};


void prepareGeometry(UnitConverter<T, DESCRIPTOR> const &converter, IndicatorF3D<T> &indicator,
                     SuperGeometry3D<T> &superGeometry) {

    OstreamManager clout(std::cout, "prepareGeometry");
    clout << "Prepare Geometry ..." << std::endl;

    // Sets material number for fluid and boundary
    superGeometry.rename(0, 2, indicator);

    T voxel = converter.getConversionFactorLength();

    Vector<T, 3> origin(T(),
                        HEIGHT / 2 + voxel,
                        HEIGHT / 2 + voxel);

    Vector<T, 3> extend(INLET_LENGTH,
                        HEIGHT / 2 + voxel,
                        HEIGHT / 2 + voxel);

    IndicatorCylinder3D<T> inletCylinder(extend, origin, INLET_RADIUS);
    superGeometry.rename(2, 1, inletCylinder);

    origin[0] = INLET_LENGTH;
    origin[1] = HEIGHT / 2 + voxel;
    origin[2] = HEIGHT / 2 + voxel;

    extend[0] = LENGTH;
    extend[1] = HEIGHT / 2 + voxel;
    extend[2] = HEIGHT / 2 + voxel;

    IndicatorCylinder3D<T> injectionTube(extend, origin, HEIGHT / 2);
    superGeometry.rename(2, 1, injectionTube);

    origin[0] = voxel;
    origin[1] = HEIGHT / 2 + voxel;
    origin[2] = HEIGHT / 2 + voxel;

    extend[0] = T();
    extend[1] = HEIGHT / 2 + voxel;
    extend[2] = HEIGHT / 2 + voxel;

    IndicatorCylinder3D<T> cylinderIN(extend, origin, INLET_RADIUS);
    superGeometry.rename(1, 3, cylinderIN);

    origin[0] = LENGTH - voxel;
    origin[1] = HEIGHT / 2 + voxel;
    origin[2] = HEIGHT / 2 + voxel;

    extend[0] = LENGTH;
    extend[1] = HEIGHT / 2 + voxel;
    extend[2] = HEIGHT / 2 + voxel;

    IndicatorCylinder3D<T> cylinderOUT(extend, origin, HEIGHT / 2);
    superGeometry.rename(1, 4, cylinderOUT);

    // Removes all not needed boundary voxels outside the surface
    superGeometry.clean();
    // Removes all not needed boundary voxels inside the surface
    superGeometry.innerClean();
    superGeometry.checkForErrors();

    superGeometry.print();

    clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(SuperLattice3D<T, DESCRIPTOR> &sLattice,
                    UnitConverter<T, DESCRIPTOR> const &converter,
                    Dynamics<T, DESCRIPTOR> &bulkDynamics,
                    SuperGeometry3D<T> &superGeometry) {


    OstreamManager clout(std::cout, "prepareLattice");
    clout << "Prepare Lattice ..." << std::endl;

    const T omega = converter.getLatticeRelaxationFrequency();

    // Material=0 -->do nothing
    sLattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>());

    // Material=1 -->bulk dynamics
    // Material=3 -->bulk dynamics (inflow)
    // Material=4 -->bulk dynamics (outflow)
    sLattice.defineDynamics(superGeometry.getMaterialIndicator({1, 3, 4}), &bulkDynamics);

    // Material=2 -->bounce back
    sLattice.defineDynamics(superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>());

    setInterpolatedVelocityBoundary<T, DESCRIPTOR>(sLattice, omega, superGeometry, 3);
    setInterpolatedPressureBoundary<T, DESCRIPTOR>(sLattice, omega, superGeometry, 4);

    clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(UnitConverter<T, DESCRIPTOR> const &converter,
                       SuperLattice3D<T, DESCRIPTOR> &lattice, SuperGeometry3D<T> &superGeometry, int iT) {

    OstreamManager clout(std::cout, "setBoundaryValues");

    if (iT == 0) {
        AnalyticalConst3D<T, T> rhoF(1);
        std::vector<T> velocity(3, T());
        AnalyticalConst3D<T, T> uF(velocity);

        // Seeding of random fluctuations and definition of the velocity field
        srand(time(nullptr));
        TurbulentVelocity3D<T, DESCRIPTOR> uSol(converter, inflowProfileMode);

        lattice.iniEquilibrium(superGeometry.getMaterialIndicator({1, 2, 4}), rhoF, uF);
        lattice.iniEquilibrium(superGeometry, 3, rhoF, uSol);

        lattice.defineU(superGeometry, 3, uSol);
        lattice.defineRho(superGeometry, 4, rhoF);

        // Make the lattice ready for simulation
        lattice.initialize();
    }
}

void getResults(SuperLattice3D<T, DESCRIPTOR> &sLattice,
                UnitConverter<T, DESCRIPTOR> const &converter, int iT,
                SuperGeometry3D<T> &superGeometry, Timer<T> &timer) {

    OstreamManager clout(std::cout, "getResults");
    SuperVTMwriter3D<T> vtmWriter("aero");

    if (iT == 0) {
        // Writes the geometry, cuboid no. and rank no. as vti file for visualization
        SuperLatticeGeometry3D<T, DESCRIPTOR> geometry(sLattice, superGeometry);
        SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(sLattice);
        SuperLatticeRank3D<T, DESCRIPTOR> rank(sLattice);
        vtmWriter.write(geometry);
        vtmWriter.write(cuboid);
        vtmWriter.write(rank);
        vtmWriter.createMasterFile();
    }

    // Writes the vtk files
    if (iT % converter.getLatticeTime(maxPhysT / OUT_STEPS) == 0) {
        // Create the data-reading functors...
        SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
        SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
        vtmWriter.addFunctor(velocity);
        vtmWriter.addFunctor(pressure);
        vtmWriter.write(iT);

        SuperEuklidNorm3D<T, DESCRIPTOR> normVel(velocity);
        BlockReduction3D2D<T> planeReduction(normVel, Vector<T, 3>({0, 1, 0}));
        // write output as JPEG
        heatmap::write(planeReduction, iT);
    }

    // Writes output on the console
    if (iT % converter.getLatticeTime(maxPhysT / OUT_STEPS) == 0) {
        timer.update(iT);
        timer.printStep();
        sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    }
}


int main(int argc, char *argv[]) {

    // === Random generators, etc.

    std::mt19937 generator(200);

    // === 1st Step: Initialization ===

    olbInit(&argc, &argv);
    singleton::directories().setOutputDir("./tmp/");
    OstreamManager clout(std::cout, "main");
    // display messages from every single mpi process
    // clout.setMultiOutput(true);

    UnitConverter<T, DESCRIPTOR> const converter(
            (T) BLOCK_SIZE / N,
            (T) BLOCK_SIZE / (M * N),
            (T) BLOCK_SIZE,
            (T) MAX_VELOCITY,
            (T) 1.5e-5,
            (T) 1.225
            );

    // Prints the converter log as console output
    converter.print();
    // Writes the converter log in a file
    converter.write("aero");

    T voxel = converter.getConversionFactorLength();

    Vector<T, 3> origin;
    Vector<T, 3> extend(LENGTH,
                        HEIGHT + 2. * voxel,
                        HEIGHT + 2. * voxel);

    IndicatorCuboid3D<T> cuboid(extend, origin);

    CuboidGeometry3D<T> cuboidGeometry(cuboid, converter.getConversionFactorLength(), singleton::mpi().getSize());
    HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

    // === 2nd Step: Prepare Geometry ===

    SuperGeometry3D<T> superGeometry(cuboidGeometry, loadBalancer, 2);
    prepareGeometry(converter, cuboid, superGeometry);

    // === 3rd Step: Prepare Lattice ===

    SuperLattice3D<T, DESCRIPTOR> sLattice(superGeometry);

    Dynamics<T, DESCRIPTOR> *bulkDynamics;
    const T omega = converter.getLatticeRelaxationFrequency();
#if defined(RLB)
    bulkDynamics = new RLBdynamics<T, DESCRIPTOR>( omega, instances::getBulkMomenta<T, DESCRIPTOR>() );
#elif defined(Smagorinsky)
    bulkDynamics = new SmagorinskyBGKdynamics<T, DESCRIPTOR>(omega, instances::getBulkMomenta<T, DESCRIPTOR>(),
                                                             0.15);
#elif defined(ShearSmagorinsky)
    bulkDynamics = new ShearSmagorinskyBGKdynamics<T, DESCRIPTOR>( omega, instances::getBulkMomenta<T, DESCRIPTOR>(),
        0.15);
#elif defined(Krause)
    bulkDynamics = new KrauseBGKdynamics<T, DESCRIPTOR>( omega, instances::getBulkMomenta<T, DESCRIPTOR>(),
        0.15);
#else //ConsitentStrainSmagorinsky
    bulkDynamics = new ConStrainSmagorinskyBGKdynamics<T, DESCRIPTOR>( omega, instances::getBulkMomenta<T, DESCRIPTOR>(),
        0.05);
#endif


    prepareLattice(sLattice, converter, *bulkDynamics, superGeometry);

    clout << "particle system within geometry" << std::endl;

    SuperParticleSystem3D<T, Particle3D> superParticleSystem(superGeometry);
    SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR> velocityProbe(sLattice, converter);
    clout << "Proceeding with particle simulation..." << std::endl;

    clout << "simulation timer" << std::endl;
    Timer<T> particleTimer(converter.getLatticeTime(maxPhysT), NUM_PARTICLES);

    clout << "attaching VTK writers" << std::endl;
    SuperParticleSysVtuWriter<T, Particle3D> particleVtuWriter(superParticleSystem, "particles",
                                                               SuperParticleSysVtuWriter<T, Particle3D>::particleProperties::velocity
                                                               |
                                                               SuperParticleSysVtuWriter<T, Particle3D>::particleProperties::mass
                                                               |
                                                               SuperParticleSysVtuWriter<T, Particle3D>::particleProperties::radius);

    clout << "generating particles" << std::endl;
    // origin: center of the cuboid, length 5 cm; radius 1 cm
    Vector<T, 3> originInlet(T(), HEIGHT / 2 + voxel,HEIGHT / 2 + voxel);
    Vector<T, 3> extendInlet(INLET_LENGTH, HEIGHT / 2 + voxel, HEIGHT / 2 + voxel);

    IndicatorCylinder3D<T> inletCylinder(extendInlet, originInlet, INLET_RADIUS);

    for (int i = 0; i < NUM_PARTICLES; i++) {

        double uni_random = uni(generator);
        // choosing the B/L/O distribution from the BLO model modal probabilities
        std::lognormal_distribution<T> dist = uni_random < prob_B ? dist_B : uni_random < prob_L ? dist_L : dist_O;
        T radius = dist(generator);
        clout << "Particle: " << radius << std::endl;
        T mass = 4. / 3. * M_PI * std::pow(radius, 3) * DENSITY;
        superParticleSystem.addParticle(inletCylinder, mass, radius);
    }

    clout << "adding gravity" << std::endl;

    std::vector<T> direction = {0, -1, 0};
    auto weightForce = make_shared
            < WeightForce3D<T, Particle3D>
            > ( direction, 9.81 );
    superParticleSystem.addForce(weightForce);

    clout << "adding Stokes drag force" << std::endl;

    auto stokesDragForce = make_shared
            < StokesDragForce3D<T, Particle3D, DESCRIPTOR>
            > ( velocityProbe, converter );
    superParticleSystem.addForce(stokesDragForce);
    superParticleSystem.setVelToFluidVel(velocityProbe);


    clout << "setting up boundary material (settle on outflow [4] and borders[2] )" << std::endl;
    std::set<int> boundMaterial = {2, 4};
    auto materialBoundary = make_shared
            < MaterialBoundary3D<T, Particle3D>
            > (superGeometry, boundMaterial);
    superParticleSystem.addBoundary(materialBoundary);

    // === 4th Step: Main Loop with Timer ===

    Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel());
    timer.start();

    for (std::size_t iT = 0; iT <= converter.getLatticeTime(maxPhysT); ++iT) {
        // === 5ath Step: Apply filter
#ifdef ADM
        SuperLatticeADM3D<T, DESCRIPTOR> admF( sLattice, 0.01, 2 );
        admF.execute( superGeometry, 1 );
#endif
        // === 5bth Step: Definition of Initial and Boundary Conditions ===
        setBoundaryValues(converter, sLattice, superGeometry, iT);

        // === 6th Step: Collide and Stream Execution ===
        sLattice.collideAndStream();

        // === 7th Step: Computation and Output of the Results ===
        getResults(sLattice, converter, iT, superGeometry, timer);
    }
    timer.stop();
    timer.printSummary();
//    delete bulkDynamics;

    clout << "go!" << std::endl;
    particleTimer.start();
    clout << "timer started" << std::endl;
    for (std::size_t iT = 0; iT <= converter.getLatticeTime(maxPhysT); ++iT) {
        superParticleSystem.simulate(converter.getConversionFactorTime());
        if (iT % converter.getLatticeTime(maxPhysT / OUT_STEPS) == 0) {
            particleTimer.update(iT);
            particleTimer.printStep();
            superParticleSystem.print();
            particleVtuWriter.write(iT);
        }
    }
    particleTimer.stop();
    particleTimer.printSummary();
}

#pragma clang diagnostic pop